module preprocessing

  integer tmin, tmax, tstp, maxn
  double precision coff
  logical preprocess

  integer numg
  integer, allocatable :: bcind(:,:), c0ci(:,:), glist(:,:)
  double complex, allocatable :: bc(:)
  private :: numg, bcind, c0ci, glist, bc

 contains

  subroutine preprocess_datacube

    use lattice, only : nlx, nly, nlz, nx, ny, nz, ax0, ay0, az0, nslice
    use wf, only : r3wf
    use param, only : bohrtonm, pi

    implicit none

    include 'fftw3.f'

    integer*8 :: plan
    double complex, allocatable :: incut(:,:,:), out(:,:,:)
    double precision bcre, bcim, bcabs

    double precision dx, dy, dz, ox, oy, oz, latx, laty, latz
    integer nat, nzd, ncelx, ncely, ncelz
    integer ibc, nbc
    integer ix, iy, iz, i, j, k, iret
    character*12 fname

    double complex                 bcval
    double precision               recip, recx, recy, recz
    integer                        ig, ncells
    integer                        hmax, kmax, lmax, hh, kk, ll, igam
    logical                        newg

    ncelx = nlx
    ncely = nly
    ncelz = nlz
    latx = ax0*bohrtonm*dble(10)
    recx = dble(2)*pi/latx
    laty = ay0*bohrtonm*dble(10)
    recy = dble(2)*pi/laty
    latz = az0*bohrtonm*dble(10)
    recz = dble(2)*pi/latz

    write(*,'(a,3f5.2,a)') ' Lattice parameters a,b,c: ',latx,laty,latz, ' Angstrom'
    write(*,'(a,3i5)') ' Size of grid   (nx,ny,nz):', nx, ny, ncelz*nslice
    write(*,'(a,3i5)') ' Cells (ncelx,ncely,ncelz):', ncelx, ncely, ncelz
    nzd = nslice
    write(*,'(i4,a)') nzd, ' slices per cell in z-direction'


  !!! here
    print *, 'Starting cycle over thickness'
    allocate(bc(maxn),bcind(maxn,3),c0ci(maxn,2),glist(maxn,3))

    call dfftw_init_threads(iret)
    print *, iret
    call dfftw_plan_with_nthreads(8)

    do k = tmin, tmax, tstp

      ncells = k/nzd

      allocate(incut(1:k,ny,nx),out(1:k,ny,nx))
      do ix=1, nx
        do iy=1, ny
          incut(1:k,iy,ix) = r3wf(ix-1,iy-1,0:k-1)
        enddo
      enddo

      print *, 'FFT...'
      call dfftw_plan_dft_3d(plan,k,ny,nx,incut,out,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,incut,out)
      call dfftw_destroy_plan(plan)
      out = out / sqrt(dble(nx*ny)) / dble(k)

      print *, 'Filtering...'
      ibc = 0
      do ix = 1, nx
        do iy = 1, ny
          do iz = 1, k
            bcval = out(iz,iy,ix)
            if(abs(bcval)>coff) then
              ibc = ibc + 1
  !            print *, ibc, bcval
              if(ibc>maxn) stop 'ERROR: enlarge the maxn'
              bc(ibc) = bcval
              bcind(ibc,1) = ix-1
              bcind(ibc,2) = iy-1
              bcind(ibc,3) = iz-1
            endif
          enddo
        enddo
      enddo
      nbc = ibc

      call quicksort(1,nbc)

  !   construct list of G-vectors and indices
      numg = 0
      do ibc = 1, nbc
        hh = bcind(ibc,1)
        if(hh>(nx/2)) hh = hh - nx
        kk = bcind(ibc,2)
        if(kk>(ny/2)) kk = kk - ny
        call gz2ggam(bcind(ibc,3),k,ncells,ll,igam)
!        print *, ig, h, k, l
        newg = isknowng(hh,kk,ll,ig)
        c0ci(ibc,1) = igam
        c0ci(ibc,2) = ig
      enddo

      if(ncells<10) then
        write(fname,'(a7,i1,a1)') 'thick__',ncells,'a'
      else if(ncells<100) then
        write(fname,'(a6,i2,a1)') 'thick_',ncells,'a'
      else
        write(fname,'(a5,i3,a1)') 'thick',ncells,'a'
      endif

      write(fname,'(f6.2)') dble(ncells)*latz/dble(10)

      open(74,file=fname)
      print *, 'Writing ', fname
      write(74,'(g10.3,a)') coff, ' # cut-off parameter'
      write(74,'(2i3,i4,a)') ncelx, ncely, ncells, ' # number of cells in x, y, z-direction; z defines gammas as igam/ncelz'
      write(74,'(i10,a)')   numg, ' # number of reciprocal lattice vectors G'
      write(74,'(i10,a)')   nbc,  ' # number of F_\tilde{g} larger than cut-off'
      write(74,*)
      write(74,'(a)') '# list of G-vectors: h, k, l and values in atomic units. WARNING: h & k are recips of the ncelx x ncely supercell'
      do ig = 1, numg ! list of G-vectors
        hh = glist(ig,1)
        kk = glist(ig,2)
        ll = glist(ig,3)
        write(74,'(3i6,3f14.8)') hh, kk, ll, recx*dble(hh)/dble(ncelx), recy*dble(kk)/dble(ncely), recz*dble(ll)
      enddo ! list of Fg's
      write(74,*)
      write(74,'(a)') '# list of igam, iG, F_\tilde{g}'
      do ibc=1, nbc
        write(74,'(2i10,2e18.10,3x,e12.4)') c0ci(ibc,1), c0ci(ibc,2), bc(ibc), abs(bc(ibc))
      enddo
      close(74)
      print *, 'Done!'

      deallocate(incut,out)

    enddo

  end subroutine preprocess_datacube


! for given "gz" (and parameters gmax, ncells) returns the "l" (from hkl) and index of gamma
! the real gamma will be then dble(igam/ncells)
  subroutine gz2ggam(gz,gmax,ncells,l,igam)

    implicit none

    integer gz, gmax, ncells, l, igam
    integer gztmp

    gztmp = gz
    if(gz>(gmax/2)) gztmp = gztmp - gmax
    l = nint(gztmp/dble(ncells))
    igam = gztmp - l*ncells

  end subroutine gz2ggam


! checks, if we already registered this G-vector
! if not yet, it returns true and adds it to the end of the list and increases numg (which is its index in array glist)
! if yes, then it returns false and sets ig to a an index of this vector in the list
  logical function isknowng(h,k,l,ig)

    implicit none

    integer h, k, l, ig
    integer i
    logical match

    match = .false.

    do i = 1, numg
      if(h==glist(i,1) .and. k==glist(i,2) .and. l==glist(i,3)) match = .true.
      if(match) exit
    enddo

    if(.not.match .and. (i>numg)) then
      numg = numg + 1
      glist(numg,1) = h
      glist(numg,2) = k
      glist(numg,3) = l
      ig = numg
    endif

    if(match) ig = i

    isknowng = match

  end function isknowng


  subroutine swapbc(i,j)

    implicit none

    integer i,j
    double complex tmpc
    integer ind(3)

    tmpc  = bc(i)
    bc(i) = bc(j)
    bc(j) = tmpc

    ind        = bcind(i,:)
    bcind(i,:) = bcind(j,:)
    bcind(j,:) = ind

  end subroutine swapbc


  integer function comparebc(i,j)

    implicit none

    integer i,j
    double precision bci, bcj

    bci = abs(bc(i))
    bcj = abs(bc(j))

    if(bci>bcj) then
      comparebc =  1
    else if(bci<bcj) then
      comparebc = -1
    else
      comparebc =  0
    endif

  end function comparebc


  recursive subroutine quicksort(first,last)

    implicit none

    integer first, last

    integer ipivot, istore, ind
    real rnd

    call random_number(rnd)
    ipivot = nint(first+rnd*(last-first))

    if(ipivot/=last) call swapbc(ipivot,last)

    istore = first
    do ind=first, last-1
      if(comparebc(ind,last)>0) then
        call swapbc(ind,istore)
        istore = istore + 1
      endif
    enddo

    call swapbc(istore,last)

    if(first<istore-1) call quicksort(first,istore-1)
    if(istore+1<last)  call quicksort(istore+1,last)

  end subroutine quicksort

end module preprocessing
