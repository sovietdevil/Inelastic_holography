subroutine sort_blochs

  use defs
  use beams
  use blochs
  use inputs
  use multislice
  use reallocate
  use struct, only: aa, bb, cc, nat, pos, nvbr, br1_dir, vbr, mult
  use hashtable
  
  implicit none
  
  integer nslic
  
  integer i, j, ig, ih, ig0, ih0
  double precision qfac, qlenfactor
  
  double complex, allocatable :: dzer(:)
  
  integer, allocatable :: nccst(:), nccln(:)
  integer ist, imh, img
  
  integer nat2, iat2                   ! nat2 -> number of atoms in unit cell, where excitation can take place -> nvbr*mult(atom)
  double precision sumall(9)           ! total scattering cross-section -> an array of up to 9 elements (nx, ny, nz, sx, sy, sz, ox, oy, oz)
  double complex sumi(3), cprod        ! term of the sum
  double precision, allocatable :: suma(:,:,:) ! scattering cross-section per atom (1st dim 9, 2nd dim -> # ats. in unit cell, 3rd dime -> ncelz)
  double precision, allocatable :: apos(:,:,:) ! list ot all atom positions (1st dim -> x,z,y, 2nd dim -> # ats. in unit cell, 3rd dim -> ncelz)
  double precision tmppos(3)           ! temporary atom position
  
  integer snum, snmax, onmax           ! contains number of s-vectors and the allocated array size
  integer, pointer :: shash(:,:,:)     ! structure for a quick look-up for new s-vectors
  integer, pointer :: svec(:,:)        ! list of s-vectors indexed by "is" (2nd index), hkl is the first index
  logical is_new_s,is_new2             ! used when generating svec
  integer is, s(3), shi, shj, shk      ! index of s and a particular s-vector and its hash indices
  integer dgam                         ! represents discretized and mod-ed [gamma^(l)-gamma^(j)]*cc*ncelz
  integer iX, nX, nXmax, nXa, nXtot    ! index and a number of X-terms, nXmax is allocated array size
  integer, pointer :: Xi(:,:)          ! maps "ix" (1st index) to a pair of "k,is" (2nd index; 1->k, 2->is)
!  integer, pointer :: Xb(:,:)          ! maps back the pair of indices "k,is" to "ix"
  double complex, pointer :: X(:)      ! contains the X-term list indexed by "ix"
  double complex Xterm                 ! particular X-term
  
  double precision, allocatable :: qvec(:,:), qlen2(:) ! list of q-vectors (actually, rather h-g+(gam^l-gam^j)n, without qk=kf-ki) and |q+qk|^2
  double precision, allocatable :: qvec2(:,:)
  double precision q(3), qk(3)   ! single q-vector, qk = kf-ki
  double precision isqrt2        ! one over square root of 2
  double precision maxgam        ! for finding maximal absolute value of gamma^(l)
  integer          imaxgam       ! the same in multiples of c* = cc/2pi
  
  type(cell), pointer :: acell
  
  external compareblochs_in, swapblochs_in
  external compareblochs_out, swapblochs_out
  external qlenfactor
  
  isqrt2 = 1.d0/sqrt(2.d0)
  
  nslic = 1000*ncelz
  
  allocate(d0d(nbeams_out*nbeams_out))
  allocate(d0di(2,nbeams_out*nbeams_out))
  allocate(dzer(nbeams_out))
  nd0d=0
  
  allocate(nccln(ncelx*ncely),nccst(ncelx*ncely))
  nccln(:) = 0
  nccst(:) = 0

!do i=0, ncelz
!  print *, i, gam_in(i)
!enddo
!stop

!THICKNESS IS FULLY DETERMINED BY INPUT FROM MULTISLICE VIA NSCZ
!* we should ignore the "tmin, tstep, tmax", they are obsolete
!* nscz then also determines the discretization grid of gammas below, perhaps best as 0..(nscz-1)?
  nat2 = nvbr*mult(atom)
  sumall = 0.d0
  allocate(suma(9,nat2,ncelz))
  suma(:,:,:) = dble(0)
  allocate(apos(3,nat2,ncelz))
  j = 0
  if(atom>1) then
    do i=1, atom-1
      j = j + mult(i)
    enddo
  endif
  iat2 = 0
  do i=1, mult(atom)
    j = j+1
    do ig=1, nvbr
      tmppos = pos(:,j)+vbr(:,ig)
      do ih=1,3
        tmppos(ih) = tmppos(ih) - dble(floor(tmppos(ig)))
      enddo
      iat2 = iat2+1
      apos(:,iat2,1) = matmul(br1_dir,tmppos)
    enddo
  enddo
  if(ncelz>1) then
    do i=2, ncelz
      do iat2=1, nat2
        apos(:,iat2,i) = apos(:,iat2,1)
        apos(3,iat2,i) = apos(3,iat2,i) + cc*dble(i-1)
!write(*,'(2i4,3f14.6)') i, iat2, apos(:,iat2,i)
      enddo
    enddo
  endif
  write(*,'(a,i6)') 'Atom list built. Total number of atoms: ', nat2*ncelz

! 1) preprocess the C0C list:
! * make nscx * nscy array, where all C0C from multislice are partitioned (a single list with indexing)
! * sort sublists one-by-one by size of C0C
! sort it by mod(h,nscx), then mod(k,nscy), and then by size
  call quicksort(1,nc0c,compareblochs_in,swapblochs_in)
  do i=1, nc0c
!print *, abs(c0c(i)), hklbeams_in(1:3,c0ci(1,i))
    img = modulo(hklbeams_in(1,c0ci(1,i)),ncelx)
    imh = modulo(hklbeams_in(2,c0ci(1,i)),ncely)
    nccln(img*ncely+imh+1) = nccln(img*ncely+imh+1) + 1
  enddo
  ist = 1
  do i=1, ncelx*ncely
    nccst(i) = ist
    ist = ist+nccln(i)
!    print *, i, nccst(i), nccln(i)
  enddo

! 2) preprocess the D0D list:
! * sort the list by size of D0D
! We multiply all D(j)_G by D(j)_0* and store products bigger than c0_min
  dzer(:) = mat_out(zerohkl_out,:)
  do i=1, nbeams_out
    do j=1, nbeams_out
      mat_out(i,j) = mat_out(i,j) * dzer(j)
!      if(dble(10)*abs(mat_out(i,j))>(c0_min/(maxfg*maxfg*maxd0d))) then
      if(abs(mat_out(i,j))>(c0_min/maxfg)) then
!write(*,'(2i4,2f14.6,4x,2(2f14.6,2x))') i, j, dzer(j), mat_out(i,j)
        nd0d = nd0d + 1
        d0d(nd0d) = mat_out(i,j)
        d0di(1,nd0d) = i
        d0di(2,nd0d) = j
      endif
    enddo
  enddo
  call quicksort(1,nd0d,compareblochs_out,swapblochs_out)
!do i=1, nd0d
!print *, abs(d0d(i)), hklbeams_out(1:3,d0di(1,i))
!enddo

  write(*,'(a,2i8)') 'Done sorting C0Cg and D0Dh terms', nc0c, nd0d

! 3) construct the ccdd lists (there are nscx * nscy of them):
! * for each of the sublists:
!   - if |C0C*D0D|>Coff we map g,h,gamma^(j),gamma^(l) to a pair s,dgam(k) such that dgam(k)=gamma^(l)-gamma^(j)=k*(c*)/N and s=h-g;
!     if |dgam(k)|>c*/2, then we update s with eventual +/- (0,0,c*) to make |dgam(k)|<=c*/2
!   - keep track of different "s"-vectors (for each sublist, probably best using integer hkl supercell indices)
!   - q(is,k) = s(is) + (0,0,dgam(k)), where "is" is integer index of "s" vector from the list defined just one line above
!   - add the C0C*D0D*exp(i*gamma^(l)*t) term to the complex array X(ix)

  nxtot = 0
  nxmax = nc0c+nd0d
  snmax = 2*(nbeams_in+nbeams_out)
  allocate(svec(3,snmax),X(nxmax),Xi(nxmax,2)) !,Xb(2*nslic,snmax))

! find maximal gam_out so that we finally fix the segmentation fault in accesses to shash
  maxgam=0.d0
  do j=1, nd0d
    if(abs(gam_out(d0di(2,j)))>abs(maxgam)) maxgam=gam_out(d0di(2,j))
  enddo
  imaxgam = nint(abs(maxgam*cc/(2.d0*pi)))+1
  allocate(shash( -hmax-meshx:hmax+meshx, -kmax-meshy:kmax+meshy, -lmax-ms_slices-imaxgam:lmax+ms_slices+imaxgam))


  do ist=1, ncelx*ncely ! cycle over subgroups
!write(*,'(i4)',advance='no') ist

   if(nccln(ist)>0) then
    snum = 0
    shash = 0
    !Xb = 0
    call hash_init(4096)
    nx = 0
    nxa = 0
    do i=nccst(ist), nccst(ist)+nccln(ist)-1 ! subgroup "ist"
      do j=1, nd0d
        Xterm = c0c(i)*conjg(d0d(j))*exp(imag*gam_out(d0di(2,j))*dble(ncelz)*cc) ! <<<<<<<<< CHECK the sign of exponential >>>>>>>>>
        if (abs(Xterm)<c0_min) exit ! Xterm too small, take next i
        nxa = nxa + 1
        call ghjl2sk(hklbeams_in(1,c0ci(1,i)), hklbeams_out(1,d0di(1,j)), gam_in(c0ci(2,i)), gam_out(d0di(2,j)), s, dgam, nslic)
!write(*,'(a,6i4,2f8.4,3i4,2x,i4)') 'hkl_in/out, gam_in/out, s, dgam:', hklbeams_in(:,c0ci(1,i)), hklbeams_out(:,d0di(1,j)), gam_in(c0ci(2,i)), gam_out(d0di(2,j)), s, dgam
        ! check if that's a new s-vector
        shi = nint(dble(s(1))/ncelx)
        shj = nint(dble(s(2))/ncely)
        shk = s(3)
!write(*,'(6i5)') shi, shj, shk, s(1:3)
!write(*,'(a,4f10.6,2i4)') 'gam_in/out, gam_diff, clat, dgam, ncelz:', gam_in(c0ci(2,i)), gam_out(d0di(2,j)), gam_out(d0di(2,j))-gam_in(c0ci(2,i)), cc, dgam, ncelz
        is = shash(shi,shj,shk)
        if(is==0) then
          is_new_s = .true.
        else
          is_new_s = .false.
        endif
        if(is_new_s) then ! new s-vector
          snum = snum+1
          if(snum>snmax) then
            onmax = snmax
            snmax = nint(1.2d0*snmax)
            call doreallocate(svec,3,snmax)
          endif
          is = snum
          svec(:,is) = s
          shash(shi,shj,shk) = is
          nx = nx+1
          if(nx>nxmax) then
            nxmax = nint(1.2d0*nxmax)
            call doreallocate(X,nxmax)
            call doreallocate(Xi,nxmax,2)
          endif
          ix = nx
          X(ix) = Xterm
          Xi(ix,1) = dgam
          Xi(ix,2) = is
          acell => hash_insert(is*nslic+dgam)
          acell%value = ix
        else ! s-vector already is in the list
          acell => hash_lookup(is*nslic+dgam)
          if(associated(acell)) then
            ix = acell%value
          else
            ix = 0
          endif
          if(ix==0) then ! ... but not for such dgam
            nx = nx+1
            if(nx>nxmax) then
              nxmax = nint(1.2d0*nxmax)
              call doreallocate(X,nxmax)
              call doreallocate(Xi,nxmax,2)
            endif
            ix = nx
            X(ix) = Xterm
            Xi(ix,1) = dgam
            Xi(ix,2) = is
            acell => hash_insert(is*nslic+dgam)
            acell%value = ix
          else ! s,dgam is already known - we just add the X-term
            X(ix) = X(ix) + Xterm
          endif
        endif
!write(*,'(a,3i4,2x,2(2f10.6,2x),i4,i4,2x,4(2f10.6,2x))') 'ix,is,dgam, X,Xterm, svec, i,j, c0c*d0d, c0c, d0d, exp(i*g*t)', ix, is, dgam,  X(ix), Xterm,  i, j, c0c(i)*conjg(d0d(j)), c0c(i), d0d(j), exp(imag*gam_out(d0di(2,j))*dble(ncelz)*cc)
      enddo ! next d0d
    enddo ! next c0c
    call hash_destroy
    write(*,'(a,i5,a,i10,a,i8,a,i12)') 'Group #:',ist,', X-term groups #:',nx,', s-vectors #:',snum,', X-terms #:', nxa
!print *, ist, nx, snum, nxa

    ! generate list of q-vectors
    allocate(qvec(3,nx),qlen2(nx),qvec2(3,nx))
    qk = vchi_out - vchi_in
!print *, "QK:", qk
    do ix=1, nx
      q(1) = svec(1,Xi(ix,2))*dble(2)*pi/aa/dble(ncelx)
      q(2) = svec(2,Xi(ix,2))*dble(2)*pi/bb/dble(ncely)
      q(3) = svec(3,Xi(ix,2))*dble(2)*pi/cc
      q(3) = q(3) + dble(Xi(ix,1))*dble(2)*pi/(cc*dble(nslic))
      qvec(:,ix) = q + qk
! JR: z-locality approximation
      qvec(3,ix) = qk(3) ! full phase factor, but modified q_z in MDFF
      qlen2(ix) = qvec(1,ix)*qvec(1,ix) + qvec(2,ix)*qvec(2,ix) + qvec(3,ix)*qvec(3,ix)
      qvec2(:,ix) = q
! -----------------------------
!write(*,'(a,4i5,4f14.6,2x,f10.6)') 'ix, svec, qvec, q^2, gam_diff:', ix, svec(1:3,Xi(ix,2)), qvec(1:3,ix), qlen2(ix), dble(Xi(ix,1))*dble(2)*pi/(cc*dble(nslic))
    enddo

    ! summation for the subgroup follows
    do iat2=1, nat2
      do i=1, ncelz
        ! Nx, Ny, Nz
        sumi(:) = 0.d0
        do ix=1, nx
!          q = qvec(:,ix)
!          write(*,'(3f14.6)') q
          sumi(:) = sumi(:) + X(ix)*exp(-imag*dot_product(qvec2(:,ix),apos(:,iat2,i))) * qvec(:,ix)/qlen2(ix)
!write(*,'(i4,4f10.6)') i, qvec2(1,ix)*apos(1,iat2,i)/(2.d0*pi), qvec2(2,ix)*apos(2,iat2,i)/(2.d0*pi), qvec2(3,ix)*apos(3,iat2,i)/(2.d0*pi)
!write(*,'(3i4,3f10.4,2x,2f10.4,2x,4f10.4,2x,2f10.4,2x,6f10.4)') iat2, i, ix, apos(1:3,iat2,i), X(ix), qvec(1:3,ix), qlen2(ix), exp(-imag*dot_product(qvec(:,ix),apos(:,iat2,i))), sumi(1:3)
        enddo
!stop
!write(*,'(6f10.6)') sumi(1:3)
        suma(1,iat2,i) = suma(1,iat2,i) + real( sumi(1)*conjg(sumi(1)) )
        suma(2,iat2,i) = suma(2,iat2,i) + real( sumi(2)*conjg(sumi(2)) )
        suma(3,iat2,i) = suma(3,iat2,i) + real( sumi(3)*conjg(sumi(3)) )
        ! off-diagonal terms
        suma(4,iat2,i) = suma(4,iat2,i) + 2.d0*real( sumi(3)*conjg(sumi(2)) )
        suma(5,iat2,i) = suma(5,iat2,i) + 2.d0*real( sumi(1)*conjg(sumi(3)) )
        suma(6,iat2,i) = suma(6,iat2,i) + 2.d0*real( sumi(2)*conjg(sumi(1)) )
        ! magnetic termssuma(1,iat2,i) + 
        suma(7,iat2,i) = suma(7,iat2,i) + 2.d0*aimag( sumi(3)*conjg(sumi(2)) )
        suma(8,iat2,i) = suma(8,iat2,i) + 2.d0*aimag( sumi(1)*conjg(sumi(3)) )
        suma(9,iat2,i) = suma(9,iat2,i) + 2.d0*aimag( sumi(2)*conjg(sumi(1)) )
!write(*,'(2i4,3f10.4,2x,9f10.4)') iat2, i, apos(:,iat2,i), suma(:,iat2,i)
!write(*,'(9f10.4)') suma(:,iat2,i)
      enddo
    enddo

!    deallocate(X,Xi,svec,qvec,qlen2,qvec2,shash)
    deallocate(qvec,qlen2,qvec2)

   endif
   nxtot = nxtot + nxa
  enddo ! next subgroup
  deallocate(c0c,d0d,c0ci,d0di,dzer)
  deallocate(X,Xi,shash) ! Xb
  write(*,'(a,i16)') 'Total number of terms summed: ', nxtot
  
  ! now comes output
  sumall(:) = 0.d0
  do iat2=1,nat2
    do i=1, ncelz
      sumall = sumall + suma(:,iat2,i)
      write(81,'(2i5,3f14.6,4x,9f20.12)') iat2, i, apos(:,iat2,i)*auinnm, suma(:,iat2,i)
    enddo
  enddo
  
  write(82,'(a)') '# combined mult-BW code, thickness, xx, yy, zz, Re[yz], Re[zx], Re[xy], Im[yz]=s_x, Im[zx]=s_y, Im[xy]=s_z'
  write(82,'(f14.6,9e20.12)') cc*auinnm*dble(ncelz), sumall(:)/dble(nat2*ncelz)

end


subroutine ghjl2sk(v1,v2,g1,g2,s,k,nsl)

  use multislice
  use struct, only: cc
  use defs, only: pi

  implicit none
  
  integer v1(3),v2(3),s(3),k,l,m,nsl
  double precision g1,g2,dg
  
  s(1) = ncelx*v2(1) - v1(1)
  s(2) = ncely*v2(2) - v1(2)
  s(3) = v2(3) - v1(3)
  dg = cc*(g2-g1)/(2.d0*pi)
  m = nint(dg)
  s(3) = s(3)+m
  k = nint((dg-m)*nsl)

end subroutine ghjl2sk



double precision function qlenfactor(i,j)
! returns (k_f-k_i)/q for q given indices to beam arrays

  use inputs
  use beams

  implicit none
  
  integer i, j
  
  double precision q0(3), q(3), cart_length2
  
  external cart_length2
  
  q0(:) = vchi_out(:) - vchi_in(:)
  q(:)  = aubeams_out(:,j) - aubeams_in(:,i) + q0(:)
  
  qlenfactor = cart_length2(q0)/cart_length2(q)

end


integer function compareblochs_in(i,j)

  use blochs
  use beams, only: hklbeams_in
  use multislice

  implicit none

  integer i, j, imh, jmh, imk, jmk

  imh = modulo(hklbeams_in(1,c0ci(1,i)),ncelx)
  jmh = modulo(hklbeams_in(1,c0ci(1,j)),ncelx)
  imk = modulo(hklbeams_in(2,c0ci(1,i)),ncely)
  jmk = modulo(hklbeams_in(2,c0ci(1,j)),ncely)

  if (imh<jmh) then
    compareblochs_in = 1
  else if (imh>jmh) then
    compareblochs_in = -1
  else if (imk<jmk) then
    compareblochs_in = 1
  else if (imk>jmk) then
    compareblochs_in = -1
  else if (abs(c0c(i))>abs(c0c(j))) then
    compareblochs_in = 1
  else if (abs(c0c(i))<abs(c0c(j))) then
    compareblochs_in = -1
  else
    compareblochs_in = 0
  endif

end

integer function compareblochs_out(i,j)

  use blochs

  implicit none

  integer i, j

  if (abs(d0d(i))>abs(d0d(j))) then
    compareblochs_out = 1
  else if (abs(d0d(i))<abs(d0d(j))) then 
    compareblochs_out = -1
  else
    compareblochs_out = 0
  endif

end

subroutine swapblochs_in(i,j)

  use blochs

  implicit none

  integer i, j, k(2)
  double complex tmp

  tmp = c0c(i)
  c0c(i) = c0c(j)
  c0c(j) = tmp

  k(:) = c0ci(:,i)
  c0ci(:,i) = c0ci(:,j)
  c0ci(:,j) = k(:)

end

subroutine swapblochs_out(i,j)

  use blochs

  implicit none

  integer i, j, k(2)
  double complex tmp

  tmp = d0d(i)
  d0d(i) = d0d(j)
  d0d(j) = tmp

  k(:) = d0di(:,i)
  d0di(:,i) = d0di(:,j)
  d0di(:,j) = k(:)

end
