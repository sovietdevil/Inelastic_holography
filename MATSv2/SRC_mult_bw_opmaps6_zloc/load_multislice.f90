subroutine load_multislice

  use defs, only: imag, pi, two, auinnm
  use beams, only: hklbeams_in, aubeams_in, zerohkl_in, nbeams_in, c0_min
  use blochs, only: nc0c, c0c, c0ci, gam_in, maxfg, maxd0d
  use struct, only: cc
  use multislice

  implicit none

  integer nfg, ifg               ! number of F_g's and cycle variable
  integer ix, iy, iz, ig         ! integers of g-vector coordinates
  double precision fgre, fgim    ! re/im values of F_g
  double complex fg, transf      ! complex value of F_g
  double precision gx, gy, gz    ! to hold the g-vector coordinates
  double precision fgcoff        ! cut-off used in preprocessing of fg's
  integer ngam, lgam, hgam, igam ! to iterate gamma array

! since the lattice is often shifted in mult calculation, we need to compensate for it by use of the shift theorem for discrete Fourier transform
! here we load the shifts
  read(75,*) shiftx, shifty
  read(75,*) meshx, meshy
  write(*,'(a,2i5)') 'Shift of potential by  ', shiftx, shifty
  write(*,'(a,2i5)') 'Mesh within unit cell  ', meshx, meshy

! load the F_g coefficients and associated g-vectors
! partition them to g = G + gamma*(0,0,1), where G is reciprocal lattice vector and gamma is the smallest possible
! conveniently, gamma(j) = (j-50)*(0,0,1)/100
! store the result as the C_G^0 * C_G^j with associated gamma(j)

  read(74,*) fgcoff
! add some check, if it is not too big compared to coff of dyndif, if so, issue warning/error

  read(74,*) ncelx, ncely, ncelz ! should be positive
  write(*,'(a,2i5)') 'Lateral supercell size ', ncelx, ncely
  write(*,'(a,2i5)') 'Thickness in cells     ', ncelz
  ngam = ncelz
  allocate(gam_in(0:ngam))
  lgam = (ngam/2)+1-ngam
  hgam = ngam/2
  do igam=lgam, hgam
    gam_in(igam-lgam+1) = dble(igam)/dble(ngam)*two*pi/cc
  enddo

  read(74,*) nbeams_in
  read(74,*) nfg

  allocate(hklbeams_in(3,nbeams_in),aubeams_in(3,nbeams_in))
  read(74,*) ! empty line
  read(74,*) ! header
  ms_slices = 0
  do ig = 1, nbeams_in
    read(74,*) hklbeams_in(1:3,ig), aubeams_in(1:3,ig)
    if(hklbeams_in(1,ig)==0 .and. hklbeams_in(2,ig)==0 .and. hklbeams_in(3,ig)==0) zerohkl_in=ig
    if(abs(hklbeams_in(3,ig))>ms_slices) ms_slices = abs(hklbeams_in(3,ig))
  enddo
  aubeams_in(:,:) = dble(10)*auinnm*aubeams_in(:,:)

  nc0c = nfg
  allocate(c0c(nc0c),c0ci(2,nc0c))
  read(74,*) ! empty line
  read(74,*) ! header
  ifg = 0
  do ix=1, nfg
    read(74,*) igam, ig, fgre, fgim
    fg = cmplx(fgre,fgim)
    ! compensate for the translated potential
    transf = exp( -two*pi*imag*( dble(hklbeams_in(1,ig)*shiftx)/dble(ncelx*meshx) + dble(hklbeams_in(2,ig)*shifty)/dble(ncely*meshy) ) )
    fg = fg * transf
!    write(91,'(2e18.10,3x,e12.4)') fg, abs(fg)
    if(ix==1) then ! the list is sorted by amplitude
      maxfg=abs(fg)
      write(*,'(a,e14.6)') 'Maximum |Fg|:  ',maxfg
    endif
    if(abs(fg)>(c0_min/(dble(10)*maxd0d*maxd0d*maxfg))) then
      ifg = ifg + 1
      c0c(ifg) = fg
      c0ci(2,ifg) = igam-lgam+1
      c0ci(1,ifg) = ig
    endif
  enddo
  nc0c = ifg

end
