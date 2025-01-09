subroutine vorticity(slice)

  use lattice
  use wf

  implicit none
  
  integer i,j
  integer slice,mslice
  double complex dxwf, dywf
  double precision vort
  character*11 fname
  
  mslice=1000000+slice
  write(fname,'(i7,a4)') mslice,'.dat'
  open(77,name=fname(2:11))
  do i=0,nx-1
    do j=0,ny-1
      call derx(i,j,dxwf)
      call dery(i,j,dywf)
      vort = imag(conjg(dxwf)*dywf)
      write(77,'(i5,i5,f20.12)') i, j, vort
    enddo
  enddo
  close(77)

end subroutine vorticity


SUBROUTINE FindOrbAngMomentum
  
  USE lattice
  USE wf
  USE angm
  
  IMPLICIT NONE
  
  integer i,j
  double complex dxwf, dywf, orb, psi

  orb = dble(0)
  DO i=0,nx-1
     DO j=0,ny-1
        call derx(i,j,dxwf)
        call dery(i,j,dywf)
        call rel(i,j)
        psi = rwf(i,j)
        orb = orb + cmplx(0,-1)*conjg(psi)*(xij*dywf-yij*dxwf)
     END DO
  END DO
  reorb = real(orb)
  imorb = imag(orb)

END SUBROUTINE FindOrbAngMomentum

subroutine derx(i,j,dxwf)

  use wf
  use lattice

  implicit none

  integer, intent(in) :: i, j
  double complex, intent(out) :: dxwf

  integer ip1, im1

  ip1 = i + 1
  if(ip1==nx) ip1 = 0
  im1 = i - 1
  if(im1==-1) im1 = nx - 1

  dxwf = (rwf(ip1,j)-rwf(im1,j))/(dble(2)*dx)

end subroutine derx

subroutine dery(i,j,dywf)

  use wf
  use lattice

  implicit none

  integer, intent(in) :: i, j
  double complex, intent(out) :: dywf

  integer jp1, jm1

  jp1 = j + 1
  if(jp1==ny) jp1 = 0
  jm1 = j - 1
  if(jm1==-1) jm1 = ny - 1

  dywf = (rwf(i,jp1)-rwf(i,jm1))/(dble(2)*dy)

end subroutine dery
