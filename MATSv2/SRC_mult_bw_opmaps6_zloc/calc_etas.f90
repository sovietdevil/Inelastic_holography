subroutine calc_etas

  use defs
  use inputs
  use beams
  use blochs

  implicit none

  integer i, j, k
  double complex vabsorp, vabs
  double precision coef

  coef= four*pi*auinnm**two/fuccon

  allocate( eta_out(nbeams_out) )
  eta_out(:) = zero

  do i=1, nbeams_out
    do j=1, nbeams_out
      vabsorp = coef * gamma_cr_out * &
                vabs(hklbeams_out(1,i)-hklbeams_out(1,j), &
                     hklbeams_out(2,i)-hklbeams_out(2,j), &
                     hklbeams_out(3,i)-hklbeams_out(3,j))
      do k=1, nbeams_out
        eta_out(k) = eta_out(k) + dble(conjg(mat_out(i,k))*mat_out(j,k)*vabsorp)
      enddo
    enddo
  enddo
  eta_out(:) = eta_out(:) / (two*chi_out)

  write(*,'(a)') 'Etas for outgoing beam:'
  write(*,'(f10.6,x,$)') eta_out(:)
  write(*,*)

end
