subroutine calc_matrices

  use inputs
  use beams
  use struct
  use defs
  use blochs
  
  implicit none
  
  double complex vhkl
  double precision hkl_length
  double precision coef, div1, div2, cart_length2, cart_scaprod2, vec(3), dwfac, dwexp
  integer i, j, k
  
  coef= four*pi*auinnm**two/fuccon
!  dwfac = 0.0034d0*(nminau/(two*pi))**two  ! fcc nickel around RT
!  dwfac = 0.0039d0*(nminau/(two*pi))**two  ! hcp cobalt, Aust. J. Phys., 1988,41,461-8
  dwfac = 0.0035d0*(nminau/(two*pi))**two ! bccFe
!  dwfac = 0.0020d0*(nminau/(two*pi))**two ! Mn3O4, JPC 7 (1974), 409
!  dwfac = 0.0050d0*(nminau/(two*pi))**two ! SrTiO3 according to http://www-hrem.msm.cam.ac.uk/hrem/local/ems/examples/bu1.html
!  dwfac = 0.0090d0*(nminau/(two*pi))**two ! fccAl
!  dwfac = 0.00666d0*(nminau/(two*pi))**two ! Ga in GaAs
! dwfac = 0.00002793d0*two * (nminau**two) ! SiC, Si
!  dwfac=0.0d0
  
  allocate( mat_out(nbeams_out,nbeams_out) )

  do i = 1, nbeams_out
    vec(:) = vchi_out(:)+aubeams_out(:,i)
    div1 = vec(3)
    mat_out(i,i) = (k_out**two - cart_length2(vec)**two)/(-two*div1)
  enddo
  do i = 1, nbeams_out-1
    vec(:) = vchi_out(:)+aubeams_out(:,i)
    div1 = vec(3)
    do j = i+1, nbeams_out
      vec(:) = vchi_out(:)+aubeams_out(:,j)
      div2 = vec(3)
      dwexp = exp(-dwfac*hkl_length(hklbeams_out(1,i)-hklbeams_out(1,j),        &
                                    hklbeams_out(2,i)-hklbeams_out(2,j),        &
				    hklbeams_out(3,i)-hklbeams_out(3,j))**two)
! the inversion of crystal would give U(hkl) -> U(-hkl)
! we employed hopefully a more natural approach, see notes in read_inputs.f90
      mat_out(i,j) = gamma_cr_out * coef * dwexp *             &
                     vhkl(hklbeams_out(1,i)-hklbeams_out(1,j), &
                          hklbeams_out(2,i)-hklbeams_out(2,j), &
                          hklbeams_out(3,i)-hklbeams_out(3,j)) &
		     / (two*sqrt(div1*div2))
      mat_out(j,i) = dconjg(mat_out(i,j))
    enddo
  enddo

  print *, 'matrix 2'

end
