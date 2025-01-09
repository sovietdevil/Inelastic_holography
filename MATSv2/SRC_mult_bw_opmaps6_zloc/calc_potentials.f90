subroutine calc_potentials

  use defs
  use struct
  use beams
  use asff
  use inputs
  
  implicit none
  
  integer i
  double precision coef, dwfac, dwexp
  double complex vhkl
  double precision hkl_length

!  dwfac = 0.0034d0*(nminau/(two*pi))**two  ! fcc nickel around RT
!  dwfac = 0.0039d0*(nminau/(two*pi))**two  ! hcp cobalt, Aust. J. Phys., 1988,41,461-8
  dwfac = 0.0035d0*(nminau/(two*pi))**two ! bccFe
!  dwfac = 0.0090d0*(nminau/(two*pi))**two ! fccAl
!  dwfac = 0.0020d0*(nminau/(two*pi))**two ! Mn3O4, JPC 7 (1974), 409
!  dwfac = 0.0050d0*(nminau/(two*pi))**two ! SrTiO3 according to http://www-hrem.msm.cam.ac.uk/hrem/local/ems/examples/bu1.html
!  dwfac = 0.00666d0*(nminau/(two*pi))**two ! Ga in GaAs
!  dwfac = 0.0027930d0*(nminau/(two*pi))**two ! SiC, Si
!  dwfac = 0.000027930d0*two * nminau**two ! SiC, Si - corrected from <u^2>
!  dwfac=0.0d0
  
!  allocate( ug_in(nbeams_in) )
  allocate( ug_out(nbeams_out) )

! prepare coefficient, so that from vhkl in Volts we get uhkl is in 1/a.u.^2
  coef= four*pi*auinnm**two/fuccon
  
!  write(*,*)
!  write(*,'(a)') "           h  k  l     Re[Ug]     Im[Ug]  "
!  write(*,'(a)') "------------------------------------------"

  do i = 1, nbeams_out
    dwexp = exp(-dwfac*hkl_length(hklbeams_out(1,i),hklbeams_out(2,i),hklbeams_out(3,i))**two)
    ug_out(i) = dwexp*coef*gamma_cr_out*vhkl(hklbeams_out(1,i),hklbeams_out(2,i),hklbeams_out(3,i))
!    write(*,'("Outgoing: ",3i3,2g14.6)') hklbeams_out(1,i),hklbeams_out(2,i),hklbeams_out(3,i), ug_out(i)
!    write(*,'("Outgoing: ",3f14.6,2g14.6)') aubeams_out(1,i),aubeams_out(2,i),aubeams_out(3,i), ug_out(i)
  enddo
 
end
