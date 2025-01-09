program dyndif

  use struct
  use defs
  use inputs
  use beams
  use asff

  implicit none
  
  double precision uzero
  double complex vhkl

  call deffile
  call init_struct
  call latgen_struct
  call read_inputs

  if (vhkltype=='WIEN ') call read_vhklwien
  if (vhkltype=='DOYLE') call init_asff
  uzero = dble(vhkl(0,0,0))
  gamma_cr_in  = one + (ebeam+uzero/dble(1000))/emass
  gamma_cr_out = one + (ebeam-(eloss-uzero)/dble(1000))/emass
  k_in  = emass/(hbarc*nminau)*sqrt(gamma_cr_in**two-one)
  k_out = emass/(hbarc*nminau)*sqrt(gamma_cr_out**two-one)
  write(*,'("Gamma in crystal for in/out: ",f12.9," /",f12.9)') gamma_cr_in, gamma_cr_out
  write(*,'("K(au) in crystal for in/out: ",f12.6," /",f12.6)') k_in, k_out
  write(*,'("l(au) in crystal for in/out: ",f12.6," /",f12.6)') two*pi/k_in, two*pi/k_out

  call generate_beams
  call calc_potentials
  call calc_ws
  call select_beams
  call sort_beams
  call calc_matrices
  call calc_blochs

  call load_multislice

  call sort_blochs

end
