SUBROUTINE Get3dWF

  USE lattice
  USE wf
  USE ioformat
  USE dft

  IMPLICIT NONE

  INTEGER :: nzr

  nzr=nslice*nlz !total number of slices in the simulation box to get the 3d wave function

!  ALLOCATE(k3wf(0:nx-1,0:ny-1,0:nzr-1))

  CALL print_3dwf_real

  !forward transform to get the 3d wf in reciprocal space
!  CALL FFT3d

  !normalize
!  k3wf=k3wf/SQRT(DBLE(nx*ny*nz))

!  CALL print_3dwf_rec

!  DEALLOCATE(r3wf)
!  DEALLOCATE(k3wf)
  
END SUBROUTINE Get3dWF
