SUBROUTINE rel(i,j)

  USE lattice

  IMPLICIT NONE

  INTEGER,INTENT(IN)::i,j
  INTEGER :: itmp, jtmp

  itmp=GridShift(i,nx)
  jtmp=GridShift(j,ny)

  xij=DBLE(itmp)*dx
  yij=DBLE(jtmp)*dy
  
END SUBROUTINE rel

SUBROUTINE reci(i,j)
  
  USE lattice
  USE param
  
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::i,j
  INTEGER :: itmp, jtmp
  
  itmp=GridShift(i,nx)
  jtmp=GridShift(j,ny)
  
  kxij=DBLE(itmp)*dkx
  kyij=DBLE(jtmp)*dky
  
END SUBROUTINE reci

SUBROUTINE prop(i,j)

  USE param
  USE beam
  USE lattice
  USE propagator

  IMPLICIT NONE

  INTEGER,INTENT(IN)::i,j
  DOUBLE PRECISION :: ksquare

  CALL reci(i,j)

  ksquare=kxij**2.D0+kyij**2.D0
  reprop=COS(PI*ksquare*lambda*dz*kbykz)
  improp=-SIN(PI*ksquare*lambda*dz*kbykz)

END SUBROUTINE prop
