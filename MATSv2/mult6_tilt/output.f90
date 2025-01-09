!*******************************************************************************
SUBROUTINE OpenFile

  IMPLICIT NONE

  INTEGER:: i,j	

  !File to write the angular momentum
  OPEN(UNIT=4,FILE='orbl.plt',STATUS='replace')
  CLOSE(4)

  !File to write the intensity in hkl space
  OPEN(UNIT=5,FILE='hklint.plt',STATUS='replace')
  WRITE(5,'(a,$)')'#d'
  DO i=0,10
	DO j=0,10
	WRITE(5,'(I3,$)')i+j
	END DO
  END DO
  WRITE(5,*)
  CLOSE(5)

END SUBROUTINE OpenFile

!******************************************************************************
SUBROUTINE WriteIntensity(k,thkn)

  USE lattice
  USE wf

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(IN)::thkn
  INTEGER, INTENT(IN)::k

  INTEGER::ktmp

  ktmp=1*unz                    !output the wavefunction every (orig. 4th) unit cell

  CALL print_intensity(thkn)
  
  CALL print_hkl_intensity(thkn)

  IF(output2d)THEN
     IF(MOD(k,ktmp).EQ.0)CALL printwfc(thkn)
  END IF

END SUBROUTINE WriteIntensity

!***************************************************************************
SUBROUTINE print_intensity(thkn)

  USE lattice
  USE wf

  IMPLICIT NONE
  
  DOUBLE PRECISION,INTENT(IN)::thkn

  INTEGER :: i,j
  DOUBLE PRECISION :: totint
  
  totint=0.D0
  DO i=0,nx-1
     DO j=0,ny-1
        totint=totint + ABS(rwf(i,j))**2.D0
     END DO
  END DO
  
  WRITE(*,*)"=========================================================================="
  WRITE(*,'(a,f16.6,a)')"Propagated by                      ", thkn, " Bohr"
  WRITE(*,'(a,f16.6)')"Total intensity in real space      ", totint
  
  totint=0.D0
  DO i=0,nx-1
     DO j=0,ny-1
        totint=totint + ABS(kwf(i,j))**2.D0
     END DO
  END DO
  
  WRITE(*,'(a,f16.6)')"Total intensity in reciprocal space", totint

END SUBROUTINE print_intensity

!***************************************************************************
SUBROUTINE print_hkl_intensity(thkn)
  
  USE lattice
  USE wf
  
  IMPLICIT NONE
  
  DOUBLE PRECISION,INTENT(IN)::thkn

  INTEGER :: open_stat
  
  INTEGER :: i,j
  DOUBLE PRECISION :: tmpi, tmpr, tmp, hklsquare

  OPEN(UNIT=5,FILE='hklint.plt',STATUS='old', &
       IOSTAT=open_stat,POSITION='append')
  IF(open_stat.NE.0)THEN
     WRITE(*,*)"hklint.plt file not present !!"
     STOP
  ELSE
     WRITE(5,'(e14.6,$)')thkn
     DO i=0,10
        DO j=0,10
           WRITE(5,'(e14.6,$)')ABS(kwf(i,j))**2.D0
        END DO
     END DO
     WRITE(5,*)
     CLOSE(5)
  END IF
  
END SUBROUTINE print_hkl_intensity

!******************************************************************************
SUBROUTINE printwf(thkn)
  
  USE lattice
  USE wf
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: thkn

  INTEGER :: i,j,itmp,jtmp,ictr,jctr
  real :: thkn1
  CHARACTER(40)::fname

  thkn1=thkn*0.5292/10

  WRITE(fname,'(f6.1,a)')thkn1,'nm.plt'
  
  ictr=0
  OPEN(UNIT=2,FILE=fname)
  DO i=nx/2+1,3*nx/2
     itmp=GridShift2(i,nx)
     jctr=0
     ictr=ictr+1
     DO j=ny/2+1,3*ny/2
        jctr=jctr+1
        jtmp=GridShift2(j,ny)
        WRITE(2,'(2I6,3e14.5)')ictr,jctr,ABS(rwf(itmp,jtmp))**2.D0,ABS(kwf(itmp,jtmp))**2.D0
     END DO
     WRITE(2,*)
  END DO
  CLOSE(2)

END SUBROUTINE printwf

!******************************************************************************
SUBROUTINE printwfc(thkn)
  
  USE lattice
  USE wf
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(IN) :: thkn

  INTEGER :: i,j,itmp,jtmp,ictr,jctr
  real :: thkn1
  CHARACTER(40)::fname

  thkn1=thkn*0.5292/10

  WRITE(fname,'(f6.1,a)')thkn1,'nm.plt'
  
  ictr=0
  OPEN(UNIT=2,FILE=fname)
  DO i=nx/2+1,3*nx/2
     itmp=GridShift2(i,nx)
     jctr=0
     ictr=ictr+1
     DO j=ny/2+1,3*ny/2
        jctr=jctr+1
        jtmp=GridShift2(j,ny)
        WRITE(2,'(2I6,4e14.5)')ictr,jctr,rwf(itmp,jtmp),kwf(itmp,jtmp)
     END DO
     WRITE(2,*)
  END DO
  CLOSE(2)

END SUBROUTINE printwfc

!***************************************************************************
SUBROUTINE print_3dwf_real

  USE lattice
  USE wf

  IMPLICIT NONE

  INTEGER :: i,j,k,nzr

  nzr=nslice*nlz !total number of slices in the simulation box to get the 3d wave function

  WRITE(70,'(a)') ' real part of 3D wavefunction in real space'
  WRITE(70,'(a)') ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
  WRITE(70,'(i5,3f12.6)')       1, 0.000000, 0.000000, 0.000000
  WRITE(70,'(i5,3f12.6)')      nx, dx,       0.000000, 0.000000
  WRITE(70,'(i5,3f12.6)')      ny, 0.000000, dy,       0.000000
  WRITE(70,'(i5,3f12.6)')     nzr, 0.000000, 0.000000, az0/dble(nslice)
  WRITE(70,'(a)')'    1   0.000000  0.000000  0.000000  0.000000'

  WRITE(71,'(a)') ' imaginary part of 3D wavefunction in real space'
  WRITE(71,'(a)') ' OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
  WRITE(71,'(i5,3f12.6)')       1, 0.000000, 0.000000, 0.000000
  WRITE(71,'(i5,3f12.6)')      nx, dx,       0.000000, 0.000000
  WRITE(71,'(i5,3f12.6)')      ny, 0.000000, dy,       0.000000
  WRITE(71,'(i5,3f12.6)')     nzr, 0.000000, 0.000000, az0/dble(nslice)
  WRITE(71,'(a)')'    1   0.000000  0.000000  0.000000  0.000000'

  DO i= 0, nx-1
     DO j= 0, ny-1
       WRITE(70,'(6e14.6)')(REAL(r3wf(i,j,k)),k=0,nzr-1)
       WRITE(71,'(6e14.6)')(IMAG(r3wf(i,j,k)),k=0,nzr-1)
       write(70,*)
       write(71,*)
     END DO
  END DO

END SUBROUTINE print_3dwf_real

!***************************************************************************
SUBROUTINE print_angm(thkn)
  
  USE angm
  
  IMPLICIT NONE
  
  DOUBLE PRECISION,INTENT(IN)::thkn
  
  INTEGER :: open_stat
  
  OPEN(UNIT=4,FILE='orbl.plt',STATUS='old', &
       IOSTAT=open_stat,POSITION='append')
  IF(open_stat.NE.0)THEN
     WRITE(*,*)"orbl.plt file not present !!"
     STOP
  ELSE	
     WRITE(4,'(3E16.6)')thkn,reorb,imorb
     CLOSE(4)
  END IF
END SUBROUTINE print_angm


!  DO i=nx/2+1,nx+nx/2
!     DO j=ny/2+1,ny+ny/2
!        IF(i.LT.nx)THEN
!           itmp=i
!        ELSE
!           itmp=i-nx
!        END IF
!        IF(j.LT.ny)THEN
!           jtmp=j
!        ELSE
!           jtmp=j-ny
!        END IF
!        CALL rel(itmp,jtmp)
