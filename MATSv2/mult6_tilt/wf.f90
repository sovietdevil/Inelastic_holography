!Set up the initial wave function*************************************************
!Output of the subroutine is initial wavefunction in real space and reciprocal space***************************
SUBROUTINE InitialWF
  
  USE beam
  USE lattice  
  USE wf  
  
  IMPLICIT NONE
  
  INTEGER :: i,j
  
  ALLOCATE(rwf(0:nx-1,0:ny-1))
  ALLOCATE(kwf(0:nx-1,0:ny-1))
  
  rwf=CMPLX(0.D0,0.D0)
  kwf=CMPLX(0.D0,0.D0)
  
  ! calculate tilt factor k/k_z that scales dz in propagator and potential term
  kbykz = 1.d0/sqrt(1.d0-(tiltx*tiltx+tilty*tilty)*lambda*lambda)
print *, kbykz
  ! TODO for large tilts (several degrees): ----------------------------
  ! 1) the k-space WF should become elliptic?
  ! 2) there will be a phase ramp at sample entrance plane in real space
  ! --------------------------------------------------------------------

  IF(ABS(beam_m).GT.99)THEN
     CALL PlaneWave
!     CALL Gaussian
  ELSE
     CALL VBeam
  END IF

END SUBROUTINE InitialWF

!*********************************************************************************
SUBROUTINE PlaneWave
  
  USE beam
  USE lattice
  USE wf
  USE dft
  
  IMPLICIT NONE
  
  INTEGER :: i,j

  !get the wf in reciprocal space
  kwf(0,0)=CMPLX(1.D0,0.D0)

  !backward transform and get the wf in real space
  CALL BFT2d
    
  !normalize
  rwf=rwf/SQRT(DBLE(nx*ny))

  WRITE(*,*)"     Plane wave calculation" 
  WRITE(*,*)"===========================================================================" 
  
END SUBROUTINE PlaneWave

!*********************************************************************************
SUBROUTINE VBeam
  
  USE beam
  USE lattice
  USE wf
  USE dft
  use param, only: pi
  
  IMPLICIT NONE
  
  INTEGER :: i,j
  DOUBLE PRECISION :: tmpr,tmpi,totint
  DOUBLE PRECISION ::  kr, phi, th, th2, th3, th4, th5, th6, chi
  
  !get wf in reciprocal space
  DO i=0,nx-1
     DO j=0,ny-1
        
        CALL reci(i,j)

        kr=SQRT((kxij-tiltx)**2.D0 + (kyij-tilty)**2.D0)
        phi=ATAN2(kyij-tilty,kxij-tiltx) + dphi*pi/180.d0
        
        tmpr=COS(DBLE(beam_m)*phi)
        tmpi=SIN(DBLE(beam_m)*phi)
        
        IF(kr.LT.qmax .and. kr.ge.qmin)THEN
           kwf(i,j)=CMPLX(tmpr,tmpi)
        ELSE
           kwf(i,j)=CMPLX(0.D0,0.D0)
        END IF
!print *, i, j, kxij, kyij

        ! aberration function according to Eq.(1) of chap.3 in Pennycook&Nellist, Springer (2011)
        ! note, lambda, kxij, kyij are in (inverse) Bohr radii -> th is in radians (=dimensionless)
        th = kr*lambda
        th2 = th*th
        chi = -th2*df/2.d0 ! df = -c10
        chi = chi + th2*( c12a*cos(2.d0*phi) + c12b*sin(2.d0*phi) )/2.d0
        th3 = th*th2
        chi = chi + th3*( c21a*cos(     phi) + c21b*sin(     phi) )/3.d0
        chi = chi + th3*( c23a*cos(3.d0*phi) + c23b*sin(3.d0*phi) )/3.d0
        th4 = th*th3
        chi = chi + th4*c3/4.d0
        chi = chi + th4*( c32a*cos(2.d0*phi) + c32b*sin(2.d0*phi) )/4.d0
        chi = chi + th4*( c34a*cos(4.d0*phi) + c34b*sin(4.d0*phi) )/4.d0
        th5 = th*th4
        chi = chi + th5*( c41a*cos(     phi) + c41b*sin(     phi) )/5.d0
        chi = chi + th5*( c43a*cos(3.d0*phi) + c43b*sin(3.d0*phi) )/5.d0
        chi = chi + th5*( c45a*cos(5.d0*phi) + c45b*sin(5.d0*phi) )/5.d0
        th6 = th*th5
        chi = chi + th6*c5/6.d0
        chi = chi + th6*( c52a*cos(2.d0*phi) + c52b*sin(2.d0*phi) )/6.d0
        chi = chi + th6*( c54a*cos(4.d0*phi) + c54b*sin(4.d0*phi) )/6.d0
        chi = chi + th6*( c56a*cos(6.d0*phi) + c56b*sin(6.d0*phi) )/6.d0
!!! just for comparison with Own et al., Ultramic., fig.8
!        th5 = th*th4
!        chi = chi + th5*( 300.d0*df*cos(5.d0*phi) )/5.d0
!!! disabled for now ------------------------------------
        ! add the phase factor due to aberrations: exp(2pi*i*chi/lambda)
        ! note: the conversion factor Bohrtonm means that C-coefficients need to be in nm
        !       which seems to be a common convention
        kwf(i,j) = kwf(i,j)*exp(dcmplx(0,2)*pi*chi/(lambda*Bohrtonm))

     END DO
  END DO

  !normalize
  totint=0.D0
  DO i=0,nx-1
     DO j=0,ny-1
        totint=totint+ABS(kwf(i,j))**2.D0
     END DO
  END DO

  kwf=kwf/SQRT(totint)
  
  !backward transform to get wf in real space
  CALL BFT2d

  !normalize
  rwf=rwf/SQRT(DBLE(nx*ny))

  WRITE(*,'(a,I4,a)')"      m = ",beam_m," calculation"
  WRITE(*,'(a,f5.2,a)')"   qmax = ",qmax," calculation"
  WRITE(*,*)"===========================================================================" 
  
END SUBROUTINE VBeam

!*********************************************************************************
SUBROUTINE Gaussian

  USE beam
  USE lattice
  USE wf
  USE param
  USE dft
  
  IMPLICIT NONE

  INTEGER :: i,j
  DOUBLE PRECISION :: tmpr,totint,tmpc

  DO i=0,nx-1
     DO j=0,ny-1
	CALL rel(i,j)
        tmpc=(xij**2.D0+yij**2.D0)/mu**2.D0
        tmpr=EXP(-tmpc)
        rwf(i,j)=CMPLX(tmpr,0.D0)
     END DO
  END DO

  !normalize
  totint=0.D0
  DO i=0,nx-1
     DO j=0,ny-1
        totint=totint+ABS(rwf(i,j))**2.D0
     END DO
  END DO
 
  rwf=rwf/SQRT(totint)

!forward transform to get wf in reciprocal space
  CALL FFT2d

  !normalize
  kwf=kwf/SQRT(DBLE(nx*ny))

  WRITE(*,*)"     Gaussian calculation" 
  WRITE(*,*)"===========================================================================" 

END SUBROUTINE Gaussian


!*********************************************************************************
!  SUBROUTINE M1

!    IMPLICIT NONE

!    INTEGER :: i,j,itmp,jtmp
!    DOUBLE PRECISION :: STRUVE_H0, STRUVE_H1
!    DOUBLE PRECISION, ALLOCATABLE ::  r(:,:)

!    ALLOCATE(r(0:nx-1,0:ny-1))
!    r=0.D0

!    DO i=0,nx-1
!       DO j=0,ny-1
          !to set up real space grid
!          IF(i.LE.nx/2)THEN
!             itmp=i
!          ELSE
!             itmp=i-nx
!          END IF
!          IF(j.LE.ny/2)THEN
!             jtmp=j
!          ELSE
!             jtmp=j-ny
!          END IF
!          x(i,j)=DBLE(itmp)*dx
!          y(i,j)=DBLE(jtmp)*dy

!          r(i,j)=DSQRT(x(i,j)**2.D0+y(i,j)**2.D0)

!          rewf(i,j)=qmax*(BesJ1(qmax*r(i,j))*STRUVE_H0(qmax*r(i,j))-BesJ0(qmax*r(i,j))*STRUVE_H1(qmax*r(i,j))) &
!               /(4.D0*r(i,j))
!          rewf(0,0)=0.D0
!          
!          imwf(i,j)=COS(phi(i,j))*rewf(i,j)          !WARNING: ALWAYS CODE IMAGINARY PART BEFORE REAL PART HERE
!          rewf(i,j)=-SIN(phi(i,j))*rewf(i,j)
          
!       END DO
!    END DO

!    WRITE(*,*)"     m = 1 calculation"
!    WRITE(*,*)"==========================================================================="           
!  END SUBROUTINE M1  
