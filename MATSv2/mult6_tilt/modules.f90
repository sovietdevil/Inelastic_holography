MODULE ioformat

  IMPLICIT NONE
  
  CHARACTER(LEN=*), PARAMETER  :: fmt1 = "(a11,f10.4,a5)"
  CHARACTER(LEN=*), PARAMETER  :: fmt2 = "(2e24.16)"
  
END MODULE ioformat
!********************************************************************************************
MODULE param

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: PI=3.1415926535D0
  DOUBLE PRECISION, PARAMETER :: EMASS=9.10938188D-31               !in the unit of kg
  DOUBLE PRECISION, PARAMETER :: RPLANCK=1.054571726D-34            !in the unit of J.s
  DOUBLE PRECISION, PARAMETER :: PLANCK=6.626068D-34                !in the unit of J.s
  DOUBLE PRECISION, PARAMETER :: LSPEED=2.99792458D8                !in the unit of m/s
  DOUBLE PRECISION, PARAMETER :: ECHARGE=1.602176565D-19            !in the unit of coulomb

  DOUBLE PRECISION, PARAMETER :: eVtoJ=1.602176565D-19              !electronvolt to joule
  DOUBLE PRECISION, PARAMETER :: mtoA=1.0D10                        !meter to angstrom
  DOUBLE PRECISION, PARAMETER :: Atom=1.0D-10                       !angstrom to meter
  DOUBLE PRECISION, PARAMETER :: JtoRy=4.587420897D17               !joul to rydberg
  DOUBLE PRECISION, PARAMETER :: mtoBohr=1.889725989D10             !meter to Bohr
  DOUBLE PRECISION, PARAMETER :: Bohrtonm=0.05291772108d0           !bohr to nm

END MODULE param

!********************************************************************************************
MODULE fftw3

  USE, INTRINSIC :: iso_c_binding

  !modern fortran specific declearations
!  INCLUDE 'fftw3.f03'
  INCLUDE 'fftw3.f'

END MODULE fftw3

!legacy fortran specific declearations
!  INCLUDE 'fftw3.f'

!********************************************************************************************
MODULE beam
  
  USE param
  
  IMPLICIT NONE
  
  INTEGER:: beam_m                         !incoming beam Lz
  DOUBLE PRECISION :: ene                  !incoming beam energy
  DOUBLE PRECISION :: qmax, qmin           !radius of aperture
  DOUBLE PRECISION :: mu		   !std deviation for gaussian
  DOUBLE PRECISION :: lambda, sigma, kbykz !wavelength, int.par., k/k_z for tilts
  
CONTAINS

  SUBROUTINE SetupEnergyScale
    
    IMPLICIT NONE
    
    DOUBLE PRECISION :: gam,e1
    
    !input energy in keV
    !ene=200.D0 ! it is an input parameter now
    e1=ene*1000.D0*eVtoJ
    
    gam=1.D0+e1/(EMASS*LSPEED**2.D0)
    lambda=(PLANCK*LSPEED)/SQRT(2.D0*e1*EMASS*LSPEED**2.D0+e1**2.D0)

    sigma=2.D0*PI*EMASS*gam*ECHARGE*lambda/PLANCK**2.D0
    sigma=sigma/mtoBohr  !unit conversion m->Bohr [sigma in the unit of 1/(V*Bohr)]
    
    lambda=lambda*mtoBohr      !unit conversion m->Bohr

    WRITE(*,*)"==========================================================================="
    WRITE(*,'(a,f12.5,a)')"      Electron beam energy ", ene, " keV"
    WRITE(*,'(a,f12.5,a)')"      Electron beam wavelength ", lambda*0.529*100, " pm"
    WRITE(*,'(a,f12.5,a)')"             Value of  m/m0", gam
  END SUBROUTINE SetupEnergyScale
  
END MODULE beam

!********************************************************************************************
MODULE lattice
  
  IMPLICIT NONE

  DOUBLE PRECISION :: ax, ax0, ay, ay0, az, az0 !lattice parameters
  INTEGER :: nlx, nly, nlz                      !number of cells in x, y and z direction
  INTEGER :: nx, ny, nz                         !grid size in total simulation box
  INTEGER :: unx, uny, unz                      !grid size in one unit cell
  INTEGER :: nslice                             !number of slices in one unit cell to generate the 3d wave function
  DOUBLE PRECISION :: xij, yij                  !real space coordinates
  DOUBLE PRECISION :: kxij, kyij                !reciprocal space coordinates
  DOUBLE PRECISION :: dx, dy, dz                !grid spacing in real space
  DOUBLE PRECISION :: dkx, dky, dkz             !grid spacing in reciprocal space
  INTEGER :: dslice
  CHARACTER(80):: strf
  
CONTAINS
  
  SUBROUTINE SetupLattice
    
    IMPLICIT NONE

    CHARACTER(100):: command
    INTEGER :: ctr,ctr1

    !lattice parameter--------------------------------------------
    open(20,file=strf)
    read(20,*)
    read(20,*)
    read(20,*)
    read(20,1020) ax0, ay0, az0
    close(20)
1020 format(6F10.7,10X,F10.7)

    !simulation box length-------------------------------------------
    ax=DBLE(nlx)*ax0
    ay=DBLE(nly)*ay0
    az=DBLE(nlz)*az0
    
    !get the  mesh size in xy plane from potential file--------------------------
    WRITE(command,'(a)')"tail -n 1 0.coul | awk '{print $1}' > junk"
    CALL SYSTEM(command)
    OPEN(UNIT=2,FILE='junk')
    READ(2,*)unx
    CLOSE(2)
    CALL UNLINK('junk')

    unx=unx+1

    WRITE(command,'(a)')"tail -n 1 0.coul | awk '{print $2}' > junk"
    CALL SYSTEM(command)
    OPEN(UNIT=2,FILE='junk')
    READ(2,*)uny
    CLOSE(2)
    CALL UNLINK('junk')

    uny=uny+1

    IF(MOD(unx,2).NE.0.OR.MOD(uny,2).NE.0)THEN
       WRITE(*,*)"Run with even grid size in plane"
       STOP
    END IF

    !Get the number of slices in one unit cell---------------
    WRITE(command,'(a)')'ls -l *.coul | wc -l > junk'
    CALL SYSTEM(command)
    OPEN(UNIT=2,FILE='junk')
    READ(2,*)unz
    CLOSE(2)
    CALL UNLINK('junk')

    dslice=unz/nslice

    IF(MOD(unz,2).NE.0)THEN
       WRITE(*,*)"Run with even grid size perpendicular to plane"
       STOP
    END IF

    IF(MOD(nslice,2).NE.0)THEN
       WRITE(*,*)"Run with even number of slices per unit cell to generate 3D wavefn"
       STOP
    END IF

    nx=nlx*unx                  !grid size in x
    ny=nly*uny                  !grid size in y
    nz=nlz*unz                  !grid size in z

    dx=ax/DBLE(nx)
    dy=ay/DBLE(ny)
    dz=az/DBLE(nz)
    dkx=1.D0/ax
    dky=1.D0/ay
    dkz=1.D0/az

    WRITE(*,*)             "   Simulation box size ", nlx, " x ", nly, " x ", nlz
    WRITE(*,*)             "  Mesh size in one box ", unx, " x ", uny, " x ", unz
    WRITE(*,'(a,3f12.5,a)')"Simulation box length ",  ax, ay, az, "  Bohr"
    WRITE(*,'(a,3f12.5,a)')"            Lat param ",ax0, ay0, az0,"  Bohr"
    WRITE(*,'(a,3f12.5,a)')" Reciprocal lat param ",1.D0/ax0, 1.D0/ay0, 1.D0/az0," 1/Bohr"
    WRITE(*,'(a,3f12.5,a)')"       Real Mesh size ", dx, dy, dz," Bohr"
    WRITE(*,'(a,3f12.5,a)')" Reciprocal Mesh size ", dkx, dky, dkz," 1/Bohr"
    
  END SUBROUTINE SetupLattice

  INTEGER FUNCTION GridShift(i,j) !Function to grid shift by one box length
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i,j     !i is any integer and j is the grid size (nx/ny/nz)
    
    IF(i.LE.j/2)THEN
       GridShift=i
    ELSE
       GridShift=i-j
    END IF
  END FUNCTION GridShift

  INTEGER FUNCTION GridShift2(i,j) !Function to grid shift by one box length
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i,j     !i is any integer and j is the grid size (nx/ny/nz)

    IF(i.LT.j)THEN
       GridShift2=i
    ELSE
       GridShift2=i-j
    END IF
  END FUNCTION GridShift2
  
END MODULE lattice

!********************************************************************************************
MODULE wf

  USE fftw3

  IMPLICIT NONE

  !modern fortran specific declearations 
!  TYPE(C_PTR) :: planfft, planbft, plan
!  logical :: firstfft = .true., firstbft = .true.
!  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: rwf(:,:), kwf(:,:)          !2d wave function in real and reciprocal space
!  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: r3wf(:,:,:), k3wf(:,:,:)    !3d wave function in real and reciprocal space
  LOGICAL :: output2d
  integer*8 planfft, planbft, plan
  logical :: firstfft = .true., firstbft = .true.
  double complex, allocatable :: rwf(:,:), kwf(:,:)          !2d wave function in real and reciprocal space
  double complex, allocatable :: r3wf(:,:,:), k3wf(:,:,:)    !3d wave function in real and reciprocal space

  double precision :: tiltx = 0.d0, tilty = 0.d0, dphi = 0.d0
  double precision df,c12a,c12b,c21a,c21b,c23a,c23b,c3,c32a,c32b,c34a,c34b ! aberration parameters according eq.(1) in chap.3
  double precision c41a,c41b,c43a,c43b,c45a,c45b,c5,c52a,c52b,c54a,c54b,c56a,c56b ! of Pennycook & Nellist, Springer (2011)
  
  integer vortout ! set to nonzero for outputting vorticity each "vortout" slices

END MODULE wf

!legacy fortran specific declearations 
!  integer*8 planfft, planbft, plan
!  logical :: firstfft = .true., firstbft = .true.
!  double complex, allocatable :: rwf(:,:), kwf(:,:)          !2d wave function in real and reciprocal space
!  double complex, allocatable :: r3wf(:,:,:), k3wf(:,:,:)    !3d wave function in real and reciprocal space

!********************************************************************************************
MODULE propagator
  
  IMPLICIT NONE
  
  DOUBLE PRECISION :: reprop,improp
  
END MODULE propagator

!********************************************************************************************
MODULE angm

  IMPLICIT NONE

  DOUBLE PRECISION :: reorb, imorb
  
END MODULE angm

!********************************************************************************************
MODULE potential

  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: rev(:,:),imv(:,:)
  integer shiftx, shifty

END MODULE potential

!********************************************************************************************
MODULE dft
  
  USE wf
  
  IMPLICIT NONE
  
CONTAINS
  
  SUBROUTINE FFT2d
    
    IMPLICIT NONE
    
    INTEGER :: n1,n2,iret

    n1=SIZE(rwf,1)
    n2=SIZE(rwf,2)
    
!    IF(firstfft) THEN
!       planfft = FFTW_PLAN_DFT_2D(n2, n1, rwf, kwf, FFTW_FORWARD, FFTW_ESTIMATE)
!       firstfft = .false.
!    END IF
!    CALL FFTW_EXECUTE_DFT(planfft, rwf, kwf)
    if(firstfft) then
      call dfftw_init_threads(iret)
      print *, iret
      call dfftw_plan_with_nthreads(16)
      call dfftw_plan_dft_2d(planfft, n1, n2, rwf, kwf, FFTW_FORWARD, FFTW_ESTIMATE)
      firstfft = .false.
    endif
    call dfftw_execute_dft(planfft, rwf, kwf)

  END SUBROUTINE FFT2d

!      call dfftw_init_threads(iret)
!      print *, iret
!      call dfftw_plan_with_nthreads(8)

!legacy fortran
!      call dfftw_plan_dft_2d(planfft, n1, n2, rwf, kwf, FFTW_FORWARD, FFTW_ESTIMATE)
!      call dfftw_execute_dft(planfft, rwf, kwf)

!      CALL FFTW_DESTROY_PLAN(planfft)

  SUBROUTINE BFT2d
    
    IMPLICIT NONE
    
    INTEGER :: n1,n2
    
    n1=SIZE(kwf,1)
    n2=SIZE(kwf,2)
    
!    IF(firstbft) THEN
!       planbft = FFTW_PLAN_DFT_2D(n2, n1, kwf, rwf, FFTW_BACKWARD, FFTW_ESTIMATE)
!       firstbft = .false.
!    END IF
!    CALL FFTW_EXECUTE_DFT(planbft, kwf, rwf)
    if(firstbft) then
      call dfftw_plan_dft_2d(planbft, n1, n2, kwf, rwf, FFTW_BACKWARD, FFTW_ESTIMATE)
      firstbft = .false.
    endif
    call dfftw_execute_dft(planbft, kwf, rwf)

  END SUBROUTINE BFT2d

!legacy fortran 
!  call dfftw_plan_dft_2d(planbft, n1, n2, kwf, rwf, FFTW_BACKWARD, FFTW_ESTIMATE)
!  call dfftw_execute_dft(planbft, kwf, rwf)

!  CALL FFTW_DESTROY_PLAN(planbft)

  SUBROUTINE FFT3d
    
    IMPLICIT NONE
    
    INTEGER :: n1,n2,n3

    n1=SIZE(r3wf,1)
    n2=SIZE(r3wf,2)
    n3=SIZE(r3wf,3)

!    plan = FFTW_PLAN_DFT_3D(n3, n2, n1, r3wf, k3wf, FFTW_FORWARD, FFTW_ESTIMATE)
!    CALL FFTW_EXECUTE_DFT(plan, r3wf, k3wf)
!    CALL FFTW_DESTROY_PLAN(plan)
    call dfftw_plan_dft_3d(plan, n1, n2, n3, r3wf, k3wf, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, r3wf, k3wf)
    call dfftw_destroy_plan(plan)

  END SUBROUTINE FFT3d

!legacy fortran
!    call dfftw_plan_dft_3d(plan, n1, n2, n3, r3wf, k3wf, FFTW_FORWARD, FFTW_ESTIMATE)
!    call dfftw_execute_dft(plan, r3wf, k3wf)
!    call dfftw_destroy_plan(plan)

END MODULE dft
