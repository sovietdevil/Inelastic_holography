!program to generate tabulated potential
  PROGRAM GenPot

    USE struct
    USE inputs
    USE beams
    USE asff

    IMPLICIT NONE

    LOGICAL :: dw
    INTEGER :: xc,yc     !shifted origin (has to be a number between 0-hmax & 0-kmax
    CHARACTER (30) :: structfile    !name of the file to read the structure

    NAMELIST/inputdata/vhkltype,hmax,kmax,lmax,dw,xc,yc,structfile
    vhkltype='WEICK'
    hmax=42 
    kmax=42
    lmax=42
    xc=0
    yc=0
    dw=.false.
    structfile='Fe.struct'
    READ(*,inputdata)

    !initializing structure and generate lattice
    CALL init_struct(structfile)
    CALL latgen_struct

    IF (vhkltype=='DOYLE') CALL init_asff

    CALL GetPot(xc,yc,dw)
    
  END PROGRAM GenPot


  !write the potential in the designated file***********************************************
  SUBROUTINE GetPot(xc,yc,dw)

    USE fftw3
    USE defs
    USE beams

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: xc,yc
    LOGICAL,INTENT(IN) :: dw

    INTEGER :: h,k,l,ht,kt,lt,ctr,i,j
    DOUBLE PRECISION :: vre,vim
    double complex :: vtmp
    CHARACTER(20) :: fname,fname1
    DOUBLE COMPLEX :: vhkl         !function to get potential for hkl plane
    DOUBLE PRECISION :: hkl_length !function to get hkl_length
    DOUBLE PRECISION, ALLOCATABLE :: repot(:,:,:),impot(:,:,:),rpot(:,:,:),ipot(:,:,:)
    DOUBLE PRECISION :: dwexp, dwfac
    !fftw variables
    TYPE(C_PTR) :: plan
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: in(:,:,:),out(:,:,:)

    IF(dw)THEN
!       dwfac = 0.0039d0*(nminau/(two*pi))**two  ! hcp cobalt, Aust. J. Phys., 1988,41,461-8
       dwfac = 0.0035d0*(nminau/(two*pi))**two ! bccFe
    ELSE
       dwfac = 0.D0
    END IF

    ALLOCATE(in(0:hmax-1,0:kmax-1,0:lmax-1))
    ALLOCATE(out(0:hmax-1,0:kmax-1,0:lmax-1))
    ALLOCATE(repot(0:hmax-1,0:kmax-1,0:lmax-1))
    ALLOCATE(impot(0:hmax-1,0:kmax-1,0:lmax-1))
    ALLOCATE(rpot(0:hmax-1,0:kmax-1,0:lmax-1))
    ALLOCATE(ipot(0:hmax-1,0:kmax-1,0:lmax-1))

    in=CMPLX(0.D0,0.D0)
    out=CMPLX(0.D0,0.D0)
    rpot=0.D0
    ipot=0.D0
    repot=0.D0
    impot=0.D0

    !structure factor and atomic form factor calculation
    DO h=0,hmax-1
print *, h
       IF(h.LE.hmax/2)THEN
          ht=h
       ELSE
          ht=h-hmax
       END IF
       DO k=0,kmax-1
       IF(k.LE.kmax/2)THEN
          kt=k
       ELSE
          kt=k-kmax
       END IF
          DO l=0,lmax-1
             IF(l.LE.lmax/2)THEN
                lt=l
             ELSE
                lt=l-lmax
             END IF
             dwexp=EXP(-dwfac*hkl_length(ht,kt,lt)**2.D0)
             vtmp = dwexp*vhkl(ht,kt,lt)
             !vre=REAL(vhkl(ht,kt,lt))*dwexp
             !vim=IMAG(vhkl(ht,kt,lt))*dwexp
             !out(h,k,l)=CMPLX(vre,vim)
             out(h,k,l)=vtmp
          END DO
       END DO
    END DO

    !backward transform to get the potential in real space
!    plan = FFTW_PLAN_DFT_3D(lmax, kmax, hmax, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
!    CALL FFTW_EXECUTE_DFT(plan, out, in)
    call dfftw_plan_dft_3d(plan, hmax, kmax, lmax, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, out, in)

!    in=in/SQRT(DBLE(hmax*lmax*kmax))

    DO h=0,hmax-1
       DO k=0,kmax-1
          DO l=0,lmax-1
             repot(h,k,l)=REAL(in(h,k,l))
             impot(h,k,l)=IMAG(in(h,k,l))
          END DO
       END DO
    END DO

    !shift the origin (unit of energy is V or J/C)
    DO l=0,lmax-1
       ht=xc
       DO h=0,hmax-1
          kt=yc
          DO k=0,kmax-1
             lt=l
             rpot(h,k,l)=repot(ht,kt,lt)
             ipot(h,k,l)=impot(ht,kt,lt)
             kt=kt+1
             IF(kt.GT.kmax-1)THEN
                kt=0
             END IF
          END DO
          ht=ht+1
          IF(ht.GT.hmax-1)THEN
             ht=0
          END IF
       END DO
    END DO
    
    !write in the coulomb potential file (unit of energy is V or J/C)
    DO l=0,lmax-1
       IF(l.LT.10)THEN
          WRITE(fname,'(I1,a)')l,'.coul'
          WRITE(fname1,'(I1,a)')l,'.plt'
       ELSE if(l.lt.100) then
          WRITE(fname,'(I2,a)')l,'.coul'
          WRITE(fname1,'(I2,a)')l,'.plt'
       else if(l.lt.1000) then
          WRITE(fname,'(I3,a)')l,'.coul'
          WRITE(fname1,'(I3,a)')l,'.plt'
       else
          WRITE(fname,'(I4,a)')l,'.coul'
          WRITE(fname1,'(I4,a)')l,'.plt'
       END IF
       OPEN(UNIT=2,FILE=fname)
!       OPEN(UNIT=3,FILE=fname1)
       DO h=0,hmax-1
          DO k=0,kmax-1
             WRITE(2,'(2I6,2E16.6)')h,k,rpot(h,k,l),ipot(h,k,l)
!             WRITE(3,'(2I6,2E16.6)')h,k,rpot(h,k,l),ipot(h,k,l)
!             IF(k.EQ.kmax-1)THEN
!                WRITE(3,*)
!             END IF
          END DO
       END DO
       CLOSE(2)
!       CLOSE(3)
    END DO

!    CALL FFTW_DESTROY_PLAN(plan)
    call dfftw_destroy_plan(plan)

  END SUBROUTINE GetPot
