SUBROUTINE MultiSliceProp
  
  USE dft
  USE beam
  USE lattice
  USE propagator
  USE wf
  USE angm
  USE potential
  USE param
  
  IMPLICIT NONE
  
  INTEGER :: i,j,k,ktmp,cnt,ui,uj,kslice,kctr,nzr
  DOUBLE PRECISION :: tmpi, tmpr, tmp, thkn

  double complex, allocatable :: aprop(:,:), vcom(:,:)
  
  ! JR: init propagator array
  allocate (vcom(0:unx-1,0:uny-1))
  allocate (aprop(0:nx-1,0:ny-1))
  do i=0,nx-1
    do j=0,ny-1
      call prop(i,j)
      aprop(i,j) = cmplx(reprop,improp)
    enddo
  enddo
  ! JR: prop. init done

  nzr=nslice*nlz !total number of slices in the simulation box to get the 3d wave function

  ALLOCATE(r3wf(0:nx-1,0:ny-1,0:nzr-1))

  ALLOCATE(rev(0:unx-1,0:uny-1))
  ALLOCATE(imv(0:unx-1,0:uny-1))

  kslice=unz/nslice !wf to be saved every kslice

  IF(MOD(unz/2,kslice).NE.0)THEN !check whether middle of cell is one of the slices (to ensure BCC symmetry)
     WRITE(*,*)"Unacceptable # of slices per unic cell to generate the 3D wave function"
     STOP
  END IF
  
  kctr=0 !counter to store the 2D wf to generate the 3D wf
  ktmp=0 !counter to propagate the wave within one unit cell

  CALL OpenFile !open some output files

  !loop to propagate the wave *************************************************
  DO k=0, nz-1
     
     thkn=k*dz !thickness
     
     !output wavefunctions before passing k th slice************************
     IF(MOD(k,unz).EQ.0)THEN
        CALL WriteIntensity(k,thkn)
        ktmp=-1
     END IF
     
     if(vortout>0 .and. mod(k,vortout).eq.0) then
       call vorticity(k)
     endif
     
     !output <Lz> before passing k th slice ***************************
     CALL FindOrbAngMomentum
     CALL print_angm(thkn)
     
     !store 3D wavefunction before passing k th slice *********************
     IF(MOD(k,kslice).EQ.0)THEN
       r3wf(:,:,kctr)=rwf(:,:)
       kctr=kctr+1
     END IF
     
     !All the outputs dumped before the wave passeing k th slice****
     !Now start to propagate the wave through the k th slice*******************
     
     ktmp=ktmp+1 !advance the counter for propagation
     
     !call the potential subroutine
     CALL Pot(ktmp)
     vcom = cmplx(rev,imv)
     
     !multiply WF and potential
     ui=0
     uj=0
     DO i=0,nx-1
        IF(ui.EQ.unx)THEN
           ui=0
        ENDIF
        DO j=0,ny-1
           IF(uj.EQ.uny)THEN
              uj=0
           ENDIF
           rwf(i,j) = rwf(i,j)*vcom(ui,uj)
           uj=uj+1
        ENDDO
        ui=ui+1
     ENDDO
     
     !forward transform
     CALL FFT2d
     
     !multiply the propagator and wf in reciprocal space
     kwf = kwf*aprop
     kwf=kwf/SQRT(DBLE(nx*ny))
     
     !Backward transform to get the propagated wf in real space
     CALL BFT2d
     
     !normalize
     rwf=rwf/SQRT(DBLE(nx*ny))
     
     !forward transform to get the wf in reciprocal space
     !CALL FFT2d
  
     
  END DO
  
  IF(kctr.NE.nzr)THEN
     WRITE(*,*)"WARNING !! 3D wavefunction grid not matching"
     STOP
  ELSE
     WRITE(*,*)"************* (-:Propagation completed :-) *************"
  END IF
  
  DEALLOCATE(rwf)
  DEALLOCATE(kwf)
  DEALLOCATE(rev)
  DEALLOCATE(imv)
  
END SUBROUTINE MultiSliceProp

!******************************************************************************************************
SUBROUTINE Pot(ktmp)
  
  USE potential
  USE lattice
  USE beam
  
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::ktmp

  INTEGER :: i,j,itmp,jtmp
  DOUBLE PRECISION :: tmp
  CHARACTER(40):: vfname

  IF(ktmp.LT.10)THEN
     WRITE(vfname,'(I1,a)')ktmp,'.coul'
  ELSE IF(ktmp.GE.10.AND.ktmp.LT.100)THEN
     WRITE(vfname,'(I2,a)')ktmp,'.coul'
  ELSE IF(ktmp.GE.100.AND.ktmp.LT.1000)THEN
     WRITE(vfname,'(I3,a)')ktmp,'.coul'
  ELSE 
     WRITE(vfname,'(I4,a)')ktmp,'.coul'
  END IF
  
  rev=0.D0
  imv=0.D0
  
!  WRITE(*,*)'Reading potential from ', TRIM(vfname)
  OPEN(UNIT=2,FILE=vfname)
  DO i=0,unx-1
     DO j=0,uny-1
        READ(2,*)itmp,jtmp,tmp    !make sure that the potential has the unit of V or J/C
!        tmp=0.D0             !no potential*********************
        itmp = modulo(i-shiftx,unx)
        jtmp = modulo(j-shifty,uny)
        rev(itmp,jtmp)=COS(tmp*sigma*dz*kbykz)
        imv(itmp,jtmp)=SIN(tmp*sigma*dz*kbykz)
     END DO
  END DO
  CLOSE(2)
  
END SUBROUTINE Pot
