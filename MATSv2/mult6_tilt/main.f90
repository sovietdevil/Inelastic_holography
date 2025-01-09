PROGRAM main
  
  USE beam
  USE lattice
  USE wf
  USE preprocessing
  use potential
  
  IMPLICIT NONE

  LOGICAL:: output3d
  
  NAMELIST/ inputdata / strf,nlx,nly,nlz,nslice,beam_m,qmax,qmin,ene,mu,output3d,output2d,vortout,shiftx,shifty,tiltx,tilty,dphi,&
                        df,c12a,c12b,c21a,c21b,c23a,c23b,c3,c32a,c32b,c34a,c34b,c41a,c41b,c43a,c43b,c45a,c45b,&
                        c5,c52a,c52b,c54a,c54b,c56a,c56b,preprocess,tmin,tmax,tstp,coff,maxn
                        ! parameters according eq.(1) in chap.3 of Pennycook & Nellist, Springer (2011)
  nlx=1
  nly=1
  nlz=1
  nslice=6
  beam_m=0
  qmax=0.4d0
  qmin=0.0d0
  ene=200.d0
  mu=1.0
  shiftx=0
  shifty=0
  tiltx=0.d0
  tilty=0.d0
  output3d=.false.
  output2d=.false.
  vortout=0
  dphi=0.0
  df=0.0
  c12a=0.0
  c12b=0.0
  c21a=0.0
  c21b=0.0
  c23a=0.0
  c23b=0.0
  c3=0.0
  c32a=0.0
  c32b=0.0
  c34a=0.0
  c34b=0.0
  c41a=0.0
  c41b=0.0
  c43a=0.0
  c43b=0.0
  c45a=0.0
  c45b=0.0
  c5=0.0
  c52a=0.0
  c52b=0.0
  c54a=0.0
  c54b=0.0
  c56a=0.0
  c56b=0.0
  preprocess=.false.
  tmin=210
  tmax=420
  tstp=210
  maxn=1000000
  READ(*,inputdata)
  
  WRITE(*,*)"*******************Code for Multislice simulation*****************************"

  CALL SetupEnergyScale     !module beam
  CALL SetupLattice         !module structure
  CALL InitialWF            !subroutine from wf.f90 file
  call printwfc(0.d0)       !subroutine from output.f90 file - prints initial wavefunction
  CALL MultiSliceProp       !subroutine from multislice.f90 file
  IF(output3d) CALL Get3dWF  !subroutine from 3dwf.f90 file
  IF(preprocess) CALL preprocess_datacube !subroutine from preprocess.f90 file

END PROGRAM main
