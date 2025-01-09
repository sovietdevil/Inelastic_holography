! We are using a convention: K = 2pi/lambda

module defs
  double precision, parameter     :: CLIGHT= 137.03599976d0
  double precision, parameter     :: HBARC = 0.197326968d0 ! hbar*c in keV.nm
  double precision, parameter     :: HC = 1.2398423007d0 ! h*c in keV.nm
  double precision, parameter     :: EMASS = 510.9989d0 ! in keV
  double precision, parameter     :: PI=     3.1415926535897932d0
  double precision, parameter     :: FUCCON= 0.47877645944d0 ! it is inverse of 2me/(4pi hbar^2) [ 1/V.nm^2 ]
  double precision, parameter     :: AUINNM= 0.05291772108d0 ! atomic unit in nm
  double precision, parameter     :: NMINAU= 18.89726125d0 ! nm in atomic units
  DOUBLE PRECISION, PARAMETER     :: JtoRy=4.587420897D17               !joul to rydberg
  double precision, parameter     :: ECHARGE=1.60217646d-19
  double precision, parameter     :: TEST=   1.D-10
  double precision, parameter     :: HALF=   5.D-1
  double precision, parameter     :: ZERO=   0.0d0
  double precision, parameter     :: ONE=    1.0d0
  double precision, parameter     :: TWO=    2.0d0
  double precision, parameter     :: THREE=  3.0d0
  double precision, parameter     :: FOUR =  4.0d0
  double precision, parameter     :: NINETY= 90.0d0
  double complex, parameter       :: ZEROC=  (zero,zero)
!  double complex, parameter       :: IMAG=   (zero,one)
end module defs


module inputs
! atom index, n- and l-quantum numbers of the edge
  integer                          atom, nc, lc
! energy loss (eV)
  double precision                 eloss, esplit
! incoming beam energy (keV)
  double precision                 ebeam
! relativistic gamma and wave-vectors (a.u.^-1) for incoming and outgoing beam outside and inside crystal
  double precision                 gamma_in, gamma_out, gamma_cr_in, gamma_cr_out
! wave-vector sizes
  double precision                 chi_in, chi_out, k_in, k_out
! wave-vectors (chi/k - outside/inside the crystal)
  double precision                 vchi_in(3), vzero_out(3), vchi_out(3), vk_in(3), vk_out(3)
! zone axis in a,b,c
  integer                          za(3)
! zone axis as normalized vector in cartesian system
  double precision                 cart_za(3)
! Laue circle center in a*,b*,c* and as vector in cartesian system
  double precision                 ilcc(3), cart_ilcc(3)
! outgoing beam angles with respect to zone axis (mrad), x given by projection of hkl_det to plane perp. to vchi_in
  double precision                 thetax, thetay
  integer                          hkl_det(3)
! Laue circle center for outgoing beam corresponding to thetax, thetay
  double precision                 olcc(3), cart_olcc(3)
! sample thickness (nm)
  double precision                 tmin, tstep
  integer                          nt
! incoming and outgoing beams and their lengths (outside the crystal) in 1/a.u.
  double precision                 ibeam(3), obeam(3), ib_len, ob_len
end module inputs


module beams
! zero order Laue zone or higher order Laue zones?
  logical                       :: holz = .true.
! modes for beam selection AUTO/FILE
  character*4                      bm_eig, bm_sum
! max |h|,|k|,|l| for beam selection
  integer                          hmax, kmax, lmax
! max |G| in Vrec^(1/3), max w=sq*ksi (dimensionless), min C0^(j)
  double precision                 g_max, w_max, c0_min, w_min
! hkl components of excited beams
  integer, pointer              :: hklbeams_in(:,:), hklbeams_out(:,:)
! beams coordinates in cartesian system (in a.u.)
  double precision, pointer     :: aubeams_in(:,:),  aubeams_out(:,:)
  integer                          nbeams_in, nbeams_out, matsize_in, matsize_out
! subset of beams selected for summation
  integer                          nb_in, nb_out
! index of the hkl = 000 beam
  integer                          zerohkl_in, zerohkl_out
! potential type ('DOYLE', 'WEICK', 'WIEN ') for REAL part
  character*5                      vhkltype
! potential type ('WEICK', 'COEFF') for IMAGINARY part
  character*5                      abstype
  double precision                 abscoeff
  logical                          doabsorp
! Fourier components of potential
  double complex, allocatable   :: ug_in(:), ug_out(:)
! excitation errors (||za), extinction distances and their products
  double precision, allocatable :: sg_in(:), ksi_in(:), sg_out(:), ksi_out(:), wg_in(:), wg_out(:)
! WIEN2k potentials (k; dis, Re, Im) & (k; h, k, l)
  double precision, allocatable :: w2k_vg(:,:)
  integer, allocatable          :: w2k_hkl(:,:)
  integer                       :: w2k_num
end module beams


module blochs
! matrices for eigenvalue problem
  double complex, allocatable   :: mat_in(:,:), mat_out(:,:)
! Bloch wave coefficients for incoming (Cgj) and outgoing (Dhl) beam
  double precision, allocatable :: cgj(:,:), dhl(:,:)
! eigenvalues (gammas)
  double precision, allocatable :: gam_in(:), gam_out(:)
! absorption (etas)
  double precision, allocatable :: eta_in(:), eta_out(:)
! number of Bloch waves (currently unused)
  integer                          nblochs_in, nblochs_out
! auxiliary arrays for sorting Bloch-coefficient products
  integer                          nc0c, nd0d
  double complex, allocatable   :: c0c(:), d0d(:)
  integer, allocatable          :: c0ci(:,:), d0di(:,:) ! 2nd index: 1->G, 2->(j)
! auxiliary arrays for sorting Bloch-coefficient double products
  integer                          ncd
  double complex, pointer       :: cd(:)
  integer, pointer              :: cdi(:,:)  ! 2nd index: 1->G, 2->(j), 3->H, 4->(l)
! arrays of sorted Bloch-coefficients quadruple products and G,j,H,l,G',j',H',l' indices
  integer                       :: nccdd
  double complex, pointer       :: ccdd(:)
  integer, pointer              :: ccddi(:,:)
  double precision              :: gam(4), eta(4) ! for final summation only
! auxiliary array to speed up sorting of cd array
  integer, pointer              :: qsv(:)
end module blochs


module qvectors
! list of q-vectors and corresponding coefficients
  integer                          nqqpr
  double precision, pointer     :: qqpr(:,:)
! multiplicities of q-vectors
  integer, allocatable          :: qvc(:), qvs(:)
  integer                          nqv
end module qvectors


module mdffs
! array of precalculated mdffs (e_index,qqpr_index,edge)
  double complex, allocatable   :: mdff(:,:,:)
! total scattering cross-section (e_index, t_index, edge) & dscs from direct terms only
! we use complex for implementation reasons
  double complex, allocatable   :: dscs(:,:,:)
! coefficients in front of the mdffs (iqqpr, it, conjg/nonconjg)
  double complex, pointer       :: coef(:,:,:)
end module mdffs


module energy
! energy mesh (eV) with respect to Fermi level (for tetra)
  double precision                 emin_mesh, estep
  integer                          ne_mesh
! actual energy mesh for calculation
  double precision                 emin, emax
  integer                          iemin, iemax, ne
! Fermi level (Ry)
  double precision                 efermi
end module energy


module asff 
! atom scattering formfactors

  integer, parameter     :: maxatn = 120
  double precision          fa(4,maxatn), fb(4,maxatn)
  character*2               atnam(maxatn)
  integer                   isour(maxatn)

 contains 

  subroutine init_asff

    implicit none

    double precision fab(8)
    integer i,j,k
    character*2 nam
    
    open(unit=90,err=66,status='old',file='atoms.sff')
 10 read(90,'(A2,2I4,8F9.3)',end=66) nam,i,k,(fab(j),j=1,8)
      isour(k) = i
      atnam(k) = nam
      do i=1,4
        fa(i,k) = fab(2*i-1)/10.0d0
        fb(i,k) = fab(2*i)/100.0d0
      enddo
    goto 10
 66 close(90)

  end subroutine init_asff

end module asff


	module reallocate
	  !     nur 1 (generischer) Name wird von aussen angesprochen
	  interface doreallocate
	    module procedure doreallocate_r8_d1
	    module procedure doreallocate_r8_d2
	    module procedure doreallocate_c16_d1
	    module procedure doreallocate_c16_d2
	    module procedure doreallocate_i4_d1
	    module procedure doreallocate_i4_d2
	    module procedure hugo     !   ;)
	  end interface
	contains

	  !     leider sind mehrere subroutines notwendig fuer verschiedene Typen
	  subroutine doreallocate_r8_d1(tf, newdimension)
	    real*8, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
	    !     nur 1 mal kopieren reicht
	    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
	    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d1(tf, newdimension)
	    complex*16, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d2(tf, newdimension1, newdimension2)
	    complex*16, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d1(tf, newdimension)
	    integer*4, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension))
            min1=min(newdimension,size(tf,1))
	    hilfsfeld(1:min1)=tf(1:min1)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
	    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2))
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
	    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	!     Es gibt auch Methoden, um das Programm unleserlich zu machen :-)
	!     das sollten wir besser vermeiden ;-)
	  subroutine hugo
	    write(*,*) " Hier ist Hugo"
	  end subroutine
	end module 

module struct
! Georg Madsen 2002
  character*80                    :: title
  character*4                     :: lattic,irel
  double precision                :: aa,bb,cc,pia(3),alpha(3)
  double precision                :: vol

  integer                         :: nat,ndif
  integer, allocatable            :: mult(:),jrj(:),iatnr(:),isplit(:)
  double precision, allocatable   :: r0(:),dx(:),rmt(:),zz(:)
  character*10, allocatable       :: aname(:)
  double precision, pointer       :: pos(:,:)

  integer                         :: nsym,iord
  integer, pointer                :: iz(:,:,:),inum(:)
  double precision, pointer       :: tau(:,:)
  logical                         :: ortho

!     ROTLOC(3,3,JATOM) 
!     ROTATE FROM THE GENERAL COORDINATION-SYSTEM INTO 
!     LOCAL SYSTEMS WITH SPECIFIED POINTSYMMETRY.
  double precision, allocatable   :: rotloc(:,:,:)

!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS INTO
!                 CARTESIAN SYSTEM                         
!     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-      
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )       
!                 TO CARTESIAN SYSTEM                                   
  double precision                :: br1_rec(3,3),br2_rec(3,3)
  double precision                :: br1_dir(3,3),br2_dir(3,3)
  double precision                :: g(3,3), grec(3,3)

  double precision                :: vbr(3,4)
  integer                         :: nvbr

CONTAINS
  subroutine init_struct(structfile)

    use reallocate
    use defs

    implicit none

    character(30),intent(in):: structfile

    integer                       :: ios
    double precision              :: test0

    integer                       :: index,i,j,j1,j2,m,jatom

    test0=1.D-5

    open(unit=20,status='old',file=structfile,iostat=ios)
    if(ios.ne.0)then
       write(*,*)'Structure file not found !!!!!!!!!'
       stop
    end if
    read (20,1000) title
    read (20,1010) lattic,nat,irel
    allocate(aname(nat),mult(0:nat),jrj(nat),r0(nat),dx(nat),rmt(nat),zz(nat),rotloc(3,3,nat),iatnr(nat),isplit(nat))
    allocate (pos(3,4096*nat))
    mult(0)=0
    read (20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
! Let's work in nm
!    aa = aa * 0.0529177d0
!    bb = bb * 0.0529177d0
!    cc = cc * 0.0529177d0
! ----------------
    if(ABS(ALPHA(1)).LT.test0) ALPHA(1)=ninety
    if(ABS(ALPHA(2)).LT.test0) ALPHA(2)=ninety
    if(ABS(ALPHA(3)).LT.test0) ALPHA(3)=ninety
    INDEX=0
    do jatom=1,NAT
       INDEX=INDEX+1
print *, jatom, nat
       read(20,1030,iostat=ios) iatnr(jatom),( pos(j,index),j=1,3 ), &
            mult(jatom),isplit(jatom) 
       if(ios /= 0) then
          write(6,*) iatnr(jatom),( pos(j,index),j=1,3 ), &
               mult(jatom),isplit(jatom) 
          write(6,*) 'ERROR IN STRUCT FILE read'
          stop
       endif
       if (mult(jatom) .EQ. 0) then
          write (6,6000) jatom, index, mult(jatom)
          stop
       endif
       do m=1,mult(jatom)-1                                     
          index=index+1                                            
          read(20,1031) iatnr(jatom),( pos(j,index),j=1,3)         
       enddo
       read(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)           
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )           
       read(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
    enddo
    ndif=index
    call doreallocate(pos, 3, ndif)
    read(20,1151) iord
    nsym=iord
    allocate(iz(3,3,nsym),tau(3,nsym),inum(nsym))
    do j=1,iord
       read(20,1101) ( (iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3 ),inum(j)
    enddo
  
1000 format(A80)                                                       
1010 format(A4,23X,I3,/,13X,A4)                                 
1020 format(6F10.7,10X,F10.7)                                          
1030 format(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I4,15X,I2)          
1031 format(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
1050 format(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
1051 format(20X,3F10.8)                                             
1101 format(3(3I2,F10.8/),I8)
1151 format(I4)
6000 format(///,3X,'ERROR IN READING STRUCT : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  end subroutine init_struct

  subroutine latgen_struct
    USE defs
    IMPLICIT NONE
    integer              :: i,j
    double precision     :: sqrt3,rvfac,sinab,sinbc,cosab,cosac,cosbc,wurzel,x1,y1,z1
    double precision     :: det, arec, brec, crec, rcosab, rcosac, rcosbc, sv, vo, vvol

    SQRT3=SQRT(3.D0)
    ALPHA(1)=ALPHA(1)*PI/180.0D0                                             
    ALPHA(2)=ALPHA(2)*PI/180.0D0                                             
    ALPHA(3)=ALPHA(3)*PI/180.0D0                                             
    PIA(1)=2.D0*PI/AA                                                 
    PIA(2)=2.D0*PI/BB                                                 
    PIA(3)=2.D0*PI/CC
    cosab=COS(alpha(3))
    cosac=COS(alpha(2))
    cosbc=COS(alpha(1))
    sinab=SIN(alpha(3))
    sinbc=SIN(alpha(1))

! Calculation of metric tensors
    sv = (alpha(1) + alpha(2) + alpha(3)) / 2.0d0
    vo = sin(sv)*sin(sv-alpha(1))*sin(sv-alpha(2))*sin(sv-alpha(3))
    vol  = 2.0d0 * aa * bb * cc * sqrt(vo)
    arec = bb * cc * sin(alpha(1)) / vol
    brec = cc * aa * sin(alpha(2)) / vol
    crec = aa * bb * sin(alpha(3)) / vol
    rcosbc = (cosac*cosab-cosbc)/(sin(alpha(2))*sin(alpha(3)))
    rcosac = (cosbc*cosab-cosac)/(sin(alpha(3))*sin(alpha(1)))
    rcosab = (cosbc*cosac-cosab)/(sin(alpha(1))*sin(alpha(2)))
    
    g(1,1)=aa*aa
    g(1,2)=aa*bb*cosab
    g(1,3)=aa*cc*cosac
    g(2,1)=g(1,2)
    g(2,2)=bb*bb
    g(2,3)=bb*cc*cosbc
    g(3,1)=g(1,3)
    g(3,2)=g(2,3)
    g(3,3)=cc*cc

    grec(1,1)=arec*arec
    grec(1,2)=arec*brec*rcosab
    grec(1,3)=arec*crec*rcosac
    grec(2,1)=grec(1,2)
    grec(2,2)=brec*brec
    grec(2,3)=brec*crec*rcosbc
    grec(3,1)=grec(1,3)
    grec(3,2)=grec(2,3)
    grec(3,3)=crec*crec
! ------------------------------   

    br1_rec=zero; br2_rec=zero

    select case (LATTIC(1:1))
    case ('H')       !.....HEXAGONAL LATTICE
       br1_rec(1,1)=2.d0/sqrt3*pia(1)                                        
       br1_rec(1,2)=1.d0/sqrt3*pia(1)                                        
       br1_rec(2,2)=pia(2)                                                   
       br1_rec(3,3)=pia(3)                                                   

       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       rvfac=2.d0/SQRT(3.d0)                                             
       ortho=.FALSE.
    case ('S','P' )
       wurzel=SQRT(sinbc**2-cosac**2-cosab**2+2*cosbc*cosac*cosab)
       br1_rec(1,1)= sinbc/wurzel*pia(1)
       br1_rec(1,2)= (-cosab+cosbc*cosac)/(sinbc*wurzel)*pia(2)
       br1_rec(1,3)= (cosbc*cosab-cosac)/(sinbc*wurzel)*pia(3)
       br1_rec(2,2)= pia(2)/sinbc
       br1_rec(2,3)= -pia(3)*cosbc/sinbc
       br1_rec(3,3)= pia(3)
       !
       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       rvfac= 1.d0/wurzel
       ortho=.TRUE.
       if(ABS(alpha(1)-pi/2.d0).GT.test) ortho=.FALSE.
       if(ABS(alpha(2)-pi/2.d0).GT.test) ortho=.FALSE.
       if(ABS(alpha(3)-pi/2.d0).GT.test) ortho=.FALSE.
    case ( 'F')
       BR1_REC(1,1)=PIA(1)                                                   
       BR1_REC(2,2)=PIA(2)                                                   
       BR1_REC(3,3)=PIA(3)                                                   

!     definitions according to column, rows convention for BR2_REC
       BR2_REC(1,1)=-PIA(1)                                                  
       BR2_REC(1,2)= PIA(1)                                                  
       BR2_REC(1,3)= PIA(1)                                                  
       BR2_REC(2,1)= PIA(2)                                                  
       BR2_REC(2,2)=-PIA(2)                                                  
       BR2_REC(2,3)= PIA(2)                                                  
       BR2_REC(3,1)= PIA(3)                                                  
       BR2_REC(3,2)= PIA(3)                                                  
       BR2_REC(3,3)=-PIA(3)                                                  
       !                                                                       
       RVFAC=4.D0                                                        
       ORTHO=.TRUE.                                             
    case ( 'B' )
       BR1_REC(1,1)=PIA(1)                                                   
       BR1_REC(2,2)=PIA(2)                                                   
       BR1_REC(3,3)=PIA(3)                                                   
       !                                                                       
       BR2_REC(1,1)= ZERO                                                    
       BR2_REC(1,2)= PIA(1)                                                  
       BR2_REC(1,3)= PIA(1)                                                  
       BR2_REC(2,1)= PIA(2)                                                  
       BR2_REC(2,2)= ZERO                                                    
       BR2_REC(2,3)= PIA(2)                                                  
       BR2_REC(3,1)= PIA(3)                                                  
       BR2_REC(3,2)= PIA(3)                                                  
       BR2_REC(3,3)= ZERO
       !
       RVFAC=2.D0
       ORTHO=.TRUE.                                             
    case( 'C' )
       select case (lattic(2:3))
       case ( 'XY' )
          BR1_REC(1,1)=PIA(1)
          BR1_REC(2,2)=PIA(2)
          BR1_REC(3,3)=PIA(3)

          BR2_REC(1,1)= PIA(1)
          BR2_REC(1,2)= PIA(1)
          BR2_REC(1,3)= ZERO                                                  
          BR2_REC(2,1)=-PIA(2)
          BR2_REC(2,2)= PIA(2)
          BR2_REC(2,3)= ZERO                                                 
          BR2_REC(3,1)= ZERO                                                  
          BR2_REC(3,2)= ZERO                                                 
          BR2_REC(3,3)= PIA(3)

          RVFAC=2.D0                                                        
          ORTHO=.TRUE.                                             
       case( 'XZ ' )
          if(ABS(ALPHA(3)-PI/2.0D0).LT.0.0001) then
             BR1_REC(1,1)=PIA(1)
             BR1_REC(2,2)=PIA(2)
             BR1_REC(3,3)=PIA(3)

             BR2_REC(1,1)= PIA(1)
             BR2_REC(1,2)= zero                                                   
             BR2_REC(1,3)= PIA(1)
             BR2_REC(2,1)= zero
             BR2_REC(2,2)= PIA(2)
             BR2_REC(2,3)= zero 
             BR2_REC(3,1)=-PIA(3)
             BR2_REC(3,2)= zero   
             BR2_REC(3,3)= PIA(3)
             RVFAC=2.0                                                         
             ORTHO=.TRUE.                                             
          else
             write(6,*) '  gamma not equal 90' !CXZ MONOCLINIC case 
             BR1_REC(1,1)= PIA(1)/SINAB 
             BR1_REC(1,2)= -PIA(2)*COSAB/SINAB
             BR1_REC(2,2)= PIA(2)
             BR1_REC(3,3)= PIA(3)

             BR2_REC(1,1)= PIA(1)/SINAB 
             BR2_REC(1,2)= -PIA(2)*COSAB/SINAB
             BR2_REC(1,3)= PIA(1)/SINAB 
             BR2_REC(2,1)= zero 
             BR2_REC(2,2)= PIA(2)
             BR2_REC(2,3)= zero
             BR2_REC(3,1)=-PIA(3)
             BR2_REC(3,2)= zero
             BR2_REC(3,3)= PIA(3)

             RVFAC=2.0/SINAB                                                   
             ORTHO=.FALSE.                                             
          endif
       case( 'YZ' )
          BR1_REC(1,1)=PIA(1)
          BR1_REC(2,2)=PIA(2)
          BR1_REC(3,3)=PIA(3)
          BR2_REC(1,1)= PIA(1)
          BR2_REC(1,2)= zero
          BR2_REC(1,3)= zero  
          BR2_REC(2,1)= zero  
          BR2_REC(2,2)= PIA(2)
          BR2_REC(2,3)= PIA(2)
          BR2_REC(3,1)= zero  
          BR2_REC(3,2)=-PIA(3)
          BR2_REC(3,3)= PIA(3)
          RVFAC=2.0d0                                                         
          ORTHO=.TRUE.
       end select
    case ( 'R' )
       BR1_REC(1,1)=1.D0/SQRT(3.D0)*PIA(1)
       BR1_REC(1,2)=1.D0/SQRT(3.D0)*PIA(1)
       BR1_REC(1,3)=-2.d0/sqrt(3.d0)*PIA(1)
       BR1_REC(2,1)=-1.0d0*PIA(2)
       BR1_REC(2,2)=1.0d0*PIA(2)
       BR1_REC(2,3)=zero
       BR1_REC(3,1)=1.0d0*PIA(3)
       BR1_REC(3,2)=1.0d0*PIA(3)
       BR1_REC(3,3)=1.0d0*PIA(3)
       
       br2_rec(1:3,1:3)=br1_rec(1:3,1:3)
       
       RVFAC=6.D0/SQRT(3.D0)
       ORTHO=.FALSE.
    case DEFAULT
       stop 'LATGEN - Wrong lattice'
    end select
    vvol=aa*bb*cc/rvfac
    call invert_struct(br1_rec,br1_dir)
    det=0.d0                                                          
    do i=1,3                                                      
       det=det+br1_dir(i,1)*br1_rec(i,1)                                       
    enddo
    br1_dir(1:3,1:3)=br1_dir(1:3,1:3)*2.d0*PI/det
    call invert_struct(br2_rec,br2_dir)
    det=0.d0                                                          
    do i=1,3                                                      
       det=det+br2_dir(i,1)*br2_rec(i,1)                                       
    enddo
    do i=1,3                                                      
       do j=1,3                                                      
          br2_dir(i,j)=br2_dir(i,j)*2.d0*PI/det
       enddo
    enddo
    write(*,'(a)')   '   convention (i,j)     (j)'
    write(*,'(a)')   '                      ax bx cx '
    write(*,'(a)')   '                  (i) ay by cy '
    write(*,'(a,/)') '                      az bz cz '
    write(*,*)  '------------BR1_REC-----------'
    write(*,106)((BR1_rec(I,J),I=1,3),J=1,3)                              
    write(*,*)  '------------BR2_REC-----------'
!br_dir are transposed br_rec**(-1)*2*pi
    write(*,106)((BR2_rec(I,J),I=1,3),J=1,3)                              
    write(*,*)  '------------BR1_DIR-----------'
    write(*,106)((br1_dir(I,J),I=1,3),J=1,3)                              
    write(*,*)  '------------BR2_DIR-----------'
    write(*,106)((br2_dir(I,J),I=1,3),J=1,3)                              
106 format(3(10X,3F10.5,/))                      

! JR: prepare basis vectors for non-primitive cells
  vbr(:,:) = zero
  select case (lattic(1:3))
  case ('P  ')
    nvbr = 1
  case ('S  ')
    nvbr = 1
  case ('H  ')
    nvbr = 1
  case ('B  ')
    nvbr = 2
    vbr(:,2) = (/ 1, 1, 1 /) / two
  case ('F  ')
    nvbr = 4
    vbr(:,2) = (/ 1, 1, 0 /) / four
    vbr(:,3) = (/ 1, 0, 1 /) / four
    vbr(:,4) = (/ 0, 1, 1 /) / four
  case ('R  ')
    nvbr = 3
    vbr(:,2) = (/ 1, 2, 2 /) / three
    vbr(:,3) = (/ 2, 1, 1 /) / three
  case ('CXY')
    nvbr = 2
    vbr(:,2) = (/ 1, 1, 0 /) / two
  case ('CXZ')
    nvbr = 2
    vbr(:,2) = (/ 1, 0, 1 /) / two
  case ('CYZ')
    nvbr = 2
    vbr(:,2) = (/ 0, 1, 1 /) / two
  case default
    write(*,'(a,a)') "ERROR: wrong lattice type: ", lattic(1:3)
  end select

  end subroutine latgen_struct

  subroutine invert_struct(mat1,inv_mat1)
    double precision  :: mat1(3,3),inv_mat1(3,3)
    inv_mat1(1,1)=mat1(2,2)*mat1(3,3)-mat1(3,2)*mat1(2,3)                 
    inv_mat1(2,1)=mat1(3,2)*mat1(1,3)-mat1(1,2)*mat1(3,3)                 
    inv_mat1(3,1)=mat1(1,2)*mat1(2,3)-mat1(2,2)*mat1(1,3)                 
    inv_mat1(1,2)=mat1(2,3)*mat1(3,1)-mat1(3,3)*mat1(2,1)                 
    inv_mat1(2,2)=mat1(3,3)*mat1(1,1)-mat1(1,3)*mat1(3,1)                 
    inv_mat1(3,2)=mat1(1,3)*mat1(2,1)-mat1(2,3)*mat1(1,1)                 
    inv_mat1(1,3)=mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2)                 
    inv_mat1(2,3)=mat1(3,1)*mat1(1,2)-mat1(1,1)*mat1(3,2)                 
    inv_mat1(3,3)=mat1(1,1)*mat1(2,2)-mat1(2,1)*mat1(1,2)                 
  end subroutine invert_struct

end module struct
