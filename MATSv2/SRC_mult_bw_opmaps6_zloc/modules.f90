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
  double precision, parameter     :: TEST=   1.D-10
  double precision, parameter     :: HALF=   5.D-1
  double precision, parameter     :: ZERO=   0.0d0
  double precision, parameter     :: ONE=    1.0d0
  double precision, parameter     :: TWO=    2.0d0
  double precision, parameter     :: THREE=  3.0d0
  double precision, parameter     :: FOUR =  4.0d0
  double precision, parameter     :: NINETY= 90.0d0
  double complex, parameter       :: ZEROC=  (zero,zero)
  double complex, parameter       :: IMAG=   (zero,one)
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
  double precision                 maxfg, maxd0d
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
  integer, pointer              :: qsv(:,:)
end module blochs


module qvectors
! list of q-vectors and corresponding coefficients
  integer                          nqqpr
  double precision, pointer     :: qqpr(:,:)
! multiplicities of q-vectors
  integer, allocatable          :: qvc(:), qvs(:)
  integer                          nqv
end module qvectors


module multislice
! number of cells in x, y, z direction
  integer                          ncelx, ncely, ncelz
  integer                          ms_slices              ! max absolute value of z-component of g-vectors in MS input file
  integer                          shiftx, shifty         ! shift of the lattice
  integer                          meshx, meshy           ! real space mesh within a unit cell
end module multislice


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


module mdffdip

  implicit none

  double precision smom(3), lmom(3), tmom(3), ls, aniorb(3,3), aniso(3,3), nh, lfin

 contains

  subroutine init_mdff

    use inputs, only : lc

    implicit none

    lfin = lc + 1

    smom   = dble(0)
    lmom   = dble(0)
    tmom   = dble(0)
    ls     = dble(0)
    aniorb = dble(0)
    aniso  = dble(0)
    nh     = dble(0)

  end subroutine init_mdff

  subroutine mdff_dip(q1,q2,mdff1,mdff2)

    implicit none

    double precision q1(3), q2(3)
  
    integer          iene, ios, i, j
    double precision qcross(3)
    double precision mdffre1, mdffim1, mdffre2, mdffim2
    double complex   mdff1, mdff2
    double precision l, l2m1, lm1, one, two, three, x, y, sq1, sq2

  ! some useful numbers
    l = dble(lfin)
    l2m1 = dble(2*lfin-1)
    lm1 = dble(l-1)
    one = dble(1)
    two = dble(2)
    three = dble(3)

  ! spin moment is 2*<S> and we will work with <S>
    smom(:) = smom(:)/two

  ! read q,q' and output corresponding S(q,q')
    sq1 = dot_product(q1,q1)
    sq2 = dot_product(q2,q2)
    qcross(1) = q1(2)*q2(3)-q1(3)*q2(2)
    qcross(2) = q1(3)*q2(1)-q1(1)*q2(3)
    qcross(3) = q1(1)*q2(2)-q1(2)*q2(1)
    mdffim1 = dot_product(qcross,lmom* l /two + (smom*l+tmom*(two*l+three))*lm1/three)
    mdffim2 = dot_product(qcross,lmom*lm1/two - (smom*l+tmom*(two*l+three))*lm1/three)
    x = ( dot_product(q1,matmul(aniorb,q2)) + dot_product(q2,matmul(aniorb,q1)) ) / two + dot_product(q1,q2)*two*l*l2m1*nh/three
    y = ( dot_product(q1,matmul(aniso,q2))  + dot_product(q2,matmul(aniso,q1))  ) / two / l - dot_product(q1,q2)*ls*l
    mdffre1 =  l *x/l2m1 + l*y/l2m1
    mdffre2 = lm1*x/l2m1 - l*y/l2m1
    mdff1 = dcmplx(mdffre1,mdffim1)/(sq1*sq2)
    mdff2 = dcmplx(mdffre2,mdffim2)/(sq1*sq2)

  end subroutine mdff_dip

end module mdffdip


module asff 
! atom scattering formfactors

  integer, parameter     :: maxatn = 120
  double precision          fa(4,maxatn), fb(4,maxatn)
  character*2               atnam(maxatn)
  integer                   isour(maxatn)

 contains 

  subroutine init_asff

    implicit none

    integer i,j

! From atoms.sff converted by Jan Rusz, 2012
! line 1: Ac, at. no. 89
    isour(89) = 0
    atnam(89) = "Ac"
    fa(1,89) = 6.278/10.0d0
    fb(1,89) = 28.323/100.0d0
    fa(2,89) = 5.195/10.0d0
    fb(2,89) = 4.949/100.0d0
    fa(3,89) = 2.321/10.0d0
    fb(3,89) = 0.557/100.0d0
    fa(4,89) = 0.000/10.0d0
    fb(4,89) = 0.000/100.0d0
! line 2: Ag, at. no. 47
    isour(47) = 1
    atnam(47) = "Ag"
    fa(1,47) = 2.036/10.0d0
    fb(1,47) = 61.497/100.0d0
    fa(2,47) = 3.272/10.0d0
    fb(2,47) = 11.824/100.0d0
    fa(3,47) = 2.511/10.0d0
    fb(3,47) = 2.846/100.0d0
    fa(4,47) = 0.837/10.0d0
    fb(4,47) = 0.327/100.0d0
! line 3: Al, at. no. 13
    isour(13) = 1
    atnam(13) = "Al"
    fa(1,13) = 2.276/10.0d0
    fb(1,13) = 72.322/100.0d0
    fa(2,13) = 2.428/10.0d0
    fb(2,13) = 19.773/100.0d0
    fa(3,13) = 0.858/10.0d0
    fb(3,13) = 3.080/100.0d0
    fa(4,13) = 0.317/10.0d0
    fb(4,13) = 0.408/100.0d0
! line 4: Am, at. no. 95
    isour(95) = 0
    atnam(95) = "Am"
    fa(1,95) = 6.378/10.0d0
    fb(1,95) = 29.156/100.0d0
    fa(2,95) = 5.495/10.0d0
    fb(2,95) = 5.102/100.0d0
    fa(3,95) = 2.495/10.0d0
    fb(3,95) = 0.565/100.0d0
    fa(4,95) = 0.000/10.0d0
    fb(4,95) = 0.000/100.0d0
! line 5: Ar, at. no. 18
    isour(18) = 1
    atnam(18) = "Ar"
    fa(1,18) = 1.274/10.0d0
    fb(1,18) = 26.682/100.0d0
    fa(2,18) = 2.190/10.0d0
    fb(2,18) = 8.813/100.0d0
    fa(3,18) = 0.793/10.0d0
    fb(3,18) = 2.219/100.0d0
    fa(4,18) = 0.326/10.0d0
    fb(4,18) = 0.307/100.0d0
! line 6: As, at. no. 33
    isour(33) = 1
    atnam(33) = "As"
    fa(1,33) = 2.399/10.0d0
    fb(1,33) = 45.718/100.0d0
    fa(2,33) = 2.790/10.0d0
    fb(2,33) = 12.817/100.0d0
    fa(3,33) = 1.529/10.0d0
    fb(3,33) = 2.280/100.0d0
    fa(4,33) = 0.594/10.0d0
    fb(4,33) = 0.328/100.0d0
! line 7: At, at. no. 85
    isour(85) = 0
    atnam(85) = "At"
    fa(1,85) = 6.133/10.0d0
    fb(1,85) = 28.047/100.0d0
    fa(2,85) = 5.031/10.0d0
    fb(2,85) = 4.957/100.0d0
    fa(3,85) = 2.239/10.0d0
    fb(3,85) = 0.558/100.0d0
    fa(4,85) = 0.000/10.0d0
    fb(4,85) = 0.000/100.0d0
! line 8: Au, at. no. 79
    isour(79) = 1
    atnam(79) = "Au"
    fa(1,79) = 2.388/10.0d0
    fb(1,79) = 42.866/100.0d0
    fa(2,79) = 4.226/10.0d0
    fb(2,79) = 9.743/100.0d0
    fa(3,79) = 2.689/10.0d0
    fb(3,79) = 2.264/100.0d0
    fa(4,79) = 1.255/10.0d0
    fb(4,79) = 0.307/100.0d0
! line 9: B, at. no. 05
    isour(05) = 1
    atnam(05) = "B"
    fa(1,05) = 0.945/10.0d0
    fb(1,05) = 46.444/100.0d0
    fa(2,05) = 1.312/10.0d0
    fb(2,05) = 14.178/100.0d0
    fa(3,05) = 0.419/10.0d0
    fb(3,05) = 3.223/100.0d0
    fa(4,05) = 0.116/10.0d0
    fb(4,05) = 0.377/100.0d0
! line 10: Ba, at. no. 56
    isour(56) = 1
    atnam(56) = "Ba"
    fa(1,56) = 7.821/10.0d0
    fb(1,56) = 117.657/100.0d0
    fa(2,56) = 6.004/10.0d0
    fb(2,56) = 18.778/100.0d0
    fa(3,56) = 3.280/10.0d0
    fb(3,56) = 3.263/100.0d0
    fa(4,56) = 1.103/10.0d0
    fb(4,56) = 0.376/100.0d0
! line 11: Be, at. no. 04
    isour(04) = 1
    atnam(04) = "Be"
    fa(1,04) = 1.250/10.0d0
    fb(1,04) = 60.804/100.0d0
    fa(2,04) = 1.334/10.0d0
    fb(2,04) = 18.591/100.0d0
    fa(3,04) = 0.360/10.0d0
    fb(3,04) = 3.653/100.0d0
    fa(4,04) = 0.106/10.0d0
    fb(4,04) = 0.416/100.0d0
! line 12: Bi, at. no. 83
    isour(83) = 1
    atnam(83) = "Bi"
    fa(1,83) = 3.841/10.0d0
    fb(1,83) = 50.261/100.0d0
    fa(2,83) = 4.679/10.0d0
    fb(2,83) = 11.999/100.0d0
    fa(3,83) = 3.192/10.0d0
    fb(3,83) = 2.560/100.0d0
    fa(4,83) = 1.363/10.0d0
    fb(4,83) = 0.318/100.0d0
! line 13: Bk, at. no. 97
    isour(97) = 0
    atnam(97) = "Bk"
    fa(1,97) = 6.502/10.0d0
    fb(1,97) = 28.375/100.0d0
    fa(2,97) = 5.478/10.0d0
    fb(2,97) = 4.975/100.0d0
    fa(3,97) = 2.510/10.0d0
    fb(3,97) = 0.561/100.0d0
    fa(4,97) = 0.000/10.0d0
    fb(4,97) = 0.000/100.0d0
! line 14: Br, at. no. 35
    isour(35) = 1
    atnam(35) = "Br"
    fa(1,35) = 2.166/10.0d0
    fb(1,35) = 33.899/100.0d0
    fa(2,35) = 2.904/10.0d0
    fb(2,35) = 10.497/100.0d0
    fa(3,35) = 1.395/10.0d0
    fb(3,35) = 2.041/100.0d0
    fa(4,35) = 0.589/10.0d0
    fb(4,35) = 0.307/100.0d0
! line 15: C, at. no. 06
    isour(06) = 1
    atnam(06) = "C"
    fa(1,06) = 0.731/10.0d0
    fb(1,06) = 36.995/100.0d0
    fa(2,06) = 1.195/10.0d0
    fb(2,06) = 11.297/100.0d0
    fa(3,06) = 0.456/10.0d0
    fb(3,06) = 2.814/100.0d0
    fa(4,06) = 0.125/10.0d0
    fb(4,06) = 0.346/100.0d0
! line 16: Ca, at. no. 20
    isour(20) = 1
    atnam(20) = "Ca"
    fa(1,20) = 4.470/10.0d0
    fb(1,20) = 99.523/100.0d0
    fa(2,20) = 2.971/10.0d0
    fb(2,20) = 22.696/100.0d0
    fa(3,20) = 1.970/10.0d0
    fb(3,20) = 4.195/100.0d0
    fa(4,20) = 0.482/10.0d0
    fb(4,20) = 0.417/100.0d0
! line 17: Cd, at. no. 48
    isour(48) = 1
    atnam(48) = "Cd"
    fa(1,48) = 2.574/10.0d0
    fb(1,48) = 55.675/100.0d0
    fa(2,48) = 3.259/10.0d0
    fb(2,48) = 11.838/100.0d0
    fa(3,48) = 2.547/10.0d0
    fb(3,48) = 2.784/100.0d0
    fa(4,48) = 0.838/10.0d0
    fb(4,48) = 0.322/100.0d0
! line 18: Ce, at. no. 58
    isour(58) = 0
    atnam(58) = "Ce"
    fa(1,58) = 5.007/10.0d0
    fb(1,58) = 28.283/100.0d0
    fa(2,58) = 3.980/10.0d0
    fb(2,58) = 5.183/100.0d0
    fa(3,58) = 1.678/10.0d0
    fb(3,58) = 0.589/100.0d0
    fa(4,58) = 0.000/10.0d0
    fb(4,58) = 0.000/100.0d0
! line 19: Cf, at. no. 98
    isour(98) = 0
    atnam(98) = "Cf"
    fa(1,98) = 6.548/10.0d0
    fb(1,98) = 28.461/100.0d0
    fa(2,98) = 5.526/10.0d0
    fb(2,98) = 4.965/100.0d0
    fa(3,98) = 2.520/10.0d0
    fb(3,98) = 0.557/100.0d0
    fa(4,98) = 0.000/10.0d0
    fb(4,98) = 0.000/100.0d0
! line 20: Cl, at. no. 17
    isour(17) = 1
    atnam(17) = "Cl"
    fa(1,17) = 1.452/10.0d0
    fb(1,17) = 30.935/100.0d0
    fa(2,17) = 2.292/10.0d0
    fb(2,17) = 9.980/100.0d0
    fa(3,17) = 0.787/10.0d0
    fb(3,17) = 2.234/100.0d0
    fa(4,17) = 0.322/10.0d0
    fb(4,17) = 0.323/100.0d0
! line 21: Cm, at. no. 96
    isour(96) = 0
    atnam(96) = "Cm"
    fa(1,96) = 6.460/10.0d0
    fb(1,96) = 28.396/100.0d0
    fa(2,96) = 5.469/10.0d0
    fb(2,96) = 4.970/100.0d0
    fa(3,96) = 2.471/10.0d0
    fb(3,96) = 0.554/100.0d0
    fa(4,96) = 0.000/10.0d0
    fb(4,96) = 0.000/100.0d0
! line 22: Co, at. no. 27
    isour(27) = 1
    atnam(27) = "Co"
    fa(1,27) = 2.367/10.0d0
    fb(1,27) = 61.431/100.0d0
    fa(2,27) = 2.236/10.0d0
    fb(2,27) = 14.180/100.0d0
    fa(3,27) = 1.724/10.0d0
    fb(3,27) = 2.725/100.0d0
    fa(4,27) = 0.515/10.0d0
    fb(4,27) = 0.344/100.0d0
! line 23: Cr, at. no. 24
    isour(24) = 1
    atnam(24) = "Cr"
    fa(1,24) = 2.307/10.0d0
    fb(1,24) = 78.405/100.0d0
    fa(2,24) = 2.334/10.0d0
    fb(2,24) = 15.785/100.0d0
    fa(3,24) = 1.823/10.0d0
    fb(3,24) = 3.157/100.0d0
    fa(4,24) = 0.490/10.0d0
    fb(4,24) = 0.364/100.0d0
! line 24: Cs, at. no. 55
    isour(55) = 1
    atnam(55) = "Cs"
    fa(1,55) = 6.062/10.0d0
    fb(1,55) = 155.837/100.0d0
    fa(2,55) = 5.986/10.0d0
    fb(2,55) = 19.695/100.0d0
    fa(3,55) = 3.303/10.0d0
    fb(3,55) = 3.335/100.0d0
    fa(4,55) = 1.096/10.0d0
    fb(4,55) = 0.379/100.0d0
! line 25: Cu, at. no. 29
    isour(29) = 1
    atnam(29) = "Cu"
    fa(1,29) = 1.579/10.0d0
    fb(1,29) = 62.940/100.0d0
    fa(2,29) = 1.820/10.0d0
    fb(2,29) = 12.453/100.0d0
    fa(3,29) = 1.658/10.0d0
    fb(3,29) = 2.504/100.0d0
    fa(4,29) = 0.532/10.0d0
    fb(4,29) = 0.333/100.0d0
! line 26: Dy, at. no. 66
    isour(66) = 0
    atnam(66) = "Dy"
    fa(1,66) = 5.332/10.0d0
    fb(1,66) = 28.888/100.0d0
    fa(2,66) = 4.370/10.0d0
    fb(2,66) = 5.198/100.0d0
    fa(3,66) = 1.863/10.0d0
    fb(3,66) = 0.581/100.0d0
    fa(4,66) = 0.000/10.0d0
    fb(4,66) = 0.000/100.0d0
! line 27: Er, at. no. 68
    isour(68) = 0
    atnam(68) = "Er"
    fa(1,68) = 5.436/10.0d0
    fb(1,68) = 28.655/100.0d0
    fa(2,68) = 4.437/10.0d0
    fb(2,68) = 5.117/100.0d0
    fa(3,68) = 1.891/10.0d0
    fb(3,68) = 0.577/100.0d0
    fa(4,68) = 0.000/10.0d0
    fb(4,68) = 0.000/100.0d0
! line 28: Eu, at. no. 63
    isour(63) = 1
    atnam(63) = "Eu"
    fa(1,63) = 6.267/10.0d0
    fb(1,63) = 100.298/100.0d0
    fa(2,63) = 4.844/10.0d0
    fb(2,63) = 16.066/100.0d0
    fa(3,63) = 3.202/10.0d0
    fb(3,63) = 2.980/100.0d0
    fa(4,63) = 1.200/10.0d0
    fb(4,63) = 0.367/100.0d0
! line 29: F, at. no. 09
    isour(09) = 1
    atnam(09) = "F"
    fa(1,09) = 0.387/10.0d0
    fb(1,09) = 20.239/100.0d0
    fa(2,09) = 0.811/10.0d0
    fb(2,09) = 6.609/100.0d0
    fa(3,09) = 0.475/10.0d0
    fb(3,09) = 1.931/100.0d0
    fa(4,09) = 0.146/10.0d0
    fb(4,09) = 0.279/100.0d0
! line 30: Fe, at. no. 26
    isour(26) = 1
    atnam(26) = "Fe"
    fa(1,26) = 2.544/10.0d0
    fb(1,26) = 64.424/100.0d0
    fa(2,26) = 2.343/10.0d0
    fb(2,26) = 14.880/100.0d0
    fa(3,26) = 1.759/10.0d0
    fb(3,26) = 2.854/100.0d0
    fa(4,26) = 0.506/10.0d0
    fb(4,26) = 0.350/100.0d0
! line 31: Fr, at. no. 87
    isour(87) = 0
    atnam(87) = "Fr"
    fa(1,87) = 6.201/10.0d0
    fb(1,87) = 28.200/100.0d0
    fa(2,87) = 5.121/10.0d0
    fb(2,87) = 4.954/100.0d0
    fa(3,87) = 2.275/10.0d0
    fb(3,87) = 0.556/100.0d0
    fa(4,87) = 0.000/10.0d0
    fb(4,87) = 0.000/100.0d0
! line 32: Ga, at. no. 31
    isour(31) = 1
    atnam(31) = "Ga"
    fa(1,31) = 2.321/10.0d0
    fb(1,31) = 65.602/100.0d0
    fa(2,31) = 2.486/10.0d0
    fb(2,31) = 15.458/100.0d0
    fa(3,31) = 1.688/10.0d0
    fb(3,31) = 2.581/100.0d0
    fa(4,31) = 0.599/10.0d0
    fb(4,31) = 0.351/100.0d0
! line 33: Gd, at. no. 64
    isour(64) = 0
    atnam(64) = "Gd"
    fa(1,64) = 5.225/10.0d0
    fb(1,64) = 29.158/100.0d0
    fa(2,64) = 4.314/10.0d0
    fb(2,64) = 5.259/100.0d0
    fa(3,64) = 1.827/10.0d0
    fb(3,64) = 0.586/100.0d0
    fa(4,64) = 0.000/10.0d0
    fb(4,64) = 0.000/100.0d0
! line 34: Ge, at. no. 32
    isour(32) = 1
    atnam(32) = "Ge"
    fa(1,32) = 2.447/10.0d0
    fb(1,32) = 55.893/100.0d0
    fa(2,32) = 2.702/10.0d0
    fb(2,32) = 14.393/100.0d0
    fa(3,32) = 1.616/10.0d0
    fb(3,32) = 2.446/100.0d0
    fa(4,32) = 0.601/10.0d0
    fb(4,32) = 0.342/100.0d0
! line 35: H, at. no. 01
    isour(01) = 0
    atnam(01) = "H"
    fa(1,01) = 0.202/10.0d0
    fb(1,01) = 30.868/100.0d0
    fa(2,01) = 0.244/10.0d0
    fb(2,01) = 8.544/100.0d0
    fa(3,01) = 0.082/10.0d0
    fb(3,01) = 1.273/100.0d0
    fa(4,01) = 0.000/10.0d0
    fb(4,01) = 0.000/100.0d0
! line 36: He, at. no. 02
    isour(02) = 1
    atnam(02) = "He"
    fa(1,02) = 0.091/10.0d0
    fb(1,02) = 18.183/100.0d0
    fa(2,02) = 0.181/10.0d0
    fb(2,02) = 6.212/100.0d0
    fa(3,02) = 0.110/10.0d0
    fb(3,02) = 1.803/100.0d0
    fa(4,02) = 0.036/10.0d0
    fb(4,02) = 0.284/100.0d0
! line 37: Hf, at. no. 72
    isour(72) = 0
    atnam(72) = "Hf"
    fa(1,72) = 5.588/10.0d0
    fb(1,72) = 29.001/100.0d0
    fa(2,72) = 4.619/10.0d0
    fb(2,72) = 5.164/100.0d0
    fa(3,72) = 1.997/10.0d0
    fb(3,72) = 0.579/100.0d0
    fa(4,72) = 0.000/10.0d0
    fb(4,72) = 0.000/100.0d0
! line 38: Hg, at. no. 80
    isour(80) = 1
    atnam(80) = "Hg"
    fa(1,80) = 2.682/10.0d0
    fb(1,80) = 42.822/100.0d0
    fa(2,80) = 4.241/10.0d0
    fb(2,80) = 9.856/100.0d0
    fa(3,80) = 2.755/10.0d0
    fb(3,80) = 2.295/100.0d0
    fa(4,80) = 1.270/10.0d0
    fb(4,80) = 0.307/100.0d0
! line 39: Ho, at. no. 67
    isour(67) = 0
    atnam(67) = "Ho"
    fa(1,67) = 5.376/10.0d0
    fb(1,67) = 28.773/100.0d0
    fa(2,67) = 4.403/10.0d0
    fb(2,67) = 5.174/100.0d0
    fa(3,67) = 1.884/10.0d0
    fb(3,67) = 0.582/100.0d0
    fa(4,67) = 0.000/10.0d0
    fb(4,67) = 0.000/100.0d0
! line 40: I, at. no. 53
    isour(53) = 1
    atnam(53) = "I"
    fa(1,53) = 3.473/10.0d0
    fb(1,53) = 39.441/100.0d0
    fa(2,53) = 4.060/10.0d0
    fb(2,53) = 11.816/100.0d0
    fa(3,53) = 2.522/10.0d0
    fb(3,53) = 2.415/100.0d0
    fa(4,53) = 0.840/10.0d0
    fb(4,53) = 0.298/100.0d0
! line 41: In, at. no. 49
    isour(49) = 1
    atnam(49) = "In"
    fa(1,49) = 3.153/10.0d0
    fb(1,49) = 66.649/100.0d0
    fa(2,49) = 3.557/10.0d0
    fb(2,49) = 14.449/100.0d0
    fa(3,49) = 2.818/10.0d0
    fb(3,49) = 2.976/100.0d0
    fa(4,49) = 0.884/10.0d0
    fb(4,49) = 0.335/100.0d0
! line 42: Ir, at. no. 77
    isour(77) = 0
    atnam(77) = "Ir"
    fa(1,77) = 5.754/10.0d0
    fb(1,77) = 29.159/100.0d0
    fa(2,77) = 4.851/10.0d0
    fb(2,77) = 5.152/100.0d0
    fa(3,77) = 2.096/10.0d0
    fb(3,77) = 0.570/100.0d0
    fa(4,77) = 0.000/10.0d0
    fb(4,77) = 0.000/100.0d0
! line 43: K, at. no. 19
    isour(19) = 1
    atnam(19) = "K"
    fa(1,19) = 3.951/10.0d0
    fb(1,19) = 137.075/100.0d0
    fa(2,19) = 2.545/10.0d0
    fb(2,19) = 22.402/100.0d0
    fa(3,19) = 1.980/10.0d0
    fb(3,19) = 4.532/100.0d0
    fa(4,19) = 0.482/10.0d0
    fb(4,19) = 0.434/100.0d0
! line 44: Kr, at. no. 36
    isour(36) = 1
    atnam(36) = "Kr"
    fa(1,36) = 2.034/10.0d0
    fb(1,36) = 29.999/100.0d0
    fa(2,36) = 2.927/10.0d0
    fb(2,36) = 9.598/100.0d0
    fa(3,36) = 1.342/10.0d0
    fb(3,36) = 1.952/100.0d0
    fa(4,36) = 0.589/10.0d0
    fb(4,36) = 0.299/100.0d0
! line 45: La, at. no. 57
    isour(57) = 0
    atnam(57) = "La"
    fa(1,57) = 4.940/10.0d0
    fb(1,57) = 28.716/100.0d0
    fa(2,57) = 3.968/10.0d0
    fb(2,57) = 5.245/100.0d0
    fa(3,57) = 1.663/10.0d0
    fb(3,57) = 0.594/100.0d0
    fa(4,57) = 0.000/10.0d0
    fb(4,57) = 0.000/100.0d0
! line 46: Li, at. no. 03
    isour(03) = 1
    atnam(03) = "Li"
    fa(1,03) = 1.611/10.0d0
    fb(1,03) = 107.638/100.0d0
    fa(2,03) = 1.246/10.0d0
    fb(2,03) = 30.480/100.0d0
    fa(3,03) = 0.326/10.0d0
    fb(3,03) = 4.533/100.0d0
    fa(4,03) = 0.099/10.0d0
    fb(4,03) = 0.495/100.0d0
! line 47: Lu, at. no. 71
    isour(71) = 0
    atnam(71) = "Lu"
    fa(1,71) = 5.553/10.0d0
    fb(1,71) = 28.907/100.0d0
    fa(2,71) = 4.580/10.0d0
    fb(2,71) = 5.160/100.0d0
    fa(3,71) = 1.969/10.0d0
    fb(3,71) = 0.577/100.0d0
    fa(4,71) = 0.000/10.0d0
    fb(4,71) = 0.000/100.0d0
! line 48: Mg, at. no. 12
    isour(12) = 1
    atnam(12) = "Mg"
    fa(1,12) = 2.268/10.0d0
    fb(1,12) = 73.670/100.0d0
    fa(2,12) = 1.803/10.0d0
    fb(2,12) = 20.175/100.0d0
    fa(3,12) = 0.839/10.0d0
    fb(3,12) = 3.013/100.0d0
    fa(4,12) = 0.289/10.0d0
    fb(4,12) = 0.405/100.0d0
! line 49: Mn, at. no. 25
    isour(25) = 1
    atnam(25) = "Mn"
    fa(1,25) = 2.747/10.0d0
    fb(1,25) = 67.786/100.0d0
    fa(2,25) = 2.456/10.0d0
    fb(2,25) = 15.674/100.0d0
    fa(3,25) = 1.792/10.0d0
    fb(3,25) = 3.000/100.0d0
    fa(4,25) = 0.498/10.0d0
    fb(4,25) = 0.357/100.0d0
! line 50: Mo, at. no. 42
    isour(42) = 1
    atnam(42) = "Mo"
    fa(1,42) = 3.120/10.0d0
    fb(1,42) = 72.464/100.0d0
    fa(2,42) = 3.906/10.0d0
    fb(2,42) = 14.642/100.0d0
    fa(3,42) = 2.361/10.0d0
    fb(3,42) = 3.237/100.0d0
    fa(4,42) = 0.850/10.0d0
    fb(4,42) = 0.366/100.0d0
! line 51: N, at. no. 07
    isour(07) = 1
    atnam(07) = "N"
    fa(1,07) = 0.572/10.0d0
    fb(1,07) = 28.847/100.0d0
    fa(2,07) = 1.043/10.0d0
    fb(2,07) = 9.054/100.0d0
    fa(3,07) = 0.465/10.0d0
    fb(3,07) = 2.421/100.0d0
    fa(4,07) = 0.131/10.0d0
    fb(4,07) = 0.317/100.0d0
! line 52: Na, at. no. 11
    isour(11) = 1
    atnam(11) = "Na"
    fa(1,11) = 2.241/10.0d0
    fb(1,11) = 108.004/100.0d0
    fa(2,11) = 1.333/10.0d0
    fb(2,11) = 24.505/100.0d0
    fa(3,11) = 0.907/10.0d0
    fb(3,11) = 3.391/100.0d0
    fa(4,11) = 0.286/10.0d0
    fb(4,11) = 0.435/100.0d0
! line 53: Nb, at. no. 41
    isour(41) = 0
    atnam(41) = "Nb"
    fa(1,41) = 4.237/10.0d0
    fb(1,41) = 27.415/100.0d0
    fa(2,41) = 3.105/10.0d0
    fb(2,41) = 5.074/100.0d0
    fa(3,41) = 1.234/10.0d0
    fb(3,41) = 0.593/100.0d0
    fa(4,41) = 0.000/10.0d0
    fb(4,41) = 0.000/100.0d0
! line 54: Nd, at. no. 60
    isour(60) = 0
    atnam(60) = "Nd"
    fa(1,60) = 5.151/10.0d0
    fb(1,60) = 28.304/100.0d0
    fa(2,60) = 4.075/10.0d0
    fb(2,60) = 5.073/100.0d0
    fa(3,60) = 1.683/10.0d0
    fb(3,60) = 0.571/100.0d0
    fa(4,60) = 0.000/10.0d0
    fb(4,60) = 0.000/100.0d0
! line 55: Ne, at. no. 10
    isour(10) = 1
    atnam(10) = "Ne"
    fa(1,10) = 0.303/10.0d0
    fb(1,10) = 17.640/100.0d0
    fa(2,10) = 0.720/10.0d0
    fb(2,10) = 5.860/100.0d0
    fa(3,10) = 0.475/10.0d0
    fb(3,10) = 1.762/100.0d0
    fa(4,10) = 0.153/10.0d0
    fb(4,10) = 0.266/100.0d0
! line 56: Ni, at. no. 28
    isour(28) = 1
    atnam(28) = "Ni"
    fa(1,28) = 2.210/10.0d0
    fb(1,28) = 58.727/100.0d0
    fa(2,28) = 2.134/10.0d0
    fb(2,28) = 13.553/100.0d0
    fa(3,28) = 1.689/10.0d0
    fb(3,28) = 2.609/100.0d0
    fa(4,28) = 0.524/10.0d0
    fb(4,28) = 0.339/100.0d0
! line 57: Np, at. no. 93
    isour(93) = 0
    atnam(93) = "Np"
    fa(1,93) = 6.323/10.0d0
    fb(1,93) = 29.142/100.0d0
    fa(2,93) = 5.414/10.0d0
    fb(2,93) = 5.096/100.0d0
    fa(3,93) = 2.453/10.0d0
    fb(3,93) = 0.568/100.0d0
    fa(4,93) = 0.000/10.0d0
    fb(4,93) = 0.000/100.0d0
! line 58: O, at. no. 08
    isour(08) = 1
    atnam(08) = "O"
    fa(1,08) = 0.455/10.0d0
    fb(1,08) = 23.780/100.0d0
    fa(2,08) = 0.917/10.0d0
    fb(2,08) = 7.622/100.0d0
    fa(3,08) = 0.472/10.0d0
    fb(3,08) = 2.144/100.0d0
    fa(4,08) = 0.138/10.0d0
    fb(4,08) = 0.296/100.0d0
! line 59: Os, at. no. 76
    isour(76) = 0
    atnam(76) = "Os"
    fa(1,76) = 5.750/10.0d0
    fb(1,76) = 28.933/100.0d0
    fa(2,76) = 4.773/10.0d0
    fb(2,76) = 5.139/100.0d0
    fa(3,76) = 2.079/10.0d0
    fb(3,76) = 0.573/100.0d0
    fa(4,76) = 0.000/10.0d0
    fb(4,76) = 0.000/100.0d0
! line 60: P, at. no. 15
    isour(15) = 1
    atnam(15) = "P"
    fa(1,15) = 1.888/10.0d0
    fb(1,15) = 44.876/100.0d0
    fa(2,15) = 2.469/10.0d0
    fb(2,15) = 13.538/100.0d0
    fa(3,15) = 0.805/10.0d0
    fb(3,15) = 2.642/100.0d0
    fa(4,15) = 0.320/10.0d0
    fb(4,15) = 0.361/100.0d0
! line 61: Pa, at. no. 91
    isour(91) = 0
    atnam(91) = "Pa"
    fa(1,91) = 6.306/10.0d0
    fb(1,91) = 28.688/100.0d0
    fa(2,91) = 5.303/10.0d0
    fb(2,91) = 5.026/100.0d0
    fa(3,91) = 2.386/10.0d0
    fb(3,91) = 0.561/100.0d0
    fa(4,91) = 0.000/10.0d0
    fb(4,91) = 0.000/100.0d0
! line 62: Pb, at. no. 82
    isour(82) = 1
    atnam(82) = "Pb"
    fa(1,82) = 3.510/10.0d0
    fb(1,82) = 52.914/100.0d0
    fa(2,82) = 4.552/10.0d0
    fb(2,82) = 11.884/100.0d0
    fa(3,82) = 3.154/10.0d0
    fb(3,82) = 2.571/100.0d0
    fa(4,82) = 1.359/10.0d0
    fb(4,82) = 0.321/100.0d0
! line 63: Pd, at. no. 46
    isour(46) = 0
    atnam(46) = "Pd"
    fa(1,46) = 4.436/10.0d0
    fb(1,46) = 28.670/100.0d0
    fa(2,46) = 3.454/10.0d0
    fb(2,46) = 5.269/100.0d0
    fa(3,46) = 1.383/10.0d0
    fb(3,46) = 0.595/100.0d0
    fa(4,46) = 0.000/10.0d0
    fb(4,46) = 0.000/100.0d0
! line 64: Pm, at. no. 61
    isour(61) = 0
    atnam(61) = "Pm"
    fa(1,61) = 5.201/10.0d0
    fb(1,61) = 28.079/100.0d0
    fa(2,61) = 4.094/10.0d0
    fb(2,61) = 5.081/100.0d0
    fa(3,61) = 1.719/10.0d0
    fb(3,61) = 0.576/100.0d0
    fa(4,61) = 0.000/10.0d0
    fb(4,61) = 0.000/100.0d0
! line 65: Po, at. no. 84
    isour(84) = 0
    atnam(84) = "Po"
    fa(1,84) = 6.070/10.0d0
    fb(1,84) = 28.075/100.0d0
    fa(2,84) = 4.997/10.0d0
    fb(2,84) = 4.999/100.0d0
    fa(3,84) = 2.232/10.0d0
    fb(3,84) = 0.563/100.0d0
    fa(4,84) = 0.000/10.0d0
    fb(4,84) = 0.000/100.0d0
! line 66: Pr, at. no. 59
    isour(59) = 0
    atnam(59) = "Pr"
    fa(1,59) = 5.085/10.0d0
    fb(1,59) = 28.588/100.0d0
    fa(2,59) = 4.043/10.0d0
    fb(2,59) = 5.143/100.0d0
    fa(3,59) = 1.684/10.0d0
    fb(3,59) = 0.581/100.0d0
    fa(4,59) = 0.000/10.0d0
    fb(4,59) = 0.000/100.0d0
! line 67: Pt, at. no. 78
    isour(78) = 0
    atnam(78) = "Pt"
    fa(1,78) = 5.803/10.0d0
    fb(1,78) = 29.016/100.0d0
    fa(2,78) = 4.870/10.0d0
    fb(2,78) = 5.150/100.0d0
    fa(3,78) = 2.127/10.0d0
    fb(3,78) = 0.572/100.0d0
    fa(4,78) = 0.000/10.0d0
    fb(4,78) = 0.000/100.0d0
! line 68: Pu, at. no. 94
    isour(94) = 0
    atnam(94) = "Pu"
    fa(1,94) = 6.415/10.0d0
    fb(1,94) = 28.836/100.0d0
    fa(2,94) = 5.419/10.0d0
    fb(2,94) = 5.022/100.0d0
    fa(3,94) = 2.449/10.0d0
    fb(3,94) = 0.561/100.0d0
    fa(4,94) = 0.000/10.0d0
    fb(4,94) = 0.000/100.0d0
! line 69: Ra, at. no. 88
    isour(88) = 0
    atnam(88) = "Ra"
    fa(1,88) = 6.215/10.0d0
    fb(1,88) = 28.382/100.0d0
    fa(2,88) = 5.170/10.0d0
    fb(2,88) = 5.002/100.0d0
    fa(3,88) = 2.316/10.0d0
    fb(3,88) = 0.562/100.0d0
    fa(4,88) = 0.000/10.0d0
    fb(4,88) = 0.000/100.0d0
! line 70: Rb, at. no. 37
    isour(37) = 1
    atnam(37) = "Rb"
    fa(1,37) = 4.776/10.0d0
    fb(1,37) = 140.782/100.0d0
    fa(2,37) = 3.859/10.0d0
    fb(2,37) = 18.991/100.0d0
    fa(3,37) = 2.234/10.0d0
    fb(3,37) = 3.701/100.0d0
    fa(4,37) = 0.868/10.0d0
    fb(4,37) = 0.419/100.0d0
! line 71: Re, at. no. 75
    isour(75) = 0
    atnam(75) = "Re"
    fa(1,75) = 5.695/10.0d0
    fb(1,75) = 28.968/100.0d0
    fa(2,75) = 4.740/10.0d0
    fb(2,75) = 5.156/100.0d0
    fa(3,75) = 2.064/10.0d0
    fb(3,75) = 0.575/100.0d0
    fa(4,75) = 0.000/10.0d0
    fb(4,75) = 0.000/100.0d0
! line 72: Rh, at. no. 45
    isour(45) = 0
    atnam(45) = "Rh"
    fa(1,45) = 4.431/10.0d0
    fb(1,45) = 27.911/100.0d0
    fa(2,45) = 3.343/10.0d0
    fb(2,45) = 5.153/100.0d0
    fa(3,45) = 1.345/10.0d0
    fb(3,45) = 0.592/100.0d0
    fa(4,45) = 0.000/10.0d0
    fb(4,45) = 0.000/100.0d0
! line 73: Rn, at. no. 86
    isour(86) = 1
    atnam(86) = "Rn"
    fa(1,86) = 4.078/10.0d0
    fb(1,86) = 38.406/100.0d0
    fa(2,86) = 4.978/10.0d0
    fb(2,86) = 11.020/100.0d0
    fa(3,86) = 3.096/10.0d0
    fb(3,86) = 2.355/100.0d0
    fa(4,86) = 1.326/10.0d0
    fb(4,86) = 0.299/100.0d0
! line 74: Ru, at. no. 44
    isour(44) = 0
    atnam(44) = "Ru"
    fa(1,44) = 4.358/10.0d0
    fb(1,44) = 27.881/100.0d0
    fa(2,44) = 3.298/10.0d0
    fb(2,44) = 5.179/100.0d0
    fa(3,44) = 1.323/10.0d0
    fb(3,44) = 0.594/100.0d0
    fa(4,44) = 0.000/10.0d0
    fb(4,44) = 0.000/100.0d0
! line 75: S, at. no. 16
    isour(16) = 1
    atnam(16) = "S"
    fa(1,16) = 1.659/10.0d0
    fb(1,16) = 36.650/100.0d0
    fa(2,16) = 2.386/10.0d0
    fb(2,16) = 11.488/100.0d0
    fa(3,16) = 0.790/10.0d0
    fb(3,16) = 2.469/100.0d0
    fa(4,16) = 0.321/10.0d0
    fb(4,16) = 0.340/100.0d0
! line 76: Sb, at. no. 51
    isour(51) = 1
    atnam(51) = "Sb"
    fa(1,51) = 3.564/10.0d0
    fb(1,51) = 50.487/100.0d0
    fa(2,51) = 3.844/10.0d0
    fb(2,51) = 13.316/100.0d0
    fa(3,51) = 2.687/10.0d0
    fb(3,51) = 2.691/100.0d0
    fa(4,51) = 0.864/10.0d0
    fb(4,51) = 0.316/100.0d0
! line 77: Sc, at. no. 21
    isour(21) = 1
    atnam(21) = "Sc"
    fa(1,21) = 3.966/10.0d0
    fb(1,21) = 88.960/100.0d0
    fa(2,21) = 2.917/10.0d0
    fb(2,21) = 20.606/100.0d0
    fa(3,21) = 1.925/10.0d0
    fb(3,21) = 3.856/100.0d0
    fa(4,21) = 0.480/10.0d0
    fb(4,21) = 0.399/100.0d0
! line 78: Se, at. no. 34
    isour(34) = 1
    atnam(34) = "Se"
    fa(1,34) = 2.298/10.0d0
    fb(1,34) = 38.830/100.0d0
    fa(2,34) = 2.854/10.0d0
    fb(2,34) = 11.536/100.0d0
    fa(3,34) = 1.456/10.0d0
    fb(3,34) = 2.146/100.0d0
    fa(4,34) = 0.590/10.0d0
    fb(4,34) = 0.316/100.0d0
! line 79: Si, at. no. 14
    isour(14) = 1
    atnam(14) = "Si"
    fa(1,14) = 2.129/10.0d0
    fb(1,14) = 57.775/100.0d0
    fa(2,14) = 2.533/10.0d0
    fb(2,14) = 16.476/100.0d0
    fa(3,14) = 0.835/10.0d0
    fb(3,14) = 2.880/100.0d0
    fa(4,14) = 0.322/10.0d0
    fb(4,14) = 0.386/100.0d0
! line 80: Sm, at. no. 62
    isour(62) = 0
    atnam(62) = "Sm"
    fa(1,62) = 5.255/10.0d0
    fb(1,62) = 28.016/100.0d0
    fa(2,62) = 4.113/10.0d0
    fb(2,62) = 5.037/100.0d0
    fa(3,62) = 1.743/10.0d0
    fb(3,62) = 0.577/100.0d0
    fa(4,62) = 0.000/10.0d0
    fb(4,62) = 0.000/100.0d0
! line 81: Sn, at. no. 50
    isour(50) = 1
    atnam(50) = "Sn"
    fa(1,50) = 3.450/10.0d0
    fb(1,50) = 59.104/100.0d0
    fa(2,50) = 3.735/10.0d0
    fb(2,50) = 14.179/100.0d0
    fa(3,50) = 2.118/10.0d0
    fb(3,50) = 2.855/100.0d0
    fa(4,50) = 0.877/10.0d0
    fb(4,50) = 0.327/100.0d0
! line 82: Sr, at. no. 38
    isour(38) = 1
    atnam(38) = "Sr"
    fa(1,38) = 5.848/10.0d0
    fb(1,38) = 104.972/100.0d0
    fa(2,38) = 4.003/10.0d0
    fb(2,38) = 19.367/100.0d0
    fa(3,38) = 2.342/10.0d0
    fb(3,38) = 3.737/100.0d0
    fa(4,38) = 0.880/10.0d0
    fb(4,38) = 0.414/100.0d0
! line 83: Ta, at. no. 73
    isour(73) = 0
    atnam(73) = "Ta"
    fa(1,73) = 5.659/10.0d0
    fb(1,73) = 28.807/100.0d0
    fa(2,73) = 4.630/10.0d0
    fb(2,73) = 5.114/100.0d0
    fa(3,73) = 2.014/10.0d0
    fb(3,73) = 0.578/100.0d0
    fa(4,73) = 0.000/10.0d0
    fb(4,73) = 0.000/100.0d0
! line 84: Tb, at. no. 65
    isour(65) = 0
    atnam(65) = "Tb"
    fa(1,65) = 5.272/10.0d0
    fb(1,65) = 29.046/100.0d0
    fa(2,65) = 4.347/10.0d0
    fb(2,65) = 5.226/100.0d0
    fa(3,65) = 1.844/10.0d0
    fb(3,65) = 0.585/100.0d0
    fa(4,65) = 0.000/10.0d0
    fb(4,65) = 0.000/100.0d0
! line 85: Tc, at. no. 43
    isour(43) = 0
    atnam(43) = "Tc"
    fa(1,43) = 4.318/10.0d0
    fb(1,43) = 28.246/100.0d0
    fa(2,43) = 3.270/10.0d0
    fb(2,43) = 5.148/100.0d0
    fa(3,43) = 1.287/10.0d0
    fb(3,43) = 0.590/100.0d0
    fa(4,43) = 0.000/10.0d0
    fb(4,43) = 0.000/100.0d0
! line 86: Te, at. no. 52
    isour(52) = 0
    atnam(52) = "Te"
    fa(1,52) = 4.785/10.0d0
    fb(1,52) = 27.999/100.0d0
    fa(2,52) = 3.688/10.0d0
    fb(2,52) = 5.083/100.0d0
    fa(3,52) = 1.500/10.0d0
    fb(3,52) = 0.581/100.0d0
    fa(4,52) = 0.000/10.0d0
    fb(4,52) = 0.000/100.0d0
! line 87: Th, at. no. 90
    isour(90) = 0
    atnam(90) = "Th"
    fa(1,90) = 6.264/10.0d0
    fb(1,90) = 28.651/100.0d0
    fa(2,90) = 5.263/10.0d0
    fb(2,90) = 5.030/100.0d0
    fa(3,90) = 2.367/10.0d0
    fb(3,90) = 0.563/100.0d0
    fa(4,90) = 0.000/10.0d0
    fb(4,90) = 0.000/100.0d0
! line 88: Ti, at. no. 22
    isour(22) = 1
    atnam(22) = "Ti"
    fa(1,22) = 3.565/10.0d0
    fb(1,22) = 81.982/100.0d0
    fa(2,22) = 2.818/10.0d0
    fb(2,22) = 19.049/100.0d0
    fa(3,22) = 1.893/10.0d0
    fb(3,22) = 3.590/100.0d0
    fa(4,22) = 0.483/10.0d0
    fb(4,22) = 0.386/100.0d0
! line 89: Tl, at. no. 81
    isour(81) = 0
    atnam(81) = "Tl"
    fa(1,81) = 5.932/10.0d0
    fb(1,81) = 29.086/100.0d0
    fa(2,81) = 4.972/10.0d0
    fb(2,81) = 5.126/100.0d0
    fa(3,81) = 2.195/10.0d0
    fb(3,81) = 0.572/100.0d0
    fa(4,81) = 0.000/10.0d0
    fb(4,81) = 0.000/100.0d0
! line 90: Tm, at. no. 69
    isour(69) = 0
    atnam(69) = "Tm"
    fa(1,69) = 5.441/10.0d0
    fb(1,69) = 29.149/100.0d0
    fa(2,69) = 4.510/10.0d0
    fb(2,69) = 5.264/100.0d0
    fa(3,69) = 1.956/10.0d0
    fb(3,69) = 0.590/100.0d0
    fa(4,69) = 0.000/10.0d0
    fb(4,69) = 0.000/100.0d0
! line 91: U, at. no. 92
    isour(92) = 1
    atnam(92) = "U"
    fa(1,92) = 6.767/10.0d0
    fb(1,92) = 85.951/100.0d0
    fa(2,92) = 6.729/10.0d0
    fb(2,92) = 15.642/100.0d0
    fa(3,92) = 4.014/10.0d0
    fb(3,92) = 2.936/100.0d0
    fa(4,92) = 1.561/10.0d0
    fb(4,92) = 0.335/100.0d0
! line 92: V, at. no. 23
    isour(23) = 1
    atnam(23) = "V"
    fa(1,23) = 3.245/10.0d0
    fb(1,23) = 76.379/100.0d0
    fa(2,23) = 2.698/10.0d0
    fb(2,23) = 17.726/100.0d0
    fa(3,23) = 1.860/10.0d0
    fb(3,23) = 3.363/100.0d0
    fa(4,23) = 0.486/10.0d0
    fb(4,23) = 0.374/100.0d0
! line 93: W, at. no. 74
    isour(74) = 0
    atnam(74) = "W"
    fa(1,74) = 5.709/10.0d0
    fb(1,74) = 28.782/100.0d0
    fa(2,74) = 4.677/10.0d0
    fb(2,74) = 5.084/100.0d0
    fa(3,74) = 2.019/10.0d0
    fb(3,74) = 0.572/100.0d0
    fa(4,74) = 0.000/10.0d0
    fb(4,74) = 0.000/100.0d0
! line 94: Xe, at. no. 54
    isour(54) = 1
    atnam(54) = "Xe"
    fa(1,54) = 3.366/10.0d0
    fb(1,54) = 35.509/100.0d0
    fa(2,54) = 4.147/10.0d0
    fb(2,54) = 11.117/100.0d0
    fa(3,54) = 2.443/10.0d0
    fb(3,54) = 2.294/100.0d0
    fa(4,54) = 0.829/10.0d0
    fb(4,54) = 0.289/100.0d0
! line 95: Y, at. no. 39
    isour(39) = 0
    atnam(39) = "Y"
    fa(1,39) = 4.129/10.0d0
    fb(1,39) = 27.548/100.0d0
    fa(2,39) = 3.012/10.0d0
    fb(2,39) = 5.088/100.0d0
    fa(3,39) = 1.179/10.0d0
    fb(3,39) = 0.591/100.0d0
    fa(4,39) = 0.000/10.0d0
    fb(4,39) = 0.000/100.0d0
! line 96: Yb, at. no. 70
    isour(70) = 0
    atnam(70) = "Yb"
    fa(1,70) = 5.529/10.0d0
    fb(1,70) = 28.927/100.0d0
    fa(2,70) = 4.533/10.0d0
    fb(2,70) = 5.144/100.0d0
    fa(3,70) = 1.945/10.0d0
    fb(3,70) = 0.578/100.0d0
    fa(4,70) = 0.000/10.0d0
    fb(4,70) = 0.000/100.0d0
! line 97: Zn, at. no. 30
    isour(30) = 1
    atnam(30) = "Zn"
    fa(1,30) = 1.942/10.0d0
    fb(1,30) = 54.162/100.0d0
    fa(2,30) = 1.950/10.0d0
    fb(2,30) = 12.518/100.0d0
    fa(3,30) = 1.619/10.0d0
    fb(3,30) = 2.416/100.0d0
    fa(4,30) = 0.543/10.0d0
    fb(4,30) = 0.330/100.0d0
! line 98: Zr, at. no. 40
    isour(40) = 0
    atnam(40) = "Zr"
    fa(1,40) = 4.105/10.0d0
    fb(1,40) = 28.492/100.0d0
    fa(2,40) = 3.144/10.0d0
    fb(2,40) = 5.277/100.0d0
    fa(3,40) = 1.229/10.0d0
    fb(3,40) = 0.601/100.0d0
    fa(4,40) = 0.000/10.0d0
    fb(4,40) = 0.000/100.0d0

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
	    allocate(hilfsfeld(newdimension),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
	    !     nur 1 mal kopieren reicht
	    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
            m1=min(newdimension,size(tf,1))
	    do i=1, m1
	      hilfsfeld(i)=tf(i)
	    enddo
	    deallocate(tf)
	    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)
	    real*8, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
            m1=min(newdimension1,size(tf,1))
            m2=min(newdimension2,size(tf,2))
	    do i1=1, m1
	      do i2=1, m2
		hilfsfeld(i1,i2) = tf(i1,i2)
	      enddo
	    enddo
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d1(tf, newdimension)
	    double complex, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
            m1=min(newdimension,size(tf,1))
	    do i=1, m1
              hilfsfeld(i) = tf(i)
            enddo
	    deallocate(tf,stat=ierr)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_c16_d2(tf, newdimension1, newdimension2)
	    complex*16, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
            m1=min(newdimension1,size(tf,1))
            m2=min(newdimension2,size(tf,2))
	    do i1=1, m1
	      do i2=1, m2
		hilfsfeld(i1,i2) = tf(i1,i2)
	      enddo
	    enddo
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d1(tf, newdimension)
	    integer*4, pointer :: hilfsfeld(:), tf(:)
	    allocate(hilfsfeld(newdimension),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
            m1=min(newdimension,size(tf,1))
	    do i=1,m1
	      hilfsfeld(i)=tf(i)
	    enddo
	    deallocate(tf)
	    tf=>hilfsfeld
	  end subroutine 

	  subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
	    integer*4, pointer :: hilfsfeld(:,:), tf(:,:)
	    allocate(hilfsfeld(newdimension1,newdimension2),stat=ierr)
	    if(ierr/=0) stop "ERROR: Not enough memory"
            m1=min(newdimension1,size(tf,1))
            m2=min(newdimension2,size(tf,2))
	    do i1=1, m1
	      do i2=1, m2
		hilfsfeld(i1,i2) = tf(i1,i2)
	      enddo
	    enddo
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
  subroutine init_struct

    use reallocate
    use defs

    implicit none

    integer                       :: ios
    double precision              :: test0

    integer                       :: index,i,j,j1,j2,m,jatom

    test0=1.D-5

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
