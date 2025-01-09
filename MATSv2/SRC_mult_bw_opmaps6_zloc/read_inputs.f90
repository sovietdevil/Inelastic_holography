subroutine read_inputs

  use inputs
  use beams
  use energy
  use defs
  use struct
  use defs

  implicit none
  
  integer           i, j, k, l
  double precision  p, q, r, s, cart_length, cart_scaprod2
  double precision  as(3), aperp(3), nchi(3), axn(3)
  character*5       modein, modeout, relout
  double precision  thxin, thyin, vxin(3), cartx(3), carty(3), vzdir(3)

  read(16,*)
  read(16,*) atom, nc, lc
  read(16,*) eloss, esplit

  read(16,*)
  read(16,*)
  read(16,*) ebeam

  gamma_in  = one + ebeam/emass
  gamma_out = one + (ebeam-eloss/dble(1000))/emass
  
  chi_in  = emass/(hbarc*nminau)*sqrt(gamma_in**two-one)
  chi_out = emass/(hbarc*nminau)*sqrt(gamma_out**two-one)

  write(*,'("Gamma in vacuum for in/out: ",f12.9," /",f12.9)') gamma_in, gamma_out
  write(*,'("K(au) in vacuum for in/out: ",f12.6," /",f12.6)') chi_in, chi_out
  write(*,'("l(au) in vacuum for in/out: ",f12.6," /",f12.6)') two*pi/chi_in, two*pi/chi_out

! zone axis (za) setup; its cartesian form (cart_za) gives the z-direction
  read(16,*) (za(i),i=1,3)
  if(za(1).eq.0 .and. za(2).eq.0 .and. za(3).eq.0) then
    write(*,*) 'ERROR: Zone axis cannot be zero'
    stop
  endif
  cart_za = matmul(br1_dir,za)
  p = sqrt(dot_product(cart_za,cart_za))
  cart_za = cart_za / p

! x-direction
  read(16,*) (vxin(i),i=1,3) ! x-direction
  cartx = matmul(br1_rec,vxin)
  ! now we make sure that the x-direction is perpendicular to ZA - note that it is a recip. vector
  p = dot_product(cartx,cart_za) ! cart_za is normalized
  if(abs(p).gt.test) then
    write(*,'(a,f10.5)')  'WARNING: Zone axis for incoming beam is not perpendicular to given x-direction, scalar product', p
    cartx = cartx - p*cart_za
    write(*,'(a,3f10.5)') '         Taking only the perpendicular part of given x-direction (in cartesian coords.):',cartx(1:3)
  endif
! WARNING: this is an interesting point to think about - should we fix it or not? influences a bit the detector position set-up
!  vxin = matmul(transpose(br1_dir),cartx)/(two*pi)
! -------------------------------------------------------
  p = sqrt(dot_product(cartx,cartx))
  cartx = cartx / p
  call cart_vecprod(cart_za,cartx,carty)


! input mode for the incoming beam
  read(16,*) modein

  if(modein=='LCC  ') then
    read(16,*) (ilcc(i),i=1,3)
    ! now we make sure that lcc is perpendicular to ZA
    p = dot_product(ilcc,za)/dot_product(za,za)
    if(abs(p).gt.test) then
      write(*,'(a,f10.5)')  'WARNING: Zone axis for incoming beam is not perpendicular to Laue circle center, scalar product', p
      ilcc = ilcc - p*za
      write(*,'(a,3f10.5)') '         Taking only the perpendicular part of lcc:',ilcc(1:3)
    endif
    cart_ilcc = matmul(br1_rec,ilcc)
    thxin = atan(dot_product(cart_ilcc,cartx)/chi_in)*dble(1000)
    thyin = atan(dot_product(cart_ilcc,carty)/chi_in)*dble(1000)
  elseif(modein=='THETA') then
    read(16,*) thxin, thyin  ! in mrad
    cart_ilcc = (tan(thxin/dble(1000))*cartx + tan(thyin/dble(1000))*carty)*chi_in
    ilcc = matmul(transpose(br1_dir),cart_ilcc)/(two*pi)
  else
    write(*,'(a,a5)') 'ERROR: Unknown mode for incoming beam direction: ', modein
    stop
  endif

! set up finally the vector of incoming beam
  p = cart_length(cart_ilcc(1),cart_ilcc(2),cart_ilcc(3))
  q = sqrt(chi_in**two-p**two)
  do i=1, 3
    vchi_in(i) = -q*cart_za(i) - cart_ilcc(i)
  enddo

! diagnostic output about the incoming beam
  write(*,'(a,3i11)')   'Zone axis (hkl):          ', za(1:3)
  write(*,'(a,3f11.5)') 'Zone axis cart. direction:', cart_za(1:3)
  write(*,'(a,3f11.5)') 'X-direction (hkl):        ', vxin(1:3)
  write(*,'(a,3f11.5)') 'X-direction in cartesian: ', cartx(1:3)
  write(*,'(a,3f11.5)') 'Y-direction in cartesian: ', carty(1:3)
  write(*,'(a,3f11.5)') 'Incoming lcc (hkl):       ', ilcc(1:3)
  write(*,'(a,3f11.5)') 'Incoming lcc in cartesian:', cart_ilcc(1:3)
  write(*,'(a,3f11.5)') 'Resulting incoming beam:  ', vchi_in(1:3)
  write(*,'(a,2f11.5)') 'Tilt of incoming beam vs zone axis (mrad): ', thxin, thyin


! skip empty line and header
  read(16,*)
  read(16,*)

! input mode for the outgoing beam
  read(16,*) relout
  read(16,*) modeout
  hkl_det = vxin

  if(relout=='RELZA') then
    vzero_out = -cart_za*chi_out
    write(*,'(a,3f11.5)') 'Outgoing reference beam (zone axis): ', vzero_out(1:3)
    ! in this case, cartx and carty are already set properly
  elseif(relout=='RELTB') then
    vzero_out = vchi_in*chi_out/chi_in
    write(*,'(a,3f11.5)') 'Outgoing reference beam (transmitted beam): ', vzero_out(1:3)
    ! in this case we need to re-set the cartx and carty vectors, so that they are perpendicular to the vzero_out
    cartx = matmul(br1_rec,vxin)
    p = dot_product(cartx,cart_za)
    cartx = cartx - p*cart_za
    p = sqrt(dot_product(cartx,cartx))
    cartx = cartx / p
    call cart_vecprod(cart_za,cartx,carty)
  else
    write(*,'(a,a5)') 'ERROR: Unknown reference beam for outgoing beam direction: ', relout
    stop
  endif
  p = sqrt(dot_product(vzero_out,vzero_out))
  vzdir = vzero_out/p

! find two-fold Bragg angle corresponding to hkl_det
  as = matmul(br1_rec,vxin)
  write(*,'(a,3i3,a,3f9.5)') 'Angles are given w.r.t. ',hkl_det(1:3),', in cart. coord.: ',as(1:3)
  axn(:) = vzero_out(:)+as(:)
  s = 1000*acos((cart_scaprod2(vzero_out,vzero_out)+cart_scaprod2(axn,axn)-cart_scaprod2(as,as)) / (two*chi_out*sqrt(cart_scaprod2(axn,axn))))
  write(*,'(a,f9.5)') 'Corresponding two-fold Bragg angle (mrad): ', s

! get detector position in mrads w.r.t. reference direction vzero_out
  if(modeout=='THETA') then
    read(16,*) thetax, thetay
  elseif(modeout=='MULTG') then
    read(16,*) thetax, thetay
    thetax = s*thetax
    thetay = s*thetay
  else
    write(*,'(a,a5)') 'ERROR: Unknown mode for outgoing beam direction: ', modeout
    stop
  endif

! set up the vector of outgoing beam
  call theta2dir(thetax,thetay,vzero_out,as,vchi_out)
  vchi_out(:) = chi_out*vchi_out(:)

! find corresponding output lcc
  p = cart_scaprod2(vchi_out,cart_za)
  cart_olcc = vchi_out - p*cart_za
  olcc = matmul(transpose(br1_dir),cart_olcc)/(two*pi)

! diagnostic output about the outgoing beam
  write(*,*)
  write(*,'(a,2f11.5)') 'Detector position (th_x, th_y) in mrad:', thetax, thetay
  write(*,'(a,2f11.5)') 'Detector position in multiples of G:   ', thetax/s, thetay/s
  write(*,'(a,3f11.5)') 'X-direction (hkl):        ', vxin(1:3)
  write(*,'(a,3f11.5)') 'X-direction in cartesian: ', cartx(1:3)
  write(*,'(a,3f11.5)') 'Y-direction in cartesian: ', carty(1:3)
  write(*,'(a,3f11.5)') 'Outgoing lcc (hkl):       ', olcc(1:3)
  write(*,'(a,3f11.5)') 'Outgoing lcc in cartesian:', cart_olcc(1:3)
  write(*,'(a,3f11.5)') 'Resulting outgoing beam:  ', vchi_out(1:3)

! sanity check
  p = olcc(1)*za(1)+olcc(2)*za(2)+olcc(3)*za(3)
  p = p/dble(maxval(za(:)))
  if(abs(p).gt.test) then
    write(*,*) 'ERROR: Zone axis for outgoing beam is not perpendicular to Laue circle center', p
    stop
  endif

! for outgoing beam the "natural" selection of lcc with respect to crystal coordinates
! needs to be inverted (propagation in -k is equivalent to propagation in k in inverted crystal)
  olcc(:)      = -olcc(:)
  cart_olcc(:) = -cart_olcc(:)
! is it still necessary?

! Later note:
! this reversion of sign is connected with construction of secular matrices - one puts in nondiagonal
! elements for outgoing beam matrix V(-g)=V*(g) instead of V(g)
! one can remove these sign inversions and obtain more 'natural' beams, but them matrix construction
! has to be changed also => performed!

  read(16,*)
  read(16,*)
  read(16,*) vhkltype
  if (vhkltype=='DOYLE') then
    write(*,'(a)') 'Doyle-Turner (Acta Cryst. A24 (1968), 390) scattering factors are used for V(hkl)'
  elseif (vhkltype=='WEICK') then
    write(*,'(a)') 'Weickenmeier-Kohl (Acta Cryst. A47 (1991), 590) scattering factors are used for V(hkl)'
  elseif (vhkltype=='WIEN ') then
    write(*,'(a)') 'First principles V(hkl) from WIEN2k are used'
  else
    write(*,'(a,a5)') 'Unrecognized V(hkl) option: ', vhkltype
    stop
  endif
  read(16,*) abstype
  read(16,*) abscoeff
  if(abstype=='COEFF') then
    write(*,'(a,f10.4)') 'Absorptive potential will be proportional to V(hkl) with parameter ', abscoeff
  elseif (abstype=='WEICK') then
    write(*,'(a)') 'Weickenmeier-Kohl (Acta Cryst. A47 (1991), 590) absorptive potential'
  else
    write(*,'(a,a5)') 'Unrecognized absorption type: ', abstype
  endif
  if (abstype=='WEICK' .or. abscoeff>test) then
    write(*,'(a)') 'Absorption calculation will be performed'
    doabsorp = .true.
  else
    write(*,'(a)') 'Absorption calculation will NOT be performed'
    doabsorp = .false.
  endif

  read(16,*)
  read(16,*)
  read(16,*) bm_eig
  if(bm_eig/='AUTO') stop
  read(16,*) hmax, kmax, lmax
  read(16,*) g_max
  read(16,*) w_max

  read(16,*)
  read(16,*)
  read(16,*) bm_sum
  if(bm_sum/='AUTO') stop
  read(16,*) w_min
  read(16,*) c0_min

  if (g_max.gt.0 .and. cart_length(cart_ilcc(1),cart_ilcc(2),cart_ilcc(3)).gt.g_max) then
    write(*,'("ERROR: G_max too small (g_max must be > |lcc_in|)")')
    stop
  endif

  if (g_max.gt.0 .and. cart_length(cart_olcc(1),cart_olcc(2),cart_olcc(3)).gt.g_max) then
    write(*,'("ERROR: G_max too small (g_max must be > |lcc_out|)")')
    stop
  endif

  read(16,*)
  read(16,*)
  read(16,*) efermi
  read(16,*) emin_mesh, estep, ne_mesh
  read(16,*) emin, emax

  iemin = 1
  p = emin_mesh
  do while (p<emin)
    iemin = iemin+1
    p = emin_mesh + (iemin-1)*estep
  enddo

  iemax = ne_mesh
  p = emin_mesh + (ne_mesh-1)*estep
  do while (p>emax)
    iemax = iemax-1
    p = emin_mesh + (iemax-1)*estep
  enddo

  ne = iemax-iemin+1
  
  write(*,'(a,f8.4,a,f8.4)') 'Will treat energy indices between Emin=',emin,' and Emax=',emax
  write(*,'(a,i4,a,i4)') 'This corresponds to indices iemin=',iemin,' and iemax=',iemax

  read(16,*)
  read(16,*) tmin, tstep, nt
  
end
