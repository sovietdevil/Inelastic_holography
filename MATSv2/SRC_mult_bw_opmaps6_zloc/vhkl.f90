double complex function vhkl(h,k,l)

  use defs
  use struct
  use inputs
  use beams
! TODO: Debye-Weller
!        use inputs
      
  implicit none
      
  double precision   trans(4,3),vr,vi,sma,smb
  integer            h,k,l
  integer            icent, ic
  integer            m, i, j, nl, ind
  double precision   xat, yat, zat, xato, yato, zato, xn, yn, zn
  double precision   sx, sx2, fh, fk, fl
  double precision   hkl_length, scafac, arg, smfac, fscat
  character*2        zsymb
  complex            fscatt
  double complex     vhklwien

! if we use WIEN2k potential, things are simpler...
  if (vhkltype=='WIEN ') then
    vhkl = vhklwien(h,k,l)
    return
  endif
! ... the rest is for potentials, which need summation over atoms in unit cell


! prepare coefficient, so that result of uhkl is in 1/a.u.^2
!  coefv= auinnm**two/(pi*fuccon)
  sx   = hkl_length(h,k,l)/two
! scafac works internally in nm without 2pi
  sx2  = (sx*nminau/(two*pi))**two
  fh   = dble(h)
  fk   = dble(k)
  fl   = dble(l)
  trans(:,:) = zero

  select case (lattic(1:3))
    case ('P  ') 
      icent = 1
    case ('H  ') 
      icent = 1
    case ('B  ') 
      icent = 2
      trans(2,1)=half
      trans(2,2)=half
      trans(2,3)=half
    case ('R  ') 
      icent = 3
      trans(2,1)=one/three
      trans(2,2)=two/three
      trans(2,3)=trans(2,2)
      trans(3,1)=trans(2,2)
      trans(3,2)=trans(2,1)
      trans(3,3)=trans(2,1)
    case ('F  ') 
      icent = 4
      trans(2,1)=half
      trans(2,2)=half
      trans(2,3)=zero
      trans(3,1)=half
      trans(3,2)=zero
      trans(3,3)=half
      trans(4,1)=zero
      trans(4,2)=half
      trans(4,3)=half
    case ('CYZ') 
      icent = 2
      trans(2,1)=zero
      trans(2,2)=half
      trans(2,3)=half
    case ('CXZ') 
      icent = 2
      trans(2,1)=half
      trans(2,2)=zero
      trans(2,3)=half
    case ('CXY') 
      icent = 2
      trans(2,1)=half
      trans(2,2)=half
      trans(2,3)=zero
  end select
        
  vhkl = zeroc
!  write(*,*) icent, nat, mult(1), ndif
  vr = zero
  vi = zero
  ind = 1
  do m = 1, nat
    sma = zero
    smb = zero
    nl  = mult(m)
    do i = 1, nl
      xat = pos(1,ind)
      yat = pos(2,ind)
      zat = pos(3,ind)
      ind = ind + 1
      do ic=1,icent
        xato = xat + trans(ic,1)
        yato = yat + trans(ic,2)
        zato = zat + trans(ic,3)
        xato = mod(xato+dble(5),one)
        xato = mod(xato+dble(5),one)
        xato = mod(xato+dble(5),one)
        arg   =-two*pi*(fh*xato+fk*yato+fl*zato)
!        write(*,'(9f12.6)') fh, fk, fl, arg, two, pi, xato, yato, zato
! TODO: Debye-Weller
!        smfac = dfloat(mult(m))/dfloat(ndif)*dexp(-btemp(ind)*sx2)
        smfac = one
        sma   = sma + cos(arg) * smfac
        smb   = smb + sin(arg) * smfac
!        write(*,'(3i3,7e18.8)') m, i, ic, sma, smb, arg, xato, yato, zato
      enddo
    enddo

!   now we evaluate the scattering factor using some of implemented approaches
    if (vhkltype=='DOYLE') then
      fscat = scafac(sx2, m)
    elseif (vhkltype=='WEICK') then
      fscat = fscatt(real(sx*two*nminau/10.0),0.0,nint(zz(m)),zsymb,200.0,.false.,.false.,.false.) &
                      / (dble(40)*pi*vol*auinnm**three/fuccon)
    else
      write(*,'(a,a5)') 'Non-recognized potential type: ', vhkltype
      stop
    endif

    vr = vr + fscat * sma
    vi = vi + fscat * smb
!    write(*,'(5g14.6)') fscat, sx2, vr, vi, vol
  enddo
  if (dabs(vr).le.1.d-10) vr=zero
  if (dabs(vi).le.1.d-10) vi=zero
!  write(*,'(3i3,2g14.6)') h, k, l, gamma_in*vr, vi
  vhkl=dcmplx(vr,vi)
  return
end
