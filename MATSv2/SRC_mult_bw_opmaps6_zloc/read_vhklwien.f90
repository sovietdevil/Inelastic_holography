subroutine read_vhklwien

  use defs
  use inputs
  use struct
  use beams

  implicit none

  integer ios, h, k, l, maxh, maxk, maxl, i, nz
  double precision dis, vre, vim, coef

  coef = four*pi*auinnm**two/fuccon

  w2k_num = 0
  nz = 0
  read(25,*,iostat=ios) h, k, l, dis, vre, vim
  do while (ios==0)
    w2k_num = w2k_num + 1
    read(25,*,iostat=ios) h, k, l, dis, vre, vim
    if (h.ne.0 .or. k.ne.0 .or. l.ne.0) nz = nz + 1
  enddo

  if (w2k_num==0) then
    write(*,'(a)') 'ERROR: Wien2k potentials are not present.'
    write(*,'(a)') '       Please run "extract_structure_factors.pl" to calculate them.'
    write(*,'(a,f8.2)') '       G(max) corresponds to s(max)=', g_max/(four*pi*auinnm*dble(10))
    stop
  endif

  allocate( w2k_vg(w2k_num+nz,3), w2k_hkl(w2k_num+nz,3) )
  rewind(25)
  k = 0
  do i=1, w2k_num
    k = k + 1
    read(25,*) w2k_hkl(k,1:3), w2k_vg(k,1:3)
    if (w2k_hkl(k,1).ne.0 .or. w2k_hkl(k,2).ne.0 .or. w2k_hkl(k,3).ne.0) then
      k = k + 1
      w2k_hkl(k,1:3) = -w2k_hkl(k-1,1:3)
      w2k_vg(k,1:3)  =  w2k_vg(k-1,1:3)
    endif
  enddo
  ! we want to have atomic units for lengths: conversion s(Ang)->G(a.u.)
  w2k_vg(:,1) = w2k_vg(:,1)*four*pi*auinnm*dble(10)
  ! we want Volts for potential
  w2k_vg(:,2) = -w2k_vg(:,2)/vol/coef
  w2k_vg(:,3) = -w2k_vg(:,3)/vol/coef
!  print *, vol, coef

  write(*,'(a,i6,a,f8.2)') &
    'Loaded ', w2k_num, ' WIEN2k V(hkl) components with max(G) = ', dis*four*pi*auinnm*dble(10)

  if(dis*four*pi*auinnm*dble(10)<g_max) then
    write(*,'(a,f8.2,a)') 'WARNING: That may be not enough, because G(max)=', g_max, ' is larger'
    write(*,'(a,f8.2)') '         G(max) corresponds to s(max)=', g_max/(four*pi*auinnm*dble(10))
  endif

end
