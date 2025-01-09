double complex function vhklwien(h,k,l)

  use defs
  use struct
  use beams

  implicit none

  integer h, k, l

  double precision hkl_length, lng
  integer i, j, kx, ky, kz, px, py, pz
  logical found

  lng = hkl_length(h,k,l)

  found = .false.
  do i=1, w2k_num
    if (abs(lng-w2k_vg(i,1))<1.0d-5) then
      kx = w2k_hkl(i,1)
      ky = w2k_hkl(i,2)
      kz = w2k_hkl(i,3)
      do j=1, nsym
        px = iz(1,1,j)*kx + iz(1,2,j)*ky + iz(1,3,j)*kz
        py = iz(2,1,j)*kx + iz(2,2,j)*ky + iz(2,3,j)*kz
        pz = iz(3,1,j)*kx + iz(3,2,j)*ky + iz(3,3,j)*kz
!        write(*,'(9i4)') h, k, l, kx, ky, kz, px, py, pz
        if (px==h .and. py==k .and. pz==l) then
          found = .true.
          exit
        endif
      enddo
    endif
    if (found) then
!      write(*,'(3i4,a,3i4)') h, k, l, ' => ', w2k_hkl(i,1:3)
      vhklwien = dcmplx(w2k_vg(i,2),w2k_vg(i,3))
      return
    endif
  enddo

!  print *, 'Not found:', h, k, l
  if(.not.found) vhklwien = zeroc

!  if(.not.found) then
!    write(*,'(a,3i4)') 'ERROR: Potential component not found for beam ', h, k, l
!    stop
!  endif

end
