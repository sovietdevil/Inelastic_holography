subroutine calc_ws

  use defs
  use inputs
  use beams
  
  implicit none
  
  integer i, j, k
  double precision v(3), ilcc2, olcc2, cart_scaprod2
  
  allocate( sg_out(nbeams_out), ksi_out(nbeams_out), wg_out(nbeams_out) )

  olcc2 = cart_scaprod2(cart_olcc,cart_olcc)
  vk_out(1) = - sqrt(k_out**two-olcc2)*cart_za(1)-cart_olcc(1)
  vk_out(2) = - sqrt(k_out**two-olcc2)*cart_za(2)-cart_olcc(2)
  vk_out(3) = - sqrt(k_out**two-olcc2)*cart_za(3)-cart_olcc(3)
  write(*,'("Outgoing |k| and k: ",4f10.5)') k_out, vk_out(1:3)

  do i = 1, nbeams_out
  
    v(:) = vk_out(:) + aubeams_out(:,i)
    sg_out(i) = k_out-sqrt(cart_scaprod2(v,v))
    if (abs(ug_out(i)).lt.test) then
      ksi_out(i) = one/test
      wg_out(i)  = one/test
    else
      ksi_out(i) = k_out*two*pi/abs(ug_out(i))
      wg_out(i)  = sg_out(i) * ksi_out(i)
    endif
!    write(*,'("Outgoing: ",3i3,3g14.6)') hklbeams_out(1,i),hklbeams_out(2,i),hklbeams_out(3,i), sg_out(i), ksi_out(i), wg_out(i)
  
  enddo

end
