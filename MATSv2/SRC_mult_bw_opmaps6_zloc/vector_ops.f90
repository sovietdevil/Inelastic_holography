subroutine hkl2cart(h,k,l,x,y,z)

  use defs
  use struct
  
  implicit none
  
  integer h, k, l
  double precision x,y,z

  x = dble(h)*br1_rec(1,1) + dble(k)*br1_rec(1,2) + dble(l)*br1_rec(1,3) 
  y = dble(h)*br1_rec(2,1) + dble(k)*br1_rec(2,2) + dble(l)*br1_rec(2,3) 
  z = dble(h)*br1_rec(3,1) + dble(k)*br1_rec(3,2) + dble(l)*br1_rec(3,3) 

end


subroutine cart2hkl(x,y,z,h,k,l)

  use defs
  use struct
  
  implicit none
  
  double precision x,y,z,h,k,l

  h = x*br1_dir(1,1) + y*br1_dir(2,1) + z*br1_dir(3,1) 
  k = x*br1_dir(1,2) + y*br1_dir(2,2) + z*br1_dir(3,2) 
  l = x*br1_dir(1,3) + y*br1_dir(2,3) + z*br1_dir(3,3) 

  h = h/(two*pi)
  k = k/(two*pi)
  l = l/(two*pi)

end


double precision function cart_length(x,y,z)

  use defs

  implicit none
  
  double precision x, y, z
  
  cart_length = sqrt(x**two+y**two+z**two)

end


double precision function cart_length2(v)

  use defs

  implicit none
  
  double precision v(3)
  
  cart_length2 = sqrt(v(1)**two+v(2)**two+v(3)**two)

end


double precision function hkl_length(h,k,l)
        
  use defs
        
  implicit none
        
  double precision    x, y, z, cart_length
  integer             h, k, l

  call hkl2cart(h,k,l,x,y,z)
  hkl_length = cart_length(x,y,z)
  
end


double precision function cart_scaprod(x1,y1,z1,x2,y2,z2)

  implicit none
  
  double precision x1, y1, z1, x2, y2, z2
  
  cart_scaprod = x1*x2 + y1*y2 + z1*z2

end


double precision function cart_scaprod2(v1,v2)

  implicit none
  
  double precision v1(3), v2(3), sp
  integer i
  
  sp = 0.0d0
  do i=1, 3
    sp = sp + v1(i)*v2(i)
  enddo
  
  cart_scaprod2 = sp

end


subroutine cart_vecprod(v1,v2,v3)

  implicit none

  double precision v1(3), v2(3), v3(3), v(3)

  v(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v(3) = v1(1)*v2(2) - v1(2)*v2(1)

  v3(:) = v(:)

end
