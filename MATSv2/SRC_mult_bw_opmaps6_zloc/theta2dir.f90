subroutine theta2dir(thx,thy,vz,vx,v)

  use defs
  use inputs

  implicit none

  double precision thx, thy, vz(3), vx(3), v(3)
  double precision vzn(3), vxn(3), vyn(3), p, q
  double precision cart_scaprod2

  p = cart_scaprod2(vz,vz)
  vzn(:) = vz(:)/sqrt(p)
  p = cart_scaprod2(vx,vzn)
  vxn(:) = vx(:)-p*vzn(:)
  p = cart_scaprod2(vxn,vxn)
  vxn(:) = vxn(:)/sqrt(p)
  call cart_vecprod(vzn,vxn,vyn)

  v(:) = vzn(:) + tan(thx/dble(1000))*vxn(:) + tan(thy/dble(1000))*vyn(:)
  p = cart_scaprod2(v,v)
  v(:) = v(:)/sqrt(p)

end
