double complex function vabs(h,k,l)

  use defs
  use struct
  use beams

  implicit none

  integer h, k, l

  double complex vhkl, vhklpot

  if (abstype=='COEFF') then
    vhklpot = vhkl(h,k,l)
    vabs = dcmplx(abscoeff)*vhklpot
  elseif (abstype=='WEICK') then
    write(*,'(a)') 'Not implemented yet...'
    stop
  else
    vabs = zeroc
  endif

end
