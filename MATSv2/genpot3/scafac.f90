double precision function scafac(sx2,m)

  use asff
  use struct
  use inputs
  use defs
        
  implicit none
        
  integer          imax, i, m
  double precision sx2, ar
        
  imax = 3
  if (isour(nint(zz(m)))==1) imax=4
  scafac = 0.0d0
  do i = 1, imax
    ar = -fb(i,nint(zz(m))) * sx2
    if (dabs(ar).le.10.) scafac = scafac + fa(i,nint(zz(m))) * dexp(ar)
  enddo

! convert scafac (Doyle, Turner, Acta Cryst. A24 (1968), 390) from nm to Volts!
  scafac = fuccon * scafac / (vol*auinnm**three)

  return

end
