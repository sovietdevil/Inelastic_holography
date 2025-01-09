logical function newbeam(h,k,l,io)

  use beams

  implicit none

  integer h, k, l, io

  integer i

  newbeam = .true.

  if (io==1) then
    do i=1, nb_in
      if (hklbeams_in(1,i)==h .and. &
          hklbeams_in(2,i)==k .and. &
          hklbeams_in(3,i)==l) newbeam=.false.
    enddo
  else
    do i=1, nb_out
      if (hklbeams_out(1,i)==h .and. &
          hklbeams_out(2,i)==k .and. &
          hklbeams_out(3,i)==l) newbeam=.false.
    enddo
  endif

end
