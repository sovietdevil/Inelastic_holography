subroutine select_beams

  use defs
  use beams

  implicit none

  integer i, j, k
  
  write(*,*)
  write(*,'("Outgoing beams, that survived:")')
  write(*,'(a)')  " h   k   l      Re[Ug]        Im[Ug]          Sg          ksi_G          w_g"
  write(*,'(a)')  "--------------------------------------------------------------------------------"

  k = 0
  do i = 1, nbeams_out
    if (abs(wg_out(i)).lt.w_max) then
      k = k + 1
      hklbeams_out(:,k) = hklbeams_out(:,i)
      aubeams_out(:,k)  = aubeams_out(:,i)
      ug_out(k)         = ug_out(i)
      sg_out(k)         = sg_out(i)
      ksi_out(k)        = ksi_out(i)
      wg_out(k)         = wg_out(i)
      if(i.eq.zerohkl_out) then
        zerohkl_out = k
      endif
      write(*,'(3i4,5g14.6)') hklbeams_out(:,k), ug_out(k), sg_out(k), ksi_out(k), wg_out(k)
    endif
  enddo
  nbeams_out = k

  write(*,'("Remaining (hkl) beams: ",i4)') nbeams_out
  write(*,*)

  write(19,*)
  write(19,'(i5)') nbeams_out
  do i = 1, nbeams_out
    write(19,'(3i4)') hklbeams_out(:,i)
  enddo

end
