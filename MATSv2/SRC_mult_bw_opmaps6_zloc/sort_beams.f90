subroutine sort_beams

  use defs
  use beams
  use inputs
  
  implicit none
  
  integer i, j
  
  external comparebeams_out, swapbeams_out
  
  print *, 'Sorting beams'

  call random_seed()

  if(nbeams_out>nb_out+1) call quicksort(nb_out+1,nbeams_out,comparebeams_out,swapbeams_out)
  
  ! incoming - just accept what we have
  nb_in = nbeams_in
  ! outgoing
  i = 1
  do while (abs(wg_out(i)).lt.w_min) 
    i = i + 1
    if(i>nbeams_out) exit
  enddo
  nb_out = i-1

  write(*,'("Number of outgoing beams with wg<wg_max: ",i4)') nb_out

  do i = 1, nb_out
    write(*,'(a,3i4,g14.6)') 'Outgoing: ', hklbeams_out(1:3,i), wg_out(i)
  enddo

  write(*,'("Number of generated beams for eigenproblem: ",i4)') nbeams_out

  do i = 1, nbeams_out
    if (hklbeams_out(1,i)==0 .and. hklbeams_out(2,i)==0 .and. hklbeams_out(3,i)==0) zerohkl_out = i
  enddo

  write(*,'("Indices of (000) beam:",i5,"/",i5)') zerohkl_in, zerohkl_out

end


integer function comparebeams_out(i,j)

  use beams
  use blochs

  implicit none

  integer i, j
  
  if (abs(wg_out(i))<abs(wg_out(j))) then
    comparebeams_out =  1
  else if (abs(wg_out(i))>abs(wg_out(j))) then
    comparebeams_out = -1
  else
    comparebeams_out =  0
  endif
  
end


subroutine swapbeams_out(i,j)

  use beams
  
  implicit none
  
  integer i, j, k(3)
  double complex tmpc
  double precision tmp, t(3)
  
  tmpc      = ug_out(i)
  ug_out(i) = ug_out(j)
  ug_out(j) = tmpc
  tmp       = sg_out(i)
  sg_out(i) = sg_out(j)
  sg_out(j) = tmp 
  tmp       = wg_out(i)
  wg_out(i) = wg_out(j)
  wg_out(j) = tmp 
  tmp        = ksi_out(i)
  ksi_out(i) = ksi_out(j)
  ksi_out(j) = tmp 
  k(1:3)              = hklbeams_out(1:3,i)
  hklbeams_out(1:3,i) = hklbeams_out(1:3,j)
  hklbeams_out(1:3,j) = k(1:3)
  t(1:3)              = aubeams_out(1:3,i)
  aubeams_out(1:3,i)  = aubeams_out(1:3,j)
  aubeams_out(1:3,j)  = t(1:3)

end
