subroutine generate_beams

  use struct
  use defs
  use inputs
  use beams
  use reallocate
  
  implicit none
  
  integer i, j, k, l, m, n, scp
  double precision len, x, y, z

  double precision cart_length
  logical newbeam

  zerohkl_out = 0
  nb_in = 0
  nb_out = 0

  n = (2*hmax+1)*(2*kmax+1)*(2*lmax+1)

  if (bm_sum=='FILE') then
  ! there are already loaded some beams
    call doreallocate(hklbeams_out,3,n)
    call doreallocate(aubeams_out,3,n)
  else
    allocate( hklbeams_out(3,n),  aubeams_out(3,n) )
  endif
  
!  nbeams_in  = nb_in
  nbeams_out = nb_out
  do i = -hmax, hmax 
    do j = -kmax, kmax
      do k = -lmax, lmax
        scp = i*za(1) + j*za(2) + k*za(3)
        call hkl2cart(i,j,k,x,y,z)
        if (cart_length(x-cart_olcc(1),y-cart_olcc(2),z-cart_olcc(3)).lt.(test+g_max) .and. newbeam(i,j,k,2) .and. (holz .or. scp==0)) then
	  nbeams_out = nbeams_out + 1
	  hklbeams_out(1,nbeams_out) = i
	  hklbeams_out(2,nbeams_out) = j
	  hklbeams_out(3,nbeams_out) = k
	  aubeams_out(1,nbeams_out) = x
	  aubeams_out(2,nbeams_out) = y
	  aubeams_out(3,nbeams_out) = z
	  if ((i==0).and.(j==0).and.(k==0)) zerohkl_out = nbeams_out
!          write(*,'("Outgoing: ",i5,":",3i4,3f10.5)') nbeams_out, i, j, k, x, y, z
	endif
      enddo
    enddo
  enddo
  write(*,'("Total generated beams (in/out): ",i5,"/",i8)') nbeams_in, nbeams_out

end
