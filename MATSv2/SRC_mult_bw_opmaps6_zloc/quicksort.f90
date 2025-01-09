recursive subroutine quicksort(first,last,comp,swap)

  implicit none
  
  integer first, last
  
  integer ipivot, istore, ind
  real rnd

  interface

    integer function comp(i,j)
      integer, intent(IN) :: i,j
    end function

    subroutine swap(i,j)
      integer, intent(IN) :: i,j
    end subroutine

  end interface

  call random_number(rnd)
  ipivot = nint(first+rnd*(last-first))

  if(ipivot/=last) call swap(ipivot,last)

  istore = first
  do ind=first, last-1
    if(comp(ind,last)>0) then
      call swap(ind,istore)
      istore = istore + 1
    endif
  enddo

  call swap(istore,last)

!  print *, first, istore, last

  if(first<istore-1) call quicksort(first,istore-1,comp,swap)
  if(istore+1<last)  call quicksort(istore+1,last,comp,swap)

end
