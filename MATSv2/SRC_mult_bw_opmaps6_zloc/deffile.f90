subroutine deffile

  implicit none

  character*80 deffn, fname
  character*11 status, form
  integer ios, iunit

  if(iargc()==1) then
    call getarg(1,deffn)
  else
    write(*,'(a)') 'WARNING: Awaiting exactly one command line parameter.'
    write(*,'(a)') '         Assuming default definition filename - dyndif.def'
    deffn='dyndif.def'
  endif
  open(1,file=deffn,status='old',iostat=ios)
  if (ios/=0) then
    write(*,'(a,a)') 'ERROR: Opening definition file failed'
    stop
  endif
  read(1,*,iostat=ios) iunit, fname, status, form
  if(ios/=0) then
    write(*,'(a,a)') 'ERROR: Reading definition file failed, ',fname
    stop
  endif
  do while(ios==0)
    open(iunit,file=fname,status=status,form=form,iostat=ios)
    if(ios/=0) then
      write(*,'(a,a,a)') 'ERROR: Opening file failed, ',fname
      stop
    endif
    read(1,*,iostat=ios) iunit, fname, status, form
    if(ios>0) then
      write(*,'(a)') 'ERROR: Reading definition file failed'
      stop
    endif
  enddo

end
