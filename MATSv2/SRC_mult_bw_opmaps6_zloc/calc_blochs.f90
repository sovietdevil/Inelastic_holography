subroutine calc_blochs

  use beams
  use blochs
  use defs
  use inputs

  implicit none
  
  integer i, j, k, m1
  double precision div1, div2, tmpd0d, vec(3), cart_scaprod, cart_scaprod2, det
  double complex orth, tmpd0
  
  double complex, allocatable   :: work(:), tmpmat(:,:)
  double precision, allocatable :: rwork(:)
  integer                          lwork, lrwork, liwork, info
  integer, allocatable          :: isuppz(:), iwork(:)

  m1 = -1

  allocate( gam_out(nbeams_out) )

  allocate( work(1), rwork(3*nbeams_out) )
!  call zheev('V','U',nbeams_out,mat_out,nbeams_out,gam_out,work,m1,rwork,info)
!  lwork = nint(abs(work(1))+1)
  lwork = 5000
  deallocate( work )
  allocate( work(lwork) )
  call zheev('V','U',nbeams_out,mat_out,nbeams_out,gam_out,work,lwork,rwork,info)
  deallocate( work, rwork )
  
  det = one
  do i=1, nbeams_out
    det = det*gam_out(i)
  enddo
  write(*,'(a,g16.8)') 'Diagonalization for outgoing beam completed. Determinant: ',det

!  write(*,'(f10.6,2X,$)') (gam_in(j),j=1,nbeams_in)
!  write(*,*)
!  write(*,'(f10.6,2X,$)') (abs(mat_in(zerohkl_in,j)),j=1,nbeams_in)
!  write(*,*)
!  write(*,'(f10.6,2X,$)') (gam_out(j),j=1,nbeams_out)
!  write(*,*)
!  write(*,'(f10.6,2X,$)') (abs(mat_out(zerohkl_out,j)),j=1,nbeams_out)
!  write(*,*)

!  print *, 'Matrix OUT'
!  do i = 1, nbeams_out
!    write(*,'(2f10.6,2X,$)') (mat_out(i,j), j=1, nbeams_out)
!    write(*,*)
!  enddo

! Calculate absorptive part of eigenvalues (eta's)
  if (doabsorp) then
    call calc_etas
    print *, 'Etas calculated'
  endif

! Backtransforming into Bloch coefficients...
! Note: eigenvectors are as columns, i.e. C^(j)_g(i) = mat_in(i,j)
  do i = 1, nbeams_out
    vec(:) = vchi_out(:) + aubeams_out(:,i)
    div1 = vec(3)
    div2 = vchi_out(3)
    do j = 1, nbeams_out
      mat_out(i,j) = mat_out(i,j) / sqrt(div1/div2)
    enddo
  enddo

! Note: if we abandon the ZOLZ, then the Bloch coefficients are not 'completely' orthogonal
!       sum(g) C^(i)_g* x C^(j)_g /= delta(i,j)
!       but usually the deviation from orthogonality is small (for not too high order of LZ's)

  print *, 'Backtransform to Blochs'

! Find maximum D0D
  maxd0d = dble(-1)
  do j=1, nbeams_out
    tmpd0 = mat_out(zerohkl_out,j)
    do i=1, nbeams_out
      tmpd0d = abs(mat_out(i,j)*tmpd0)
      if(tmpd0d>maxd0d) maxd0d = tmpd0d
    enddo
  enddo
  write(*,'(a,e14.6)') 'Maximum |d0d|: ', maxd0d

!  print *, 'Matrix OUT'
!  do i = 1, nbeams_out
!    write(*,'(2f10.6,2X,$)') (mat_out(i,j), j=1, nbeams_out)
!    write(*,*)
!  enddo

!  print *, '(Non-)orthogonality of Bloch-wave coefficients'
  
!  print *, 'Matrix OUT'
!  do i = 1, nbeams_out
!    do j = 1, nbeams_out
!      orth = zeroc
!      do k = 1, nbeams_out
!        orth = orth + conjg(mat_out(k,i))*mat_out(k,j)
!      enddo
!      write(*,'(2f10.6,2X,$)') orth
!    enddo
!    write(*,*)
!  enddo
  
end
