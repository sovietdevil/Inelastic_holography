module hashtable
! rewritten for fortran, following Blog post by Jeff Preshing:
! http://preshing.com/20130107/this-hash-table-is-faster-than-a-judy-array/

  type cell
    integer key
    integer value
  end type cell
  
  type(cell), pointer :: m_cells(:)
  integer m_arraySize
  integer m_population
  logical m_zeroUsed
  type(cell), target :: m_zeroCell

 contains
  
  ! init hash table, initial size must be power of 2
  subroutine hash_init(init_size)
  
    integer init_size, i
    
    m_arraySize = init_size
    if(iand(m_arraySize,m_arraySize-1)/=0) then
      print *, 'ERROR: initial size of array must be power of two, but was given', m_arraySize
      stop
    endif
    allocate(m_cells(0:m_arraySize-1))
    do i=0, init_size-1
      m_cells(i)%key = 0
      m_cells(i)%value = 0
    enddo
    m_population = 0
    
    m_zeroUsed = .false.
    m_zeroCell%key = 0
    m_zeroCell%value = 0
    
  end subroutine hash_init
  
  ! destructor of hash table
  subroutine hash_destroy()
    deallocate(m_cells)
  end subroutine hash_destroy
  
  ! lookup for a key and return the cell
  function hash_lookup(key)
    type(cell), pointer :: hash_lookup

    integer key, icell
    type(cell), pointer :: curr_cell
    
    if(key/=0) then
      icell = iand(integerHash(key),m_arraySize-1)
      do while(.true.)
        curr_cell => m_cells(icell)
!if(key==786432*3) print *, key, icell, curr_cell%key, curr_cell%value
        if(curr_cell%key==0) then
!print *, 'Failed', key, icell
          hash_lookup => null()
          exit
        else
          if(curr_cell%key==key) then
            hash_lookup => curr_cell
            exit
          else
            icell = icell + 1
            if(icell==m_arraySize) icell = 0
          endif
        endif
      enddo
    else
      if(m_zeroUsed) then
        hash_lookup => m_zeroCell
      else
        hash_lookup => null()
      endif
    endif
  
  end function hash_lookup

  ! insert key, return the cell
  function hash_insert(key)
    type(cell), pointer :: hash_insert
    
    integer key, icell
    type(cell), pointer :: curr_cell => null()
    
    if(key/=0) then
      icell = iand(integerHash(key),m_arraySize-1)
      do while(.true.)
        curr_cell => m_cells(icell)
!if(key==786432*3) print *, key, icell, curr_cell%key, curr_cell%value
        if(curr_cell%key==key) then
          hash_insert => curr_cell
          exit
        else
          if(curr_cell%key==0) then
            m_population = m_population + 1
            curr_cell%key = key
            if(m_population/3*4>=m_arraySize) then
              call hash_repopulate(m_arraySize*2)
              hash_insert => hash_lookup(key)
              exit
            else
              hash_insert => curr_cell
              exit
            endif
          endif
        endif
        icell = icell + 1
        if(icell==m_arraySize) icell = 0
      enddo
    else
      if(.not.m_zeroUsed) then
        m_zeroUsed = .true.
        m_population = m_population + 1
        if(m_population/3*4>=m_arraySize) call hash_repopulate(m_arraySize*2)
      endif
      hash_insert => m_zeroCell
    endif
  
  end function hash_insert
  
  !
  subroutine hash_repopulate(new_size)
  
    integer new_size, icell, jcell, key
    type(cell), save, pointer :: new_cells(:)
    type(cell), pointer :: curr_cell
    
    if(new_size<4) new_size=4
    
    if(iand(new_size,new_size-1)/=0) then
      print *, 'ERROR: in repopulation, the desired new size is not power of two', new_size
      stop
    endif
    if(m_population/3*4>new_size) then
      print *, 'ERROR: in repopulation, the desired new size is too small', new_size, m_population
      stop
    endif
    
!print *, 'Allocating new cells, size:', new_size
    allocate(new_cells(0:new_size-1))
!print *, 'Allocated new cells, zeroing them...'
    do icell=0, new_size-1
      new_cells(icell)%key = 0
      new_cells(icell)%value = 0
    enddo
!print *, 'New cells ready', new_size
    
    do icell=0, m_arraySize-1
      key = m_cells(icell)%key
      if(key/=0) then
        jcell = iand(integerHash(key),new_size-1)
        do while(.true.)
          curr_cell => new_cells(jcell)
          if(curr_cell%key==0) then
            curr_cell = m_cells(icell)
            exit
          else
            jcell = jcell + 1
            if(jcell==new_size) jcell=0
          endif
        enddo
      endif
    enddo
    
    m_arraySize = new_size
    deallocate(m_cells)
    m_cells => new_cells
  
  end subroutine hash_repopulate

  ! clearing the hash
  subroutine hash_clear
    
    integer icell
    
    do icell=0, m_arraySize-1
      m_cells(icell)%key = 0
      m_cells(icell)%value = 0
    enddo
    m_population = 0
    m_zeroUsed = .false.
    m_zeroCell%value = 0
    
  end subroutine hash_clear
  
  ! TODO: compacting the hash
  subroutine hash_compact
  
    integer k, l
    
    if(m_population<1024) then
      k = (4*m_population+3)/3
    else
      k = m_population/3*4+3
    endif
    k = k - 1
    l = ishft(k,-1)
    k = ior(k,l)
    l = ishft(k,-2)
    k = ior(k,l)
    l = ishft(k,-4)
    k = ior(k,l)
    l = ishft(k,-8)
    k = ior(k,l)
    l = ishft(k,-16)
    k = ior(k,l)
    k = k + 1
    call hash_repopulate(k)
  
  end subroutine hash_compact
  
  ! erase a key/value pair from the hash
  subroutine hash_delete(key)
  
    integer key, icell, jcell, ideal, off_idc, off_idn
    type(cell), pointer :: curr_cell, neighbor
    
    if(key/=0) then
      icell = iand(integerHash(key),m_arraySize-1)
      ! first we find the cell (icell)
      do while(.true.)
        curr_cell => m_cells(icell)
        if(curr_cell%key==0) then
          print *, 'ERROR: attempted to delete key, which was not found', key
          stop
        else
          if(curr_cell%key==key) then
            exit ! we found the element to delete
          else
            icell = icell + 1
            if(icell==m_arraySize) icell = 0
          endif
        endif
      enddo
      ! now we erase it, reshuffling the nearby elements, if needed
      jcell = icell
      do while(.true.)
        jcell = jcell + 1
        if(jcell==m_arraySize) jcell = 0
        neighbor => m_cells(jcell)
        if(neighbor%key==0) then ! no need of (more) shuffling, just delete
          curr_cell%key = 0
          curr_cell%value = 0
          m_population = m_population - 1
          exit
        endif
        ideal = iand(integerHash(neighbor%key),m_arraySize-1)
!write(*,'(8i10)') key, icell, curr_cell%key, curr_cell%value, jcell, neighbor%key, neighbor%value, ideal
        off_idc = icell - ideal
        if(icell<ideal) off_idc = off_idc + m_arraySize
        off_idn = jcell - ideal
        if(jcell<ideal) off_idn = off_idn + m_arraySize
        if(off_idc<off_idn) then
          curr_cell = neighbor
          curr_cell => m_cells(jcell)
          icell = jcell
        endif
      enddo
    else
      if(m_zeroUsed) then
        m_zeroUsed = .false.
        curr_cell%value = 0
        m_population = m_population - 1
      else
        print *, 'ERROR: attempted to delete the unused zero cell'
        stop
      endif
    endif
    
  end subroutine hash_delete
  
  ! based on code.google.com/p/smhasher/wiki/MurmurHash3
  integer function integerHash(h)
    
    integer h, k, l
    
    k = ishft(h,-16)
    l = ieor(h,k)
!    l = l*Z'85ebca6b' ! hexadecimal constant
    l = l*B'10000101111010111100101001101011' ! binary constant
    k = ishft(l,-13)
    l = ieor(l,k)
!    l = l*Z'c2b2ae35' ! hexadecimal constant
    l = l*B'11000010101100101010111000110101' ! binary constant
    k = ishft(l,-16)
    l = ieor(l,k)

    integerHash = l
    
  end function integerHash

  ! NOTE: the Iterator part has been skipped for the moment (TODO?)

end module hashtable
