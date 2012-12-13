module util_search
  use util_runtime,only:err_exit
  use util_memory,only:reallocate
  use util_string,only:integer_to_string
  implicit none
  private
  public::LookupTable,initiateTable,addToTable,indexOf

  type LookupTable
     integer::size !number of element in the table
     integer,allocatable::list(:) !correspondence table
  end type LookupTable

CONTAINS
  !> \brief initiate table with specified initial size
  !>
  !> The table size will increase automatically when adding a new item will cause an out-of-bound error
  subroutine initiateTable(table,initialSize)
    type(LookupTable),intent(inout)::table
    integer,intent(in)::initialSize

    integer::jerr

    table%size=0
    allocate(table%list(1:initialSize),stat=jerr)
    if (jerr.ne.0) call err_exit(TRIM(__FILE__)//":"//integer_to_string(__LINE__))
  end subroutine initiateTable

  !> \brief Add item to table, ignore if exist, increase table size by 2 if table is full and expand is set to .true., signal an error otherwise
  integer function addToTable(table,item,expand)
    type(LookupTable),intent(inout)::table
    integer,intent(in)::item
    logical,intent(in),optional::expand

    integer::i
    logical::lexpand

    if (present(expand)) then
       lexpand=expand
    else
       lexpand=.false.
    end if

    do i=1,table%size
       if (table%list(i).eq.item) exit
    end do

    if (i.gt.table%size) then
       ! item not found
       if (i.gt.size(table%list)) then
          if (lexpand) then
             call reallocate(table%list,1,2*size(table%list))
          else
             call err_exit(TRIM(__FILE__)//":"//integer_to_string(__LINE__))
          end if
       end if
       table%size=i
       table%list(i)=item
    end if
    addToTable=i
  end function addToTable

  !> \brief look up item in table
  !>
  !> Return value: index of item in table; 0 if not found
  integer function indexOf(table,item)
    type(LookupTable),intent(inout)::table
    integer,intent(in)::item

    integer::i

    do i=1,table%size
       if (table%list(i).eq.item) exit
    end do
    if (i.gt.table%size) i=0
    indexOf=i
  end function indexOf
end module util_search
