module Strings
  !Auxiliary string processing procedures.

  implicit none

  private

  public :: upcase, downcase, count_multispaces, itoa, tokenize

contains

  function upcase(string) result(res)
    character(*),intent(in) :: string
    character(len(string))  :: res
    integer :: j

    do j = 1,len(string)
      if(string(j:j) >= "a" .and. string(j:j) <= "z") then
           res(j:j) = achar(iachar(string(j:j)) - 32)
      else
           res(j:j) = string(j:j)
      end if
    end do
  end function upcase

  function downcase(string) result(res)
    character(*),intent(in) :: string
    character(len(string))  :: res
    integer :: j

    do j = 1,len(string)
      if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
           res(j:j) = achar(iachar(string(j:j)) + 32)
      else
           res(j:j) = string(j:j)
      end if
    end do
  end function downcase

  function count_multispaces(string) result(res)
    !count groups of spaces between groups of other characters
    integer :: res
    character(*),intent(in) :: string
    character(len(string))  :: loc_string
    integer :: i
    logical :: lastspace

    loc_string = adjustl(string)

    lastspace = .false.
    res = 0
    do i=1,len_trim(loc_string)
      if (loc_string(i:i) == ' ') then
        if (.not.lastspace) res = res + 1
        lastspace = .true.
      else
        lastspace = .false.
      end if
    end do

  end function count_multispaces
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function


  subroutine tokenize(stream,res,ierr)
    use str_lists, only: str_list => list, char_len
    type(str_list),intent(out) :: res
    character(*),intent(in),target :: stream
    character,pointer :: s
    integer, intent(out) :: ierr
    character(char_len) :: token
    integer :: stream_len, i, pos
    logical :: in_string

    character, parameter :: tab = achar(9)

    in_string = .false.
    pos = 1
    token = ''
    
    i = 0

    stream_len = len_trim(stream)

    do

      i = i + 1
      if (i > stream_len) exit

      s => stream(i:i)

      if (in_string) then
        !we are inside a string opened previously by "
        
        ! " finishes the string and is not counted to it if not escaped by '\'
        if (s=='"'.and.stream(i-1:i-1)/='\') then
          token(pos:pos) = s
          call send
          in_string = .false.
        ! the token length exceeded the fixed limit (the limit exists due to gfortran limitations)
        else if (pos > len(token)) then
          !TODO: remove when moving to allocatable strings
          !consider enabling strings across multiple tokens
          ierr = 2
          return
        ! add the character to the string
        else
          token(pos:pos) = s
          pos = pos + 1
        end if
        
      else
        ! parantheses are tokens of length 1
        if (s=='('.or.s==')'.or.s=='['.or.s==']') then
          if (pos>1) call res%add(token)
          token = s
          call send
        ! , and = are tokens of length 1
        else if (s==','.or.s=='=') then
          if (pos>1) call res%add(token)
          token = s
          call send
        ! " delimits a string
        else if (s=='"') then
          if (pos>1) then
            ierr = 1
            return
          else
            token(1:1) = s
            pos = pos + 1
            in_string = .true.
          end if
        ! whitespace ends a token
        else if (s==' ' .or. s==tab .or. s==new_line('a')) then
          if (pos>1) call send
        ! comments
        else if (s=='!'.or.s=='#') then
          if (pos>1) call send
          call skip_line
        ! add the character to the current token
        else
          token(pos:pos) = s
          pos = pos + 1
        end if
        
      end if

    end do

    if (pos>1) call send

    ierr = 0
    
  contains

    subroutine send
      call res%add(token)
      token = ''
      pos = 1
    end subroutine

    subroutine skip_line
      do
        if ((stream(i:i)/=new_line('a')) .and. (i<stream_len)) then
          i = i + 1
        else
          exit
        end if
      end do
    end subroutine

  end subroutine tokenize

end module Strings




module ParseTrees_Fields
  use Strings

  use str_lists, only: str_list => list, char_len

  type field_names
    character(char_len) :: name
    !not nullified bucause of a bug in GCC 5.3
    class(*), pointer :: var
  end type

  type field_names_a
    character(char_len) :: name
    !not nullified bucause of a bug in GCC 5.3
    class(*), pointer :: var(:)
  end type

  type field_names_a_int_alloc
    character(char_len) :: name
    !not nullified bucause of a bug in GCC 5.3
    integer, allocatable :: var(:)
  end type

  type field_names_str
    character(char_len) :: name
    !workaround of BUG https://gcc.gnu.org/bugzilla/show_bug.cgi?id=60359 fixed in GCC 4.9
    character(char_len), pointer :: var => null()
  end type

  interface field_names
    procedure field_names_init
  end interface

  interface field_names_a
    procedure field_names_init
  end interface

  interface field_names_a_int_alloc
    procedure field_names_a_int_alloc_init
  end interface

contains

  function field_names_init(name, var) result(res)
    type(field_names) :: res
    character(*) :: name
    class(*), target, intent(in) :: var

    res%name = name
    res%var => var
  end function

  function field_names_a_init(name, var) result(res)
    type(field_names_a) :: res
    character(*) :: name
    class(*), target, intent(in) :: var(:)

    res%name = name
    res%var => var
  end function

  function field_names_a_int_alloc_init(name) result(res)
    type(field_names_a_int_alloc) :: res
    character(*) :: name

    res%name = name
  end function


end module ParseTrees_Fields




module ParseTrees

  use ParseTrees_Fields

  implicit none

  !HACK: private must be disabled due to a bug in Intel Fortran 2021.1
  !   private

  public field_names, field_names_a, field_names_a_int_alloc, field_names_str
  !HACK: must be called directly, the overloaded constructor causes https://gcc.gnu.org/bugzilla/show_bug.cgi?id=49213
  
  public field_names_init, field_names_a_init, field_names_a_int_alloc_init

  public tree_object, tree_object_field, tree_object_fields, tree_object_ptr, &
         parse_file, find_object_get_field_values, get_object_field_values, &
         char_len

  type tree_object_field
    character(char_len) :: name
    logical :: is_array = .false.
    logical :: is_object = .false.
    character(char_len) :: value
    type(tree_object), pointer :: object_value => null()
    character(char_len), allocatable :: array_value(:)
  end type

  type tree_object_fields
    type(tree_object_field), allocatable :: array(:)
  end type

  type tree_object
    character(char_len) :: name
    type(tree_object_fields) :: fields
  contains
    procedure :: move_from => tree_object_move
    procedure :: finalize => tree_object_finalize
  end type
  
  type tree_object_ptr
    type(tree_object), pointer :: ptr => null()
  end type

contains

  recursive subroutine tree_object_finalize(obj)
    class(tree_object), intent(inout) :: obj
    integer :: i

    if (allocated(obj%fields%array)) then
      do i = 1, size(obj%fields%array)
        associate(field => obj%fields%array(i))
          if (associated(field%object_value)) then
            call field%object_value%finalize
            deallocate(field%object_value)
          end if
        end associate
      end do
      deallocate(obj%fields%array)
    end if
  end subroutine


  recursive subroutine tree_object_move(lhs, rhs)
    !shallow copy
    class(tree_object), intent(out)   :: lhs
    type(tree_object), intent(inout) :: rhs
    
    lhs%name = rhs%name
    call move_alloc(rhs%fields%array, lhs%fields%array)
  end subroutine


  subroutine to_array(ch_list, ch_array)
    type(str_list) :: ch_list
    character(char_len), allocatable :: ch_array(:)
    allocate(ch_array(0))
    call ch_list%for_each(aux)
  contains
    subroutine aux(item)
      character(char_len) :: item

      ch_array = [ch_array, item]
    end subroutine
  end subroutine

  subroutine parse_file(tree, fname, stat)
    use Strings, only: itoa, downcase
    type(tree_object), allocatable, intent(out) :: tree(:)
    character(*), intent(in) :: fname
    integer, intent(out) :: stat

    character(:), allocatable :: stream
    type(tree_object), allocatable :: tmp(:)
    type(str_list) :: token_list
    character(char_len), allocatable :: tokens(:)
    type(tree_object) :: object

    integer :: pos, io, unit
    character(256) :: line


    open(newunit=unit, file=fname,status="old",action="read",iostat=stat)
       
    if (stat/=0) return

    stream = ''

    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      stream = stream // ' ' // trim(adjustl(line)) // new_line('a')
    end do

    call tokenize(stream, token_list, stat)

    if (stat/=0) then
      select case (stat)
        case (2)
          write(*,*) "Error, unterminated string or string too long in file ", trim(fname)
      end select

      call token_list%finalize
      return
    end if

    call to_array(token_list, tokens)
    call token_list%finalize

    allocate(tree(0))

    if (size(tokens)<=0) then
      !the file is empty or contains just comments
      stat = 0
      return
    end if
    
    pos = 1
    do
      call get_object(object, pos, stat)
      if (stat /=0) exit

      call move_alloc(tree, tmp)
      allocate(tree(size(tmp)+1))
      tree(:size(tmp)) = tmp

      call tree(size(tmp)+1)%move_from(object)

      if (pos > size(tokens)) exit
    end do

  contains

    recursive subroutine get_array(array, pos, stat)
      character(char_len), allocatable :: array(:)
      integer, intent(inout) :: pos
      integer, intent(out) :: stat
      character(char_len) :: str

      if (pos+1 > size(tokens)) then
        !need 2 tokens for an empty array
        stat = 4
        return
      end if

      allocate(array(0))
      pos = pos + 1
      do
        if (tokens(pos) == ']') exit

        call get_string(str, pos, stat)

        array = [array, str]

        if (stat /= 0) return

        if (tokens(pos) /= ']' .and. tokens(pos) /= ',') then
          write(*,*) "Error in '" // &
                     fname // &
                     "', expected '" // &
                     ",' or ']" // &
                     "' read '" // &
                     downcase(tokens(pos)) // &
                     "' instead."
          stat = 1
        else if (tokens(pos) == ',') then
          pos = pos + 1
        end if
        if (stat /= 0) return
      end do

      call check_string(']', pos, stat)

    end subroutine

    recursive subroutine get_object(object, pos, stat)
      type(tree_object) :: object
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      if (pos+1 > size(tokens)) then
        !need 2 tokens for an empty anonymous object
        stat = 2
        return
      end if

      if (tokens(pos)=='(') then
        !anonymous object
        object%name = ""
        pos = pos + 1
        call get_object_fields(object%fields, pos, stat)
        if (stat /= 0) return
        call check_string(')', pos, stat)
      else if (tokens(pos+1)=='(') then
        call get_string(object%name, pos, stat)
        if (stat /= 0) return
        pos = pos + 1 !checked above
        call get_object_fields(object%fields, pos, stat)
        if (stat /= 0) return
        call check_string(')', pos, stat)
      else
        stat = 3
        return
      end if
    end subroutine

    subroutine check_string(str, pos, stat)
      character(*), intent(in) :: str
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      if (downcase(str)==downcase(tokens(pos))) then
        stat = 0
        pos = pos + 1
      else
        write(*,*) "Error in '" // &
                   fname // &
                   "', expected '" // &
                   str // &
                   "' read '" // &
                   downcase(tokens(pos)) // &
                   "' instead."
        stat = 1
      end if
    end subroutine

    recursive subroutine get_object_fields(fields, pos, stat)
      type(tree_object_fields) :: fields
      integer, intent(inout) :: pos
      integer, intent(out) :: stat
      type(tree_object_field), allocatable :: tmp(:)

      allocate(fields%array(0))
      do
        if (tokens(pos) == ')') return

        call move_alloc(fields%array, tmp)
        allocate(fields%array(size(tmp)+1))
        fields%array(:size(tmp)) = tmp

        call get_field(fields%array(size(tmp)+1), pos, stat)
        if (stat /= 0) return

        if (tokens(pos) /= ')' .and. tokens(pos) /= ',') then
          write(*,*) "Error in '" // &
                     fname // &
                     "', expected '" // &
                     ",' or ')" // &
                     "' read '" // &
                     downcase(tokens(pos)) // &
                     "' instead."
          stat = 1
        else if (tokens(pos) == ',') then
          pos = pos + 1
        end if
        if (stat /= 0) return

      end do
    end subroutine

    recursive subroutine get_field(field, pos, stat)
      type(tree_object_field) :: field
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      call get_string(field%name, pos, stat)
      call check_string('=', pos, stat)

      if (tokens(pos) == '(' .or. tokens(pos+1) == '(') then
        field%is_object = .true.
        allocate(field%object_value)
        call get_object(field%object_value, pos, stat)
      else if (tokens(pos) == '[') then
        field%is_array = .true.
        call get_array(field%array_value, pos, stat)
      else
        call get_string(field%value, pos, stat)
      end if
    end subroutine

    recursive subroutine get_string(str, pos, stat)
      character(char_len) :: str
      integer, intent(inout) :: pos
      integer, intent(out) :: stat
      if (pos>size(tokens)) then
        stat = 1
        return
      else
        str = tokens(pos)
        pos = pos + 1
      end if
    end subroutine

  end subroutine parse_file
  
  
  subroutine find_object_get_field_values(tree, object_name, stat, &
                                     fields, &
                                     fields_a, &
                                     fields_a_int_alloc, &
                                     fields_str)
    type(tree_object), intent(in) :: tree(:)
    character(*), intent(in) :: object_name
    integer, intent(out) :: stat
    !stat 0    .. success
    !stat < 0  .. multiple definitions of object named object_name
    !stat 1    .. object named object_name not found
    !stat 2    .. inconsistent number of components in an array in fields_a
    !stat 11   .. unexpected type of variable in fields
    !stat 12   .. unexpected type of variable in fields_a
    type(field_names), intent(inout), optional :: fields(:)
    type(field_names_a), intent(inout), optional :: fields_a(:)
    type(field_names_a_int_alloc), intent(inout), optional :: fields_a_int_alloc(:)
    type(field_names_str), intent(inout), optional :: fields_str(:)

    integer :: iobj

    stat = 1

    do iobj = 1, size(tree)
    
      if (downcase(tree(iobj)%name)==object_name) then
        stat = stat - 1
        call get_object_field_values(tree(iobj), stat, &
                                     fields, &
                                     fields_a, &
                                     fields_a_int_alloc, &
                                     fields_str)
        if (stat>1) return
      end if
      
    end do

  end subroutine find_object_get_field_values




  subroutine get_object_field_values(obj, stat, &
                        fields, &
                        fields_a, &
                        fields_a_int_alloc, &
                        fields_str)
    !To extract values of variables from the parse tree
    use iso_fortran_env, only: real32, real64, int32, int64
    type(tree_object), intent(in) :: obj
    integer, intent(out) :: stat
    !stat 0    .. success
    !stat < 0  .. multiple definitions of object named object_name
    !stat 2    .. inconsistent number of components in an array in fields_a
    !stat 11   .. unexpected type of variable in fields
    !stat 12   .. unexpected type of variable in fields_a
    type(field_names), intent(inout), optional :: fields(:)
    type(field_names_a), intent(inout), optional :: fields_a(:)
    type(field_names_a_int_alloc), intent(inout), optional :: fields_a_int_alloc(:)
    type(field_names_str), intent(inout), optional :: fields_str(:)

    integer :: i, j


    if (allocated(obj%fields%array)) then

      associate(obj_fields => obj%fields%array)

fields_do:  do j = 1, size(obj_fields)
         
          !workaround of 60359, see also below
          if (present(fields_str)) then
            do i = 1, size(fields_str)
              if (obj_fields(j)%name == fields_str(i)%name) then                
                read(obj_fields(j)%value, *) fields_str(i)%var
                cycle fields_do
              end if
            end do
          end if


          if (present(fields)) then
            do i = 1, size(fields)
              if (obj_fields(j)%name == fields(i)%name) then
                select type (var => fields(i)%var)
                  type is (tree_object_ptr)
                    if (.not.obj_fields(j)%is_object .or. &
                        .not.associated(obj_fields(j)%object_value)) then
                      stat = 13
                      return
                    end if
                    var%ptr => obj_fields(j)%object_value
                  type is (integer)
                    read(obj_fields(j)%value, *) var
                  type is (real(real32))
                    read(obj_fields(j)%value, *) var
                  type is (real(real64))
                    read(obj_fields(j)%value, *) var
                  type is (logical)
                    read(obj_fields(j)%value, *) var
                  type is (character(*))
                    read(obj_fields(j)%value, *) var
                  class default
                    stat = 11
                    return
                end select
                cycle fields_do
              end if
            end do
          end if

         if (present(fields_a)) then
           do i = 1, size(fields_a)
              if (obj_fields(j)%name == fields_a(i)%name) then
                if (.not.obj_fields(j)%is_array .or. &
                    .not.allocated(obj_fields(j)%array_value)) then
  ! uncomment to get details when debugging
                  write(*,*) "Error, expecting an array in '",trim(fields_a(i)%name),"'."
                  stat = 2
                  return
                end if
                if (size(fields_a(i)%var)/=size(obj_fields(j)%array_value)) then
  ! uncomment to get details when debugging
                  write(*,*) "Error, expecting", &
                             size(fields_a(i)%var), &
                             "vector components, but", &
                             size(obj_fields(j)%array_value), &
                             "components present."
                  stat = 3
                  return
                end if

                !workaround of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90498
                ! when fixes in the wild, use just
                ! select type (var => fields_a(i)%var)
                block
                  class(*), pointer :: var(:)
                  var => fields_a(i)%var
                  
                  select type (var)
                    type is (integer)
                      read(obj_fields(j)%array_value, *) var
                    type is (real(real32))
                      read(obj_fields(j)%array_value, *) var
                    type is (real(real64))
                      read(obj_fields(j)%array_value, *) var
                    type is (logical)
                      read(obj_fields(j)%array_value, *) var
                    class default
                      stat = 12
                      return
                  end select
                end block
                
                cycle fields_do
              end if
            end do
          end if

          
          if (present(fields_a_int_alloc)) then
            do i = 1, size(fields_a_int_alloc)
              if (obj_fields(j)%name == fields_a_int_alloc(i)%name) then
                allocate(fields_a_int_alloc(i)%var(size(obj_fields(j)%array_value)))
                if (size(fields_a_int_alloc(i)%var)>0) &
                  read(obj_fields(j)%array_value, *) fields_a_int_alloc(i)%var
              end if
            end do
          end if

        end do fields_do

      end associate

    end if
  end subroutine get_object_field_values
  


end module ParseTrees
