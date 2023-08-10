module Pthreads
  implicit none

  interface
    subroutine pthread_create_opaque(threadptr, procptr, dataptr, err) bind(C,name="pthread_create_opaque")
      use iso_c_binding
      type(c_ptr) :: threadptr
      type(c_funptr),value :: procptr
      type(c_ptr),value :: dataptr
      integer(c_int),intent(out) :: err
    end subroutine

    subroutine pthread_join_opaque(thread, err) bind(C,name="pthread_join_opaque")
      use iso_c_binding
      type(c_ptr),value :: thread
      integer(c_int),intent(out) :: err
    end subroutine
  end interface
end module Pthreads
