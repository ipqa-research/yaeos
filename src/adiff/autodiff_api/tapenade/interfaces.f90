module yaeos__tapenade_interfaces
    !! # Tapenade runtime stack interfaces
    !! Fortran interfaces to the push/pop routines used by Tapenade-generated adjoint code.
    !!
    !! # Description
    !! These external procedures exchange scalar and array values with the
    !! Tapenade AD stack. They are implementation support for generated
    !! derivative routines and are not normally called directly by model users.
    !!
    !! # Examples
    !! The interfaces are invoked from Tapenade-generated code during reverse
    !! mode evaluation.
    use yaeos__constants, only: pr
    implicit none

   interface
      subroutine pushinteger4(i)
         !! Push a 4-byte integer value onto the AD stack.
         integer :: i !! Integer value to push.
      end subroutine
      
      subroutine popinteger4(i)
         !! Pop a 4-byte integer value from the AD stack.
         integer :: i !! Returned integer value.
      end subroutine

      subroutine pushreal8array(a, n)
         !! # pushreal8array
         !!
         !! # Description
         !! Push a real array onto the AD stack.
         import pr
         real(pr), dimension(n) :: a !! Array values to push.
         integer :: n !! Number of array elements.
      end subroutine

      subroutine pushreal8(a)
         !! # pushreal8
         !!
         !! # Description
         !! Push one real value onto the AD stack.
         import pr
         real(pr) :: a !! Value to push.
      end subroutine

      subroutine POPREAL8(a)
         !! # POPREAL8
         !!
         !! # Description
         !! Pop one real value from the AD stack.
         import pr
         real(pr) :: a !! Returned value.
      end subroutine

      subroutine POPREAL8ARRAY(a, n)
         !! # POPREAL8ARRAY
         !!
         !! # Description
         !! Pop a real array from the AD stack.
         import pr
         real(pr), dimension(n) :: a !! Returned array values.
         integer :: n !! Number of array elements.
      end subroutine
   end interface
end module
