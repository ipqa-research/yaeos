module yaeos__tapenade_interfaces
    use yaeos__constants, only: pr
    implicit none

   interface 
      subroutine pushinteger4(i)
         integer :: i
      end subroutine
      
      subroutine popinteger4(i)
         integer :: i
      end subroutine

      subroutine pushreal8array(a, n)
         import pr
         real(pr), dimension(n) :: a
         integer :: n
      end subroutine

      subroutine pushreal8(a)
         import pr
         real(pr) :: a
      end subroutine

      subroutine POPREAL8(a)
         import pr
         real(pr) :: a
      end subroutine

      subroutine POPREAL8ARRAY(a, n)
         import pr
         real(pr), dimension(n) :: a
         integer :: n
      end subroutine
   end interface
end module