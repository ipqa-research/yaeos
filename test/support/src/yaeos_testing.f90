module testing_aux
   implicit none

   character(len=*), parameter :: red_start = char(27) // "[31m"
   character(len=*), parameter :: green_start = char(27) // "[32m"
   character(len=*), parameter :: underline_start = char(27) // "[4m"

   character(len=*), parameter :: style_reset = char(27) // "[0m"

   logical :: compact = .true.

contains

   function test_title(str)
      character(*), intent(in) :: str
      character(100) :: test_title
      test_title = new_line("a") // underline_start // "TESTING: " // str // style_reset
   end function test_title

   function test_ok(str)
      character(*), intent(in) :: str
      character(100) :: test_ok
      if (compact) then
         test_ok = green_start // "." // style_reset
      else
         test_ok = green_start // "OK: " // str // style_reset
      end if
   end function test_ok

   function test_notok(str)
      character(*), intent(in) :: str
      character(100) :: test_notok
      test_notok = red_start // "ERROR: " // str // style_reset
   end function test_notok

   subroutine assert(expr, test_name)
      logical, intent(in) :: expr
      character(*), intent(in) :: test_name
      if (.not. expr) then
         error stop test_notok(test_name)
      else
         if (compact) then
            write(*,  "(A)", advance="no") trim(test_ok(test_name))
         else
            write(*,  *) trim(test_ok(test_name))
         end if
      end if
   end subroutine assert
end module testing_aux
