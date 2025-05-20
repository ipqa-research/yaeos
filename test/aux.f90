module testing_aux
   implicit none

   character(len=*), parameter :: red_start = char(27) // "[31m"
   character(len=*), parameter :: green_start = char(27) // "[32m"
   character(len=*), parameter :: underline_start = char(27) // "[4m"

   character(len=*), parameter :: style_reset = char(27) // "[0m"

contains

   function test_title(str)
      character(*), intent(in) :: str
      character(100) :: test_title
      test_title = underline_start // "TESTING: " // str // style_reset
   end function

   function test_ok(str)
      character(*), intent(in) :: str
      character(100) :: test_ok
      test_ok = green_start // "OK: " // str // style_reset
   end function
   
   function test_notok(str)
      character(*), intent(in) :: str
      character(100) :: test_notok
      test_notok = red_start // "ERROR: " // str // style_reset
   end function

   subroutine assert(expr, test_name)
      logical, intent(in) :: expr
      character(*), intent(in) :: test_name
      if (.not. expr) then
         error stop test_notok(test_name)
      else
         print *, test_ok(test_name)
      end if
   end subroutine
end module