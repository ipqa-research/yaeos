module testing_aux
   implicit none

contains

   function test_title(str)
      use stdlib_ansi, only: style_underline, style_reset, operator(//)
      character(*), intent(in) :: str
      character(100) :: test_title
      test_title = style_underline // "TESTING: " // str // style_reset
   end function

   function test_ok(str)
      use stdlib_ansi, only: fg_color_green, style_reset, operator(//)
      character(*), intent(in) :: str
      character(100) :: test_ok
      test_ok = fg_color_green // "OK: " // str // style_reset
   end function
   
   function test_notok(str)
      use stdlib_ansi, only: fg_color_red, style_reset, operator(//)
      character(*), intent(in) :: str
      character(100) :: test_notok
      test_notok = fg_color_red // "ERROR: " // str // style_reset
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