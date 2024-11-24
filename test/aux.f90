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
end module