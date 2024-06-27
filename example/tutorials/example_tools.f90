! Here there are some functions defined to simplify some examples
module yaeos__example_tools
   use yaeos

contains

   type(CubicEoS) function methane_butane_pr76() result(model)
      !! Binary mixture of methane/butane with PR76
      integer, parameter :: nc = 2
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)

      integer :: i

      ! Methane/ Butane mixture
      tc = [190.564, 425.12]     ! Critical temperatures
      pc = [45.99, 37.96]        ! Critical pressures
      w = [0.0115478, 0.200164]  ! Acentric factors

      model = PengRobinson76(tc, pc, w)
   end function methane_butane_pr76
end module yaeos__example_tools
