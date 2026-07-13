module tests_cmr
   use yaeos
   use testing_aux, only: assert
   integer, parameter :: nc = 3
contains
   subroutine test_aijk_numdiff
      use yaeos__models_ar_cubic_mixing_base, only: CMR_aijk
      use yaeos__extra_fluids, only: co2_h2o_isop, FixtureFluid
      real(pr) :: k(nc, nc, nc), dkdt(nc, nc, nc), dkdt2(nc, nc, nc)
      type(FixtureFluid) :: fluid

      integer :: i, j, l
      real(pr) :: T, dt=1e-3
      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
      real(pr) :: df, df2

      fluid = co2_h2o_isop()

      k(1, 1, 2) = 0.1
      k(1, 2, 1) = 0.1
      k(2, 1, 1) = 0.1

      k(1, 2, 3) = 0.3
      k(2, 1, 3) = 0.3
      k(3, 2, 1) = 0.3

      k(2, 2, 1) = 0.2
      k(2, 1, 2) = 0.2
      k(1, 2, 2) = 0.2
      dkdt = 0
      dkdt2 = 0

      T = 200

      associate(model => fluid%ar_model)
         select type(model)
          type is (CubicEoS)
            call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
            call CMR_aijk(&
               ai, daidt, daidt2, k, dkdt, dkdt2, &
               a, dadt, dadt2 &
               )
         end select
      end associate

      do i=1,nc
         do j=1,nc
            do l=1,nc
               df = (f(T+dT, i, j, l) - f(T-dT, i, j, l))/(2*dt)

               df2 = (f(T+dT, i, j, l) - 2*f(T, i, j, l) + f(T-dT, i, j, l))/(dT**2)
               call assert(abs((df - dadt(i, j, l))/df) < 1e-5, "First order analytic and numerical derivatives should be similar")
               call assert(abs((df2 - dadt2(i, j, l))/df) < 1e-5, "Second order analytic and numerical derivatives should be similar")
            end do
         end do
      end do

   contains

      real(pr) function f(T, i, j, l)
         real(pr), intent(in) :: T
         integer, intent(in) :: i, j, l
         real(pr) :: ai(nc), daidt(nc), daidt2(nc)
         real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
         real(pr) :: eps=1e-5_pr

         associate(model => fluid%ar_model)
            select type(model)
             type is (CubicEoS)
               call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
               call CMR_aijk(&
                  ai, daidt, daidt2, k, dkdt, dkdt2, &
                  a, dadt, dadt2 &
                  )
            end select
         end associate
         f = a(i, j, l)
      end function f
   end subroutine test_aijk_numdiff
end module tests_cmr

program test_cmr
   use tests_cmr, only: test_aijk_numdiff
   use testing_aux, only: test_title

   print *, test_title("Cubic Mixing Rules")
   call test_aijk_numdiff
end program test_cmr
