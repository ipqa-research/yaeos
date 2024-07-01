module test_fitting
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("FitKijLij", test_fit_kij_lij) &
                   ]
   end subroutine collect_suite

   subroutine test_fit_kij_lij(error)
      use yaeos__fitting, only: optimize, error_function
      use yaeos__fitting_fit_kij_lij, only: FitKijLij
      use yaeos__models, only: ArModel, SoaveRedlichKwong
      use yaeos__equilibria, only: EquilibriaState
      type(error_type), allocatable, intent(out) :: error
      class(ArModel), allocatable :: model
      type(EquilibriaState) :: exp_point

      real(pr) :: Tc(2) = [126.2, 568.7]
      real(pr) :: pc(2) = [33.98, 24.90]
      real(pr) :: w(2) = [3.7e-2_pr, 0.397_pr]
      real(pr) :: X(2) = [0, 0]
      real(pr) :: err0, err_kij, err_kij_lij

      type(FitKijLij) :: fitting_problem

      exp_point = EquilibriaState( &
                  kind="bubble", T=344.5_pr, P=23.9_pr, &
                  x=[0.0309_pr, 1 - 0.0309_pr], y=[0.9883_pr, 1 - 0.9883_pr], &
                  Vx=0._pr, Vy=0._pr, beta=0.0_pr &
                  )

      fitting_problem%model = SoaveRedlichKwong(tc, pc, w)
      fitting_problem%experimental_points = [exp_point]
      fitting_problem%parameter_step = [0.1, 0.1] 

      fitting_problem%fit_kij = .true.
      fitting_problem%verbose = .false.

      err0 = error_function(x, func_data=fitting_problem)
      err_kij = optimize(X, fitting_problem)

      call check(error, err_kij < err0)
      
      fitting_problem%fit_lij = .true.
      X = 0
      err_kij_lij = optimize(X, fitting_problem)

      call check(error, err_kij_lij < err_kij)
   end subroutine
end module test_fitting
