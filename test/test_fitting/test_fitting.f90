module test_fitting
   use yaeos, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose, rel_error
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [&
                     new_unittest("FitKijLij", test_fit_kij_lij), &
                     new_unittest("FitMHVNRTL", test_fit_MHV_NRTL) &
                   ]
   end subroutine collect_suite

   subroutine test_fit_kij_lij(error)
      use yaeos__fitting, only: optimize, error_function
      use yaeos__fitting_fit_kij_lij, only: FitKijLij
      use yaeos__models, only: ArModel, SoaveRedlichKwong
      use yaeos__equilibria, only: EquilibriumState
      use yaeos__optimizers_powell_wrap, only: PowellWrapper
      type(error_type), allocatable, intent(out) :: error
      class(ArModel), allocatable :: model
      type(EquilibriumState) :: exp_points
      type(PowellWrapper) :: opt

      real(pr) :: Tc(2) = [126.2, 568.7]
      real(pr) :: pc(2) = [33.98, 24.90]
      real(pr) :: w(2) = [3.7e-2_pr, 0.397_pr]
      real(pr) :: X(2) = [0, 0]
      real(pr) :: err0, err_kij, err_kij_lij

      type(FitKijLij) :: fitting_problem

      exp_points = &
            EquilibriumState( &
                  kind="bubble", T=344.5_pr, P=23.9_pr, &
                  x=[0.0309_pr, 1 - 0.0309_pr], y=[0.9883_pr, 1 - 0.9883_pr], &
                  Vx=0._pr, Vy=0._pr, beta=0.0_pr &
                  )

      fitting_problem%model = SoaveRedlichKwong(tc, pc, w)
      fitting_problem%experimental_points = [exp_points]
      opt%parameter_step = [0.1, 0.1] 

      fitting_problem%fit_kij = .true.
      fitting_problem%verbose = .false.

      call error_function(x, err0, func_data=fitting_problem)
      err_kij = optimize(X, opt, fitting_problem)

      call check(error, err_kij < err0)
      
      fitting_problem%fit_lij = .true.
      fitting_problem%verbose = .true.
      X = 0
      err_kij_lij = optimize(X, opt, fitting_problem)
      call check(error, err_kij_lij < err_kij)
   end subroutine

   subroutine test_fit_mhv_nrtl(error)
      use yaeos__optimizers_powell_wrap, only: PowellWrapper
      use yaeos__optimizers_nelder_mead, only: NelderMead
      use yaeos__fitting, only: optimize, error_function
      use yaeos__fitting_fit_nrtl_mhv, only: FitMHVNRTL
      use yaeos__models, only: CubicEoS, GeModel, NRTL, SoaveRedlichKwong, MHV
      use yaeos__equilibria, only: EquilibriumState
      use yaeos, only: PTEnvel2, pt_envelope_2ph, saturation_pressure, saturation_temperature
      type(error_type), allocatable, intent(out) :: error
      type(CubicEoS) :: model
      type(NRTL) :: ge_model
      type(MHV) :: mixrule
      type(EquilibriumState) :: exp_point

      real(pr) :: Tc(2) = [647.13, 514.0]
      real(pr) :: pc(2) = [220.55, 61.37] 
      real(pr) :: w(2) = [0.344861, 0.643558]
      real(pr), allocatable :: X(:)
      real(pr) :: err0, err_lij, err_ge, err_ge_lij

      type(FitMHVNRTL) :: fitting_problem
      type(PowellWrapper) :: opt
      type(PTEnvel2) :: env

      exp_point = EquilibriumState( &
                  kind="bubble", T=437._pr, P=25.36_pr, &
                  x=[0.21_pr, 1 - 0.21_pr], y=[0.4_pr, 0.6_pr], &
                  Vx=0._pr, Vy=0._pr, beta=0.0_pr &
                  )

     
      ! Provide initials for the NRTL model
      ! reshape (a11 a12 a21 a22)
      ge_model%a=reshape([0.0, 3.458, -0.80, 0.0], [2, 2])
      ge_model%b=reshape([0.0, -586.0, 246.0, 0.0], [2, 2])
      ge_model%c=reshape([0.0, 0.3, 0.3, 0.0], [2, 2])

      model = SoaveRedlichKwong(tc, pc, w)
      mixrule = MHV(ge=ge_model, q=-0.593_pr, b=model%b)
      deallocate(model%mixrule)
      model%mixrule = mixrule
      
      fitting_problem%experimental_points = [exp_point]
      fitting_problem%verbose = .true.

      X = [3.458, -0.8, -586.0, 246.0, 0.3, 0.3, -0.15]
      opt%parameter_step = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001]
      fitting_problem%model = model
      fitting_problem%fit_lij = .true.
      call  error_function(x, err0, func_data=fitting_problem)
      err_lij = optimize(X, opt, fitting_problem)
      
      fitting_problem%verbose = .false.
      X = [3.458, -0.8, -586.0, 246.0, 0.3, 0.3, 0.0]
      fitting_problem%model = model
      fitting_problem%fit_lij = .false.
      fitting_problem%fit_nrtl = .true.
      err_ge = optimize(X, opt, fitting_problem)
      
      X = [3.458, -0.8, -586.0, 246.0, 0.3, 0.3, 0.0]
      fitting_problem%model = model
      fitting_problem%fit_lij = .true.
      fitting_problem%fit_nrtl = .true.
      err_ge_lij = optimize(X, opt, fitting_problem)

      print *, err0, err_lij, err_ge, err_ge_lij

      call check(error, err_lij < err0)
      call check(error, err_ge_lij < err_lij)
      call check(error, err_ge < err_lij)
      
   end subroutine
end module test_fitting
