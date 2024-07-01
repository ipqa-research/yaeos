module yaeos__fitting_fit_nrtl
   use yaeos__constants, only: pr
   use yaeos__fitting, only: FittingProblem
   use yaeos__models, only: ArModel, NRTL, CubicEoS, MHV
   use forsus, only: Substance
   implicit none

   integer, parameter :: nc = 2

   type, extends(FittingProblem) :: FitMHVNRTL
   contains
      procedure :: get_model_from_X => model_from_X
   end type

contains

   subroutine init_model(problem, sus)
      use yaeos, only: R, ArModel, CubicEoS, PengRobinson78, RKPR, SoaveRedlichKwong
      class(FitMHVNRTL), intent(in out) :: problem
      type(Substance), intent(in) :: sus(2)
      type(MHV) :: mixrule
      type(NRTL) :: ge
      real(pr) :: tc(nc), pc(nc), w(nc), vc(nc), zc(nc)
      real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc), bi(nc)

      a=0; b=0; c=0

      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5
      w = sus%critical%acentric_factor%value
      vc = sus%critical%critical_volume%value
      zc = pc*vc/(R*tc)

      ge = NRTL(a, b, c)

      ! problem%model = RKPR(tc, pc, w, zc)
      allocate(CubicEoS :: problem%model)
      problem%model = SoaveRedlichKwong(tc, pc, w)

      associate(m => problem%model)
         select type(m)
         type is (CubicEoS)
            bi = m%b
            mixrule = MHV(ge=ge, q=-0.593_pr, b=bi)
            deallocate(m%mixrule)
            m%mixrule = mixrule
         end select
      end associate
   end subroutine init_model

   function model_from_X(problem, X) result(model)
      use yaeos, only: R, RKPR, PengRobinson78, ArModel, QMR, CubicEoS
      use yaeos__models_ar_cubic_quadratic_mixing, only: RKPR_D1mix
      real(pr), intent(in) :: X(:)
      class(FitMHVNRTL), intent(in) :: problem
      class(ArModel), allocatable :: model
      type(NRTL) :: ge

      real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc)

      a=0; b=0; c=0

      a(1, 2) = x(1)
      a(2, 1) = x(2)

      b(1, 2) = x(3)
      b(2, 1) = x(4)

      c(1, 2) = x(5)
      c(2, 1) = x(6)

      ge = NRTL(a, b, c)

      associate (pm => problem%model)
         select type(pm)
         type is (CubicEoS)
            model = pm
         end select
      end associate

      ! model = problem%model

      select type(model)
      class is (CubicEoS)
         associate(mr => model%mixrule)
            select type (mr)
            class is (MHV)
               mr%l(1, 2) = x(7)
               mr%l(2, 1) = x(7)
               mr%ge = ge
               model%del1 = x(8:)
               model%del2 = (1._pr - model%del1)/(1._pr + model%del1)
            end select
         end associate
      end select
   end function model_from_X
end module

program main
   !! Binary system parameter optimization
   use yaeos, only: EquilibriaState, pr, ArModel, PengRobinson78, CubicEoS, saturation_pressure
   use forsus, only: Substance, forsus_dir
   use yaeos__fitting, only: FittingProblem, fobj, optimize
   use yaeos__fitting_fit_nrtl, only: FitMHVNRTL, init_model
   integer, parameter :: nc = 2, np=7 + nc
   integer :: i, infile, iostat

   type(EquilibriaState), allocatable :: exp_points(:)
   type(EquilibriaState) :: point

   type(FitMHVNRTL) :: prob
   type(Substance) :: sus(2)

   class(ArModel), allocatable :: model

   real(pr) :: T, P, x1, y1, kij, X(np), told, error
   character(len=14) :: kind

   ! ==========================================================================
   ! Setup components and read data file
   ! --------------------------------------------------------------------------
   forsus_dir = "build/dependencies/forsus/data/json"
   sus(1) = Substance("nitrogen", only=["critical"])
   sus(2) = Substance("n-octane", only=["critical"])

   allocate (exp_points(0))
   open (newunit=infile, file="fit_case", iostat=iostat)
   do
      read (infile, *, iostat=iostat) kind, t, p, x1, y1
      if (iostat /= 0) exit
      select case (kind)
      case ("bubble", "dew", "liquid-liquid")
         point = EquilibriaState( &
                 kind=kind, T=T, P=P, x=[x1, 1 - x1], y=[y1, 1 - y1], &
                 Vx=0._pr, Vy=0._pr, iters=0, beta=0._pr &
                 )
      end select
      exp_points = [exp_points, point]
   end do
   close (infile)

   ! ==========================================================================
   ! Setup optimization problem and call the optimization function
   ! --------------------------------------------------------------------------
   call init_model(prob, sus)

   prob%experimental_points = exp_points

   X = 0
   X(1:2) = [0.1, 0.3]
   X(5:6) = [0.1, 0.2]
   X(8:) = [1, 1]

   prob%parameter_step = [(0.5_pr, i=1,size(x))]
   prob%solver_tolerance = 1e-7
   prob%verbose = .true.

   print *, "X0:", X
   error = optimize(X, prob)
   print *, "FO:", error
   print *, "Xf:", X

   if (allocated(model)) deallocate (model)
   model = prob%get_model_from_X(X)

   ! ===========================================================================
   ! Write out results and experimental values
   ! ---------------------------------------------------------------------------
   told = exp_points(1)%T
   do i = 1, size(exp_points)
      point = saturation_pressure( &
              model, exp_points(i)%x, exp_points(i)%t, kind="bubble", &
              p0=exp_points(i)%p, y0=exp_points(i)%y &
              )
      if (told /= point%t) write (*, "(/)")
      print *, exp_points(i)%x(1), exp_points(i)%y(1), exp_points(i)%P, &
         point%x(1), point%y(1), point%P
      told = point%T
   end do
end program main
