program main
   !! Binary system parameter optimization
   use yaeos, only: EquilibriumState, pr, ArModel, SoaveRedlichKwong, CubicEoS, saturation_pressure, PengRobinson78
   use forsus, only: Substance, forsus_dir
   use yaeos__fitting, only: FittingProblem, error_function, optimize
   use yaeos__optimizers_powell_wrap, only: PowellWrapper
   integer, parameter :: nc = 2
   integer :: i, infile, iostat

   type(EquilibriumState), allocatable :: exp_points(:)
   type(EquilibriumState) :: point

   type(Substance) :: sus(nc)

   type(PowellWrapper) :: opt

   class(ArModel), allocatable :: model

   real(pr) :: T, P, x1, y1, told, error
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
         point = EquilibriumState( &
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

   ! First set up the model to be used
   model = PengRobinson78(&
      sus%critical%critical_temperature%value, &
      sus%critical%critical_pressure%value/1e5, &
      sus%critical%acentric_factor%value &
      )

   call fit_kij_lij(model, exp_points)

   call exit

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

contains

   subroutine fit_kij_lij(model, exp_points)
      use yaeos__fitting_fit_kij_lij, only: FitKijLij
      class(ArModel), intent(in) :: model
      type(EquilibriumState), intent(in) :: exp_points(:)
      type(FitKijLij) :: prob
      integer, parameter :: np=2
      real(pr) :: X(np)
      ! Set up the optimization problem settings
      prob = FitKijLij(&
         model=model, experimental_points=exp_points, &
         verbose=.true., fit_kij=.true., fit_lij=.true. &
      )
      
      ! Set up the experimental points
      prob%experimental_points = exp_points

      ! Fit Kij and Lij
      prob%fit_kij = .true.
      prob%fit_lij = .true.
      prob%verbose = .true.

      ! Initial X value
      X = 0.01
      opt%parameter_step = [0.01, 0.01]
      print *, "X0:", X
      error = optimize(X, opt, prob)
      print *, "FO:", error
      print *, "Xf:", X

      call prob%get_model_from_X(X)
   end subroutine fit_kij_lij
end program main
