
program main
   !! Binary system parameter optimization
   use yaeos, only: EquilibriaState, pr, ArModel, PengRobinson78, CubicEoS, saturation_pressure
   use forsus, only: Substance, forsus_dir
   use yaeos__fitting, only: FittingProblem, fobj, optimize
   use yaeos__fitting_fit_kij_lij, only: FitKijLij, init_model
   integer, parameter :: nc=2
   integer :: i, infile, iostat

   type(EquilibriaState), allocatable :: exp_points(:)
   type(EquilibriaState) :: point

   type(FitKijLij) :: prob
   type(Substance) :: sus(2)

   class(ArModel), allocatable :: model

   real(pr) :: T, P, x1, y1, kij, X(4), told, error
   character(len=14) :: kind

   ! ==========================================================================
   ! Setup components and read data file
   ! --------------------------------------------------------------------------
   forsus_dir = "build/dependencies/forsus/data/json"
   sus(1) = Substance("nitrogen", only=["critical"])
   sus(2) = Substance("n-octane", only=["critical"])

   allocate(exp_points(0))
   open(newunit=infile, file="fit_case", iostat=iostat)
   do
      read(infile, *, iostat=iostat) kind, t, p, x1, y1
      if (iostat /= 0) exit
      select case(kind)
       case("bubble", "dew", "liquid-liquid")
         point = EquilibriaState(&
            kind=kind, T=T, P=P, x=[x1, 1-x1], y=[y1,1-y1], &
            Vx=0._pr, Vy=0._pr, iters=0, beta=0._pr &
            )
      end select
      exp_points = [exp_points, point]
   end do
   close(infile)

   ! ==========================================================================
   ! Setup optimization problem and call the optimization function
   ! --------------------------------------------------------------------------
   call init_model(prob, sus)
   prob%fit_kij = .true.
   prob%fit_lij = .true.

   prob%experimental_points = exp_points

   X = [0.01, 0.01, 1., 1.]
   prob%parameter_step = [0.01, 0.01, 0.1, -0.1]
   print *, "X0:", X
   error = optimize(X, prob)
   print *, "FO:", error
   print *, "Xf:", X

   if (allocated(model)) deallocate(model)
   model = prob%get_model_from_X(X)

   ! ==========================================================================
   ! Write out results and experimental values
   ! --------------------------------------------------------------------------
   told = exp_points(1)%T
   do i=1, size(exp_points)
      point = saturation_pressure(&
         model, exp_points(i)%x, exp_points(i)%t, kind="bubble", &
         p0=exp_points(i)%p, y0=exp_points(i)%y &
         )
      if (told /= point%t) write(*, "(/)")
      print *, exp_points(i)%x(1), exp_points(i)%y(1), exp_points(i)%P, &
         point%x(1), point%y(1), point%P
      told = point%T
   end do
end program main
