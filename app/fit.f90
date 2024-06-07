program main
    !! Binary system parameter optimization
    use yaeos, only: EquilibriaState, pr, ArModel, PengRobinson78, CubicEoS, saturation_pressure
    use forsus, only: Substance, forsus_dir
    use yaeos__fitting, only: FittingProblem, fobj, optimize
    integer, parameter :: nc=2
    integer :: i, infile, iostat

    type(EquilibriaState), allocatable :: exp_points(:)
    type(EquilibriaState) :: point
    type(FittingProblem) :: prob
    type(Substance) :: sus(2)

    class(ArModel), allocatable :: model

    real(pr) :: T, P, x1, y1, kij, X(4), told, error
    character(len=14) :: kind

    ! ==========================================================================
    ! Setup componentes and read data file
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
    prob%experimental_points = exp_points
    prob%get_model_from_X => model_from_X

    X = [0.01, 0.01, 1., 1.]
    prob%parameter_step = [0.01, 0.01, 0.1, -0.1]
    print *, "X0:", X
    error = optimize(X, prob)
    print *, "FO:", error
    print *, "Xf:", X

    model = prob%get_model_from_X(X, prob)

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

contains

    function model_from_X(X, problem) result(model)
        use yaeos, only: R, RKPR, PengRobinson76
        real(pr), intent(in) :: X(:)
        class(FittingProblem), intent(in) :: problem
        class(ArModel), allocatable :: model
        real(pr) :: kij(nc, nc), lij(nc, nc)

        real(pr) :: tc(nc), pc(nc), w(nc), vc(nc), zc(nc), mult(nc)

        tc = sus%critical%critical_temperature%value
        pc = sus%critical%critical_pressure%value/1e5
        w = sus%critical%acentric_factor%value
        vc = sus%critical%critical_volume%value

        zc = pc * vc / (R * tc)

        kij = 0
        kij(1, 2) = X(1)
        kij(2, 1) = kij(1, 2)

        lij = 0
        lij(1, 2) = X(2)
        lij(2, 1) = X(2)

        model = RKPR(tc, pc, w, zc=zc, kij=kij, lij=lij)
        model = PengRobinson78(tc, pc, w, kij=kij, lij=lij)

        mult = X(3:)
        ! select type(model)
        !     type is (CubicEoS)
        !         model%del1 = mult * model%del1
        !         model%del2 = (1._pr - model%del1)/(1._pr + model%del1)
        ! end select
    end function
end program