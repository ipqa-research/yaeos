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
    
    real(pr) :: T, P, x1, y1, kij, X(4), told
    character(len=14) :: kind

    ! ==========================================================================
    ! Read file
    ! --------------------------------------------------------------------------
    allocate(exp_points(0))
    open(newunit=infile, file="fit_case", iostat=iostat)
    do
        read(infile, *, iostat=iostat) kind, t, p, x1, y1
        if (iostat /= 0) exit
        point = EquilibriaState(&
            kind=kind, T=T, P=P, x=[x1, 1-x1], y=[y1,1-y1], &
            Vx=0._pr, Vy=0._pr, iters=0, beta=0._pr &
        )
        exp_points = [exp_points, point]
    end do
    close(infile)
    ! ==========================================================================

    prob%experimental_points = exp_points
    prob%get_model_from_X => model_from_X

    forsus_dir = "build/dependencies/forsus/data/json"
    sus(1) = Substance("nitrogen")
    sus(2) = Substance("n-octane")

    X = 0.01
    print *, optimize(X, prob)

    call model_from_x(prob, X, model)

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

    subroutine model_from_X(problem, X, model)
        use yaeos, only: R, RKPR
        class(FittingProblem), intent(in) :: problem
        real(pr), intent(in) :: X(:)
        class(ArModel), allocatable, intent(out) :: model
        real(pr) :: kij(nc, nc), lij(nc, nc)

        real(pr) :: tc(nc), pc(nc), w(nc), vc(nc), zc(nc), del1(nc)

        tc = sus%critical%critical_temperature%value
        pc = sus%critical%critical_pressure%value/1e5
        w = sus%critical%acentric_factor%value
        vc = sus%critical%critical_volume%value

        zc = pc * vc / (R * tc)

        kij = 0
        kij(1, 2) = X(1)
        kij(2, 1) = kij(1, 2)

        lij = 0
        ! lij(1, 2) = X(2)
        ! lij(2, 1) = X(2)

        del1 = X(3:)

        model = RKPR(tc, pc, w, zc=zc, kij=kij, lij=lij)

        ! select type(model)
        !     type is (CubicEoS)
        !         model%del1 = del1
        !         model%del2 = (1._pr - model%del1)/(1._pr + model%del1)
        ! end select
    end subroutine
end program