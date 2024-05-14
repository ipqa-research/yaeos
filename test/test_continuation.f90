program test_continuation
    use yaeos_constants, only: pr
    use yaeos__math_continuation, only: continuation
    use fixtures_models, only: binary_PR76

    integer, parameter :: max_points=1000
    real(pr) :: XS(max_points, 2)
    integer :: i

    XS = continuation(&
        circle, &
        X0=[0.1_pr, 0.9_pr], ns0=2, S0=0.9_pr, dS0=0.01_pr, &
        max_points=max_points, solver_tol=1e-5_pr)
contains
    subroutine circle(X, ns, S, F, dF, dFdS)
        real(pr), intent(in) :: X(:)
        integer, intent(in) :: ns
        real(pr), intent(in) :: S
        real(pr), intent(out) :: F(:), dF(:, :), dFdS(:)

        f = 0
        df = 0

        F(1) = X(1)**2 + X(2)**2 - 1
        F(2) = X(ns) - S

        df(1, :)  = [2*X(1), 2*X(2)]
        df(2, ns) = 1

        dFdS = 0
        dFdS(2) = -1
    end subroutine
end program