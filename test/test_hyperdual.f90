program main
    use yaeos_autodiff
    use yaeos_constants, only: pr
    implicit none
    type(hyperdual) :: x, y, z

    x = 2.0_pr
    y = 3.0_pr

    x%f1 = 1
    y%f2 = 1

    z = f(x, y)

    if (abs(z%f0 - 2.4494897427831779) > 1e-5)   error stop 1
    if (abs(z%f1 - 0.61237243569579447) > 1e-5)  error stop 1
    if (abs(z%f2 - 0.40824829046386302) > 1e-5)  error stop 1
    if (abs(z%f12 - 0.10206207261596575) > 1e-5) error stop 1

contains

    function f(x, y)
        type(hyperdual) :: x, y, f
        f = sqrt(x*y)
    end function f
end program main
