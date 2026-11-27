module auxiliar_functions
    !! # Auxiliary numerical functions
    !!
    !! # Description
    !! Relative-error and array-comparison helpers used in tests and
    !! consistency checks. The relative-error normalization is by the first
    !! argument, so callers should avoid a zero reference value.
    !!
    !! # Examples
    !! ```fortran
    !! if (allclose(reference, calculated, 1.0e-8_pr)) print *, "consistent"
    !! ```
    use yaeos__constants, only: pr
contains
    elemental function rel_error(x, y)
        !! # rel_error
        !! Return the absolute relative error with `x` as reference.
        !!
        !! # Description
        !! Computes \(|x-y|/|x|\). The reference `x` must be nonzero.
        real(pr), intent(in) :: x !! Nonzero reference value.
        real(pr), intent(in) :: y !! Value compared against the reference.
        real(pr) :: rel_error

        rel_error = abs(x - y)/abs(x)
    end function rel_error

    function allclose(x, y, rtol)
        !! Check whether all elements of two arrays agree within `rtol`.
        real(pr), intent(in) :: x(:) !! First array
        real(pr), intent(in) :: y(:) !! Second array
        real(pr), intent(in) :: rtol !! Relative tolerance

        logical :: allclose
        allclose = maxval(rel_error(x, y)) < rtol
        if (.not. allclose) then
            print *, "Max relative error:", maxval(rel_error(x, y))
            print *, "x:", x
            print *, "y:", y
        end if
    end function allclose
end module auxiliar_functions
