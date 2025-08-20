module auxiliar_functions
    !! # Test Auxiliary Functions
    !!
    !! This module provides utility functions for testing and validation
    !! purposes. These functions are primarily used in the test suite
    !! to compare numerical results and assess convergence.
    !!
    !! ## Main Features
    !!
    !! - **Relative error calculation**: `rel_error` function
    !! - **Array comparison**: `allclose` function for comparing arrays
    !! - **Numerical tolerance testing**: Built-in tolerance checking
    !!
    !! ## Usage Examples
    !!
    !! ```fortran
    !! use auxiliar_functions
    !!
    !! real(pr) :: expected(3) = [1.0, 2.0, 3.0]
    !! real(pr) :: calculated(3) = [1.001, 1.998, 3.002]
    !! logical :: test_passed
    !!
    !! ! Check if arrays are close within tolerance
    !! test_passed = allclose(calculated, expected, 1e-2_pr)
    !! ```
    !!
    !! These functions are essential for ensuring numerical accuracy
    !! in yaeos calculations and for regression testing.
    use yaeos__constants, only: pr
contains
    elemental function rel_error(x, y)
        real(pr), intent(in) :: x, y
        real(pr) :: rel_error

        rel_error = abs(x - y)/abs(x)
    end function rel_error

    function allclose(x, y, rtol)
        real(pr), intent(in) :: x(:)
        real(pr), intent(in) :: y(:)
        real(pr), intent(in) :: rtol

        logical :: allclose
        allclose = maxval(rel_error(x, y)) < rtol
        if (.not. allclose) then
            print *, "Max relative error:", maxval(rel_error(x, y))
            print *, "x:", x
            print *, "y:", y
        end if
    end function allclose
end module auxiliar_functions
