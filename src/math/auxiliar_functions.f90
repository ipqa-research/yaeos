module auxiliar_functions
    use yaeos__constants, only: pr
contains
    elemental function rel_error(x, y)
        real(pr), intent(in) :: x, y
        real(pr) :: rel_error

        rel_error = abs(x - y)/abs(x)
    end function rel_error

    function allclose(x, y, atol)
        real(pr), intent(in) :: x(:)
        real(pr), intent(in) :: y(:)
        real(pr), intent(in) :: atol

        logical :: allclose
        allclose = maxval(rel_error(x, y)) < atol
    end function allclose
end module auxiliar_functions
