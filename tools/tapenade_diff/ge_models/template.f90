module tapenade_ge_model_template
    use yaeos__constants, only: pr, R
    use yaeos__tapenade_ge_api, only: GeModelTapenade
    implicit none

   
    <parameters here>

    type(GeModelTapenade) :: model

contains

    subroutine setup(<parameters>)
        <parameters definition>

        model%ge => excess_gibbs
        model%ge_b => excess_gibbs_b
        model%ge_d => excess_gibbs_d
        model%ge_d_b => excess_gibbs_d_b
        model%ge_d_d => excess_gibbs_d_b
    end subroutine

    subroutine excess_gibbs(n, T, ge)
        real(8), intent(in) :: n(:)
        real(8), intent(in) :: T
        real(8), intent(out) :: ge

        <your Ge Model implementation here>
    end subroutine

end module