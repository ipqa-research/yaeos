module yaeos__adiff_hyperdual_ar_api
    !! Module that contains the automatic differentiation logic for an Ar model.
    !! 
    !! All that is needed to define an Ar model that uses automatic
    !! differentiation with hyperdual numbers is to define a new derived type
    !! that overloads the method to the Ar function that you want to use.
    !! A minimal example follows:
    !!
    !! ```fortran
    !! module newmodel
    !! use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
    !!
    !! type, extends(ArModelAdiff) :: YourNewModel
    !!       type(Substances) :: composition
    !!       real(8) :: parameters(:)
    !!     contains
    !!       procedure :: Ar => arfun
    !!       procedure :: get_v0 => v0
    !! end type
    !! contains
    !! subroutine arfun(self, n, v, t, Ar)
    !!    class(YourNewModel), intent(in) :: self
    !!    type(hyperdual), intent(in) :: n(:) ! Number of moles
    !!    type(hyperdual), intent(in) :: v ! Volume [L]
    !!    type(hyperdual), intent(in) :: t ! Temperature [K]
    !!    type(hyperdual), intent(out) :: ar_value ! Residual Helmholtz Energy
    !! 
    !!    ! A very complicated residual helmholtz function of a mixture
    !!    Ar = sum(n) * v * t
    !! end subroutine
    !! 
    !! function v0(self, n, p, t)
    !!    class(YourNewModel), intent(in) :: self
    !!    real(pr), intent(in) :: n(:) ! Number of moles
    !!    real(pr), intent(in) :: p ! Pressure [bar]
    !!    real(pr), intent(in) :: t ! Temperature [K]
    !!    real(pr) :: v0
    !! 
    !!    v0 = self%parameters(3)
    !! end function
    !! ```
    !! 
    !! A complete implementation of the PR76 Equation of State can me found in
    !! `example/adiff/adiff_pr76.f90`
    !!
    use yaeos__constants, only: pr
    use yaeos__models_ar, only: ArModel
    use hyperdual_mod

    implicit none

    type, abstract, extends(ArModel) :: ArModelAdiff
    contains
        procedure(hyperdual_ar), deferred :: Ar
        procedure :: residual_helmholtz => residual_helmholtz
    end type

    abstract interface
        type(hyperdual) function hyperdual_Ar(self, n, v, t)
            import hyperdual, ArModelAdiff
            class(ArModelAdiff) :: self
            type(hyperdual), intent(in) :: n(:)
            type(hyperdual), intent(in) :: v
            type(hyperdual), intent(in) :: t
        end function
    end interface
contains

    subroutine residual_helmholtz(&
        self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
    )
        class(ArModelAdiff), intent(in) :: self
        real(pr), intent(in) :: n(:)
        real(pr), intent(in) :: v, t
        real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
        real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
        real(pr), optional, intent(out) :: Arn2(size(n), size(n))

        type(hyperdual) :: d_v, d_t, d_n(size(n))
        type(hyperdual) :: d_Ar

        if (present(ArV)) then
            if (present(ArV2)) call get_dardv2
            if (present(ArVn)) call get_dardvn
            if (present(ArTV)) call get_dardvt
            if (.not. (present(ArV2) .and. present(ArVn) .and. present(ArTV))) &
                call get_dardv
        end if

        if (present(ArT)) then
            if (present(ArT2)) call get_dardt2
            if (present(ArTn)) call get_dardtn
            if (.not. (present(ArT2) .and. present(ArTn))) call get_dardt
        end if

        if (present(Arn)) then
            if (present(Arn2)) then
                call get_dardn2
            else
                call get_dardn
            end if
        end if

        if (present(Ar)) Ar = d_Ar%f0

    contains

        subroutine get_dardn()
            integer :: i, j

            do i=1,size(n)
                call reset_vars
                d_n(i)%f1 = 1
                d_Ar = self%Ar(d_n, d_v, d_t)
                Arn(i) = d_Ar%f1
            end do
        end subroutine
        
        subroutine get_dardn2()
            integer :: i, j

            do i=1,size(n)
                do j=i,size(n)
                    call reset_vars
                    d_n(i)%f1 = 1
                    d_n(j)%f2 = 1
                    
                    d_Ar = self%Ar(d_n, d_v, d_t)

                    Arn(i) = d_Ar%f1
                    Arn2(i, j) = d_Ar%f12
                    Arn2(j, i) = d_Ar%f12
                end do
            end do
        end subroutine
       
        subroutine get_dardvn()
            integer :: i

            do i=1,size(n)
                call reset_vars
                d_n(i)%f1 = 1
                d_v%f2 = 1
                d_Ar = self%Ar(d_n, d_v, d_t)
                Arn(i) = d_Ar%f1
                ArV = d_Ar%f2
                ArVn(i) = d_Ar%f12
            end do
        end subroutine
        
        subroutine get_dardtn()
            integer :: i

            do i=1,size(n)
                call reset_vars
                d_n(i)%f1 = 1
                d_t%f2 = 1
                d_Ar = self%Ar(d_n, d_v, d_t)
                Arn(i) = d_Ar%f1
                ArT = d_Ar%f2
                ArTn(i) = d_Ar%f12
            end do
        end subroutine

        subroutine get_dardv()
            call reset_vars
            d_v%f1 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            ArV = d_Ar%f1
        end subroutine
        
        subroutine get_dardt()
            call reset_vars
            d_t%f1 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            ArT = d_Ar%f1
        end subroutine
        
        subroutine get_dardt2()
            call reset_vars
            d_t%f1 = 1
            d_t%f2 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            ArT = d_Ar%f1
            ArT2 = d_Ar%f12
        end subroutine
        
        subroutine get_dardv2()
            call reset_vars
            d_v%f1 = 1
            d_v%f2 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            ArV = d_Ar%f1
            ArV2 = d_Ar%f12
        end subroutine

        subroutine get_dardvt()
            call reset_vars
            d_v%f1 = 1
            d_t%f2 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            ArV = d_Ar%f1
            ArTV = d_Ar%f12
        end subroutine

        subroutine reset_vars()
            d_n = n
            d_v = v
            d_t = t
        end subroutine
    end subroutine

end module