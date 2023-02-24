module pengrobinson76
    use cubic_eos

    implicit none

    type, extends(CubicEOS) :: PR76
    end type

    interface setup_model
        module procedure :: setup_pr76
    end interface

    interface attractive_parameter
        module procedure :: a_pr
    end interface

    interface repulsive_parameter
        module procedure :: b_pr
    end interface
    
    interface del1_parameter
        module procedure :: del1_pr
    end interface
    
    interface del2_parameter
        module procedure :: del2_pr
    end interface

contains

    subroutine setup_pr76(model, n, tc, pc, w)
        type(PR76) :: model
        integer :: n
        real(pr) :: tc(n), pc(n), w(n)

        call alloc(model%CubicEoS, n)

        model%del1 = 1 + sqrt(2.0_pr)
        model%del2 = 1 - sqrt(2.0_pr)

        model%tc = tc
        model%pc = pc
        model%w = w

        model%ac = 0.45723553 * R**2 * tc**2 / pc
        model%b = 0.07779607*R*tc/pc
        model%k = 0.37464 + 1.54226*w - 0.26993*w**2
        model%c = 0
    end subroutine

    subroutine a_pr(model, z, v, t, a)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: a(size(model))

        call a_classic(model, z, v, t, a)
    end subroutine
    
    subroutine b_pr(model, z, p, v, t, b)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: b(size(model))

        call b_classic(model, z, v, t, b)
    end subroutine

    subroutine del1_pr(model, z, v, t, del1)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: p, v, t
        type(hyperdual) :: del1(size(model))

        del1 = 1 + sqrt(2.0_pr)
    end subroutine
    
    subroutine del2_pr(model, z, v, t, del2)
        type(PR76) :: model
        type(hyperdual) :: z(size(model))
        type(hyperdual) :: v, t
        type(hyperdual) :: del2(size(model))

        del2 = 1 - sqrt(2.0_pr)
    end subroutine

end module