#:include 'cubic/general.fypp'
#:include 'cubic/alphas.fypp'
#:include 'cubic/mixing_rules.fypp'
#:include 'cubic/ar_funs.fypp'
! =============================================================================
! -> Setup or select your model here

#:set alpha=soave
#:set mixingrule=qmr
#:set ar=generic_cubic
#:set modelname="PR"

! =============================================================================
module model
    use yaeos___constants, only: pr
    use yaeos___models_ar, only: ArModel
    implicit none

    type :: ${modelname}$
        @:critical()
        @:ar(params=True)
        @:alpha(params=True)
        @:mixingrule(params=True)
    end type

contains

    subroutine ar(self, n, v, t, arval)
        class(${modelname}$) :: self
        real(8), intent(in) :: n(:)
        real(8), intent(in) :: v, t
        real(8), intent(out) :: arval
        @:ar(vars=True)
        @:alpha(vars=True)
        @:mixingrule(vars=True)
        
        @:alpha(eqs=True)
        @:mixingrule(eqs=True)
        @:ar(eqs=True)
    end subroutine

    pure function volume_initalizer(n, p, t) result(v0)
        real(8), intent(in) :: n(:)
        real(8), intent(in) :: p
        real(8), intent(in) :: t
        real(8) :: v0
        v0 = sum(n*b)/sum(b)
    end function
end module