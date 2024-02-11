module yaeos_tapenade_ar_api
   !! Module that wraps tapenade generated routines to calculate !
   !! Ar and derivatives.
   use yaeos_constants, only: pr
   use yaeos_models_ar, only: ArModel
   implicit none

   private

   public :: ArModelTapenade

   type, extends(ArModel) :: ArModelTapenade
      procedure(tapenade_ar), pointer, nopass :: ar
      procedure(tapenade_ar_d), pointer, nopass :: ar_d
      procedure(tapenade_ar_b), pointer, nopass :: ar_b
      procedure(tapenade_ar_d_b), pointer, nopass :: ar_d_b
      procedure(tapenade_ar_d_d), pointer, nopass :: ar_d_d
      procedure(tapenade_v0), pointer, nopass :: v0
   contains
      procedure :: residual_helmholtz => residual_helmholtz
      procedure :: get_v0 => get_v0
   end type

   abstract interface
      subroutine tapenade_ar(n, v, t, arval)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(out) :: arval
      end subroutine

      subroutine tapenade_ar_d(n, nd, v, vd, t, td, arval, arvald)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald
      end subroutine

      subroutine tapenade_ar_b(n, nb, v, vb, t, tb, arval, arvalb)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arvalb
         real(pr) :: nb(:), vb, tb
         real(pr) :: arval
      end subroutine

      subroutine tapenade_ar_d_b(n, nb, v, vb, t, tb, arval, arvalb)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arvalb
         real(pr) :: nb(:), vb, tb
         real(pr) :: arval
      end subroutine

      subroutine tapenade_ar_d_d(n, nd, v, vd0, vd, t, td0, td, &
            arval, arvald0, arvald, arvaldd)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: vd0, td0
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald0, arvald, arvaldd
      end subroutine

      function tapenade_v0(n, p, t)
         import pr
         real(pr), intent(in) :: n(:), p, t
         real(pr) :: tapenade_v0
      end function
   end interface

contains
   subroutine residual_helmholtz(&
         self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz model generic interface
      class(ArModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: v, t
      real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
      real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))

      if(present(Arn2)) then
      end if
   end subroutine

   function get_v0(self, n, p, t)
      class(ArModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:), p, t
      real(pr) :: get_v0

      get_v0 = self%v0(n, p, t)
   end function
end module