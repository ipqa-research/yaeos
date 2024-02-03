module yaeos_tapenade_ar_api
   !! Module that wraps tapenade generated routines to calculate !
   !! Ar and derivatives.
   use yaeos_constants, only: pr
   implicit none

   private

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
   end interface

contains
   subroutine residual_helmholtz( &
         n, v, t, &
         ar, dardn, dardv, dardt, &
         dardn2, darnv, dardnt, dardvt, dardv2, dardt2 &
      )
      real(pr), intent(in) :: n(:), v, t
      real(pr), optional, intent(out) :: ar, dardn(size(n)), dardv, dardt
      real(pr), optional, intent(out) :: dardn2(size(n), size(n)), darnv(size(n)), dardnt(size(n))
      real(pr), optional, intent(out) :: dardvt, dardv2, dardt2
   end subroutine
end module