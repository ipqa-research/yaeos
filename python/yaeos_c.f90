module yaeos_c
   !! C interface to yaeos subroutines
   use iso_c_binding, only: c_double, c_int
   use yaeos, only: ArModel, PengRobinson76, fugacity_vt, CubicEoS
   implicit none

   class(ArModel), allocatable, private :: model

   integer(c_int) :: running_model=0

contains

   subroutine pr76(tc, pc, w, kij, lij)
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      model = PengRobinson76(tc, pc, w, kij, lij)
   end subroutine

   subroutine fugacity(n, v, t, lnfug, dlnphidp, dlnphidt, dlnphidn)
      real(c_double), intent(in) :: n(:), v, t
      real(c_double), intent(out) :: lnfug(size(n))
      real(c_double) :: p

      real(c_double), optional, intent(in out) :: dlnphidp(size(n)), dlnphidt(size(n)), dlnphidn(size(n), size(n))

      call fugacity_vt(&
         model, &
         n, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
      )
   end subroutine
end module
