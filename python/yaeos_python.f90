module yaeos_python
   use yaeos, only: ArModel, PengRobinson76, fugacity_vt, CubicEoS
   implicit none

   class(ArModel), allocatable, private :: model

   integer(8) :: running_model=0

contains

   subroutine pr76(tc, pc, w, kij, lij)
      real(8), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      model = PengRobinson76(tc, pc, w, kij, lij)
   end subroutine

   subroutine fugacity(n, v, t, lnfug, dlnphidp, dlnphidt, dlnphidn)
      real(8), intent(in) :: n(:), v, t
      real(8), intent(out) :: lnfug(size(n))
      real(8) :: p

      real(8), optional, intent(in out) :: dlnphidp(size(n)), dlnphidt(size(n)), dlnphidn(size(n), size(n))

      call fugacity_vt(&
         model, &
         n, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
      )
   end subroutine
end module
