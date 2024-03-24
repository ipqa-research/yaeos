module yaeos_c
   !! C interface to yaeos subroutines
   use iso_c_binding, only: c_double, c_int
   use yaeos, only: &
                    & ArModel, & ! Ar Model
                    ! Cubic Models
                    & SoaveRedlichKwong, PengRobinson76, PengRobinson78, &
                    ! Thermodynamic properties
                    & fugacity_vt
   implicit none

   type :: ModelContainer
      integer :: id
      class(ArModel), pointer :: model
   end type

   class(ArModel), allocatable, private :: model
   class(ModelContainer), allocatable :: models(:)

   integer(c_int) :: running_model=0

contains

   subroutine pr76(tc, pc, w, kij, lij)
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      model = PengRobinson76(tc, pc, w, kij, lij)
   end subroutine
   
   subroutine srk(tc, pc, w, kij, lij)
      real(c_double), intent(in) :: tc(:), pc(:), w(:), kij(:, :), lij(:, :)
      model = SoaveRedlichKwong(tc, pc, w, kij, lij)
   end subroutine

   subroutine fug_vt(n, v, t, lnfug, dlnphidp, dlnphidt, dlnphidn)
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
