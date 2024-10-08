module yaeos__models_ge_gc_td
   use yaeos__constants, only: pr
   use yaeos__models_ge_group_contribution_groups, only: Groups
   implicit none

   ! ===========================================================================
   ! PsiFunction that defines the temperature dependence of a UNIFAC-like model
   ! ---------------------------------------------------------------------------
   type, abstract :: PsiFunction
      !! # \(\psi(T)\) function
      !! UNIFAC \(\psi(T)\) functions abstract type
      !!
      !! # Description
      !! Abstract derived type for UNIFAC models temperature dependent functions
      !!
   contains
      procedure(temperature_dependence), deferred :: psi
   end type PsiFunction

   abstract interface
      subroutine temperature_dependence(&
         self, systems_groups, T, psi, dpsi_dt, dpsi_dt2&
         )
         !! # temperature_dependence interface
         !! Interface subroutine for UNIFAC models temperature dependent
         !! functions
         !!
         import pr, PsiFunction, Groups
         class(PsiFunction) :: self
         !! PsiFunction type variable
         class(Groups) :: systems_groups
         !! Groups type variable containig all the system's groups. See the
         !! `groups_stew` variable on the `UNIFAC` documentation.
         real(pr), intent(in) :: T
         !! Temperature [K]
         real(pr), optional, intent(out) :: psi(:, :)
         !! \(\psi(T)\)
         real(pr), optional, intent(out) :: dpsi_dt(:, :)
         !! \(\frac{d \psi (T)}{dT}\)
         real(pr), optional, intent(out) :: dpsi_dt2(:, :)
         !! \(\frac{d^2 \psi (T)}{dT^2}\)
      end subroutine temperature_dependence
   end interface

   ! ===========================================================================
   ! Implementations
   ! ---------------------------------------------------------------------------
   type, extends(PsiFunction) :: UNIFACPsi
      !! # Original UNIFAC \(\psi\) function
      !! \[
      !!    \psi_{ij}(T) = \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! \[
      !!    \frac{d \psi_{ij}(T)}{dT} = \frac{A_{ij}}{T^2}
      !!    \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! \[
      !!    \frac{d^2 \psi_{ij}(T)}{dT^2} =
      !!    \frac{Aij (Aij - 2T)}{T^4} \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.ddbst
      !! .com/published-parameters-unifac.html)
      !! 2. Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975).
      !! Group‐contribution estimation of activity coefficients in nonideal liquid
      !! mixtures. AIChE Journal, 21(6), 1086–1099.
      !! [https://doi.org/10.1002/aic.690210607](https://doi.org/10.1002/aic.690210607)
      !! 3. Skjold-Jorgensen, S., Kolbe, B., Gmehling, J., & Rasmussen, P. (1979).
      !! Vapor-Liquid Equilibria by UNIFAC Group Contribution. Revision and
      !! Extension. Industrial & Engineering Chemistry Process Design and
      !! Development, 18(4), 714–722.
      !! [https://doi.org/10.1021/i260072a024](https://doi.org/10.1021/i260072a024)
      !! 4. Gmehling, J., Rasmussen, P., & Fredenslund, A. (1982). Vapor-liquid
      !! equilibriums by UNIFAC group contribution. Revision and extension. 2.
      !! Industrial & Engineering Chemistry Process Design and Development, 21(1),
      !! 118–127.
      !! [https://doi.org/10.1021/i200016a021](https://doi.org/10.1021/i200016a021)
      !! 5. Macedo, E. A., Weidlich, U., Gmehling, J., & Rasmussen, P. (1983).
      !! Vapor-liquid equilibriums by UNIFAC group contribution. Revision and
      !! extension. 3. Industrial & Engineering Chemistry Process Design and
      !! Development, 22(4), 676–678.
      !! [https://doi.org/10.1021/i200023a023](https://doi.org/10.1021/i200023a023)
      !! 6. Tiegs, D., Rasmussen, P., Gmehling, J., & Fredenslund, A. (1987).
      !! Vapor-liquid equilibria by UNIFAC group contribution. 4. Revision and
      !! extension. Industrial & Engineering Chemistry Research, 26(1), 159–161.
      !! [https://doi.org/10.1021/ie00061a030](https://doi.org/10.1021/ie00061a030)
      !! 7. Hansen, H. K., Rasmussen, P., Fredenslund, A., Schiller, M., &
      !! Gmehling, J. (1991). Vapor-liquid equilibria by UNIFAC group
      !! contribution. 5. Revision and extension. Industrial & Engineering
      !! Chemistry Research, 30 (10), 2352–2355.
      !! [https://doi.org/10.1021/ie00058a017](https://doi.org/10.1021/ie00058a017)
      !! 8. Wittig, R., Lohmann, J., & Gmehling, J. (2003). Vapor−Liquid Equilibria
      !! by UNIFAC Group Contribution. 6. Revision and Extension. Industrial &
      !! Engineering Chemistry Research, 42(1), 183–188.
      !! [https://doi.org/10.1021/ie020506l](https://doi.org/10.1021/ie020506l)
      !! 9. [SINTEF - Thermopack](https://github.com/thermotools/thermopack)
      !!
      real(pr), allocatable :: Aij(:, :)
   contains
      procedure :: psi => UNIFAC_temperature_dependence
   end type UNIFACPsi

   type, extends(PsiFunction) :: QuadraticPsi
      real(pr), allocatable :: Aij(:, :)
      real(pr), allocatable :: Bij(:, :)
      real(pr), allocatable :: Cij(:, :)
   contains
      procedure :: psi => Quadratic_temperature_dependence
   end type QuadraticPsi

contains

   subroutine UNIFAC_temperature_dependence(&
      self, systems_groups, T, psi, dpsi_dt, dpsi_dt2 &
      )
      !! # UNIFAC temperature dependence
      !! Implementation of the \(\psi(T) \) function of the UNIFAC model.
      !!
      !! \[
      !!    \psi_{ij}(T) = \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! \[
      !!    \frac{d \psi_{ij}(T)}{dT} = \frac{A_{ij}}{T^2}
      !!    \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! \[
      !!    \frac{d^2 \psi_{ij}(T)}{dT^2} =
      !!    \frac{Aij (Aij - 2T)}{T^4} \exp(-\frac{A_{ij}}{T})
      !! \]
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.ddbst
      !! .com/published-parameters-unifac.html)
      !! 2. Fredenslund, A., Jones, R. L., & Prausnitz, J. M. (1975).
      !! Group‐contribution estimation of activity coefficients in nonideal liquid
      !! mixtures. AIChE Journal, 21(6), 1086–1099.
      !! [https://doi.org/10.1002/aic.690210607](https://doi.org/10.1002/aic.690210607)
      !! 3. Skjold-Jorgensen, S., Kolbe, B., Gmehling, J., & Rasmussen, P. (1979).
      !! Vapor-Liquid Equilibria by UNIFAC Group Contribution. Revision and
      !! Extension. Industrial & Engineering Chemistry Process Design and
      !! Development, 18(4), 714–722.
      !! [https://doi.org/10.1021/i260072a024](https://doi.org/10.1021/i260072a024)
      !! 4. Gmehling, J., Rasmussen, P., & Fredenslund, A. (1982). Vapor-liquid
      !! equilibriums by UNIFAC group contribution. Revision and extension. 2.
      !! Industrial & Engineering Chemistry Process Design and Development, 21(1),
      !! 118–127.
      !! [https://doi.org/10.1021/i200016a021](https://doi.org/10.1021/i200016a021)
      !! 5. Macedo, E. A., Weidlich, U., Gmehling, J., & Rasmussen, P. (1983).
      !! Vapor-liquid equilibriums by UNIFAC group contribution. Revision and
      !! extension. 3. Industrial & Engineering Chemistry Process Design and
      !! Development, 22(4), 676–678.
      !! [https://doi.org/10.1021/i200023a023](https://doi.org/10.1021/i200023a023)
      !! 6. Tiegs, D., Rasmussen, P., Gmehling, J., & Fredenslund, A. (1987).
      !! Vapor-liquid equilibria by UNIFAC group contribution. 4. Revision and
      !! extension. Industrial & Engineering Chemistry Research, 26(1), 159–161.
      !! [https://doi.org/10.1021/ie00061a030](https://doi.org/10.1021/ie00061a030)
      !! 7. Hansen, H. K., Rasmussen, P., Fredenslund, A., Schiller, M., &
      !! Gmehling, J. (1991). Vapor-liquid equilibria by UNIFAC group
      !! contribution. 5. Revision and extension. Industrial & Engineering
      !! Chemistry Research, 30 (10), 2352–2355.
      !! [https://doi.org/10.1021/ie00058a017](https://doi.org/10.1021/ie00058a017)
      !! 8. Wittig, R., Lohmann, J., & Gmehling, J. (2003). Vapor−Liquid Equilibria
      !! by UNIFAC Group Contribution. 6. Revision and Extension. Industrial &
      !! Engineering Chemistry Research, 42(1), 183–188.
      !! [https://doi.org/10.1021/ie020506l](https://doi.org/10.1021/ie020506l)
      !! 9. [SINTEF - Thermopack](https://github.com/thermotools/thermopack)
      !!
      class(UNIFACPsi) :: self
      !! \(\psi\) function
      class(Groups) :: systems_groups
      !! Groups in the system
      real(pr), intent(in) :: T
      !! Temperature [K]
      real(pr), optional, intent(out) :: psi(:, :)
      !! \(\psi\)
      real(pr), optional, intent(out) :: dpsi_dt(:, :)
      !! \(\frac{d \psi}{dT}\)
      real(pr), optional, intent(out) :: dpsi_dt2(:, :)
      !! \(\frac{d^2 \psi}{dT^2}\)

      integer :: i, j
      integer :: ngroups

      real(pr) :: Aij
      real(pr) :: Eij

      ngroups = size(systems_groups%groups_ids)

      do concurrent(i=1:ngroups, j=1:ngroups)
         Aij = self%Aij(i, j)
         Eij = exp(-Aij / T)
         if (present(psi)) &
            psi(i, j) = Eij
         if (present(dpsi_dt)) &
            dpsi_dt(i, j) = Aij * Eij / T**2
         if (present(dpsi_dt2)) &
            dpsi_dt2(i, j) = Aij * (Aij - 2_pr*T) * Eij / T**4
      end do
   end subroutine UNIFAC_temperature_dependence

   subroutine Quadratic_temperature_dependence(&
      self, systems_groups, T, psi, dpsi_dt, dpsi_dt2 &
      )
      !! # Quadratic temperature dependence
      class(QuadraticPsi) :: self
      !! \(\psi\) function
      class(Groups) :: systems_groups
      !! Groups in the system
      real(pr), intent(in) :: T
      !! Temperature [K]
      real(pr), optional, intent(out) :: psi(:, :)
      !! \(\psi\)
      real(pr), optional, intent(out) :: dpsi_dt(:, :)
      !! \(\frac{d \psi}{dT}\)
      real(pr), optional, intent(out) :: dpsi_dt2(:, :)
      !! \(\frac{d^2 \psi}{dT^2}\)

      integer :: i, j
      integer :: ngroups

      real(pr) :: u, dudt, dudt2
      real(pr) :: a, b, c

      ngroups = size(systems_groups%groups_ids)

      do concurrent(i=1:ngroups, j=1:ngroups)
         a = self%Aij(i, j)
         b = self%Bij(i, j)
         c = self%Cij(i, j)

         u = -(A + B*T + C*T**2)/T
         dudt = a / T**2 - c
         dudt2 = -2._pr * a / T**3

         if (present(psi)) then
            psi(i, j) = exp(u)
         end if

         if (present(dpsi_dt)) then
            dpsi_dt(i, j) = dudt * exp(u)
         end if

         if (present(dpsi_dt2)) then
            dpsi_dt2(i, j) = (dudt2 + dudt**2)*exp(u)
         end if

      end do
   end subroutine Quadratic_temperature_dependence
end module yaeos__models_ge_gc_td
