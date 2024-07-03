module yaeos__consistency_gemodel
   !! # yaeos__consistency_gemodel
   !! Consistency checks of Helmholtz free energy models ([[GeModel]]).
   !!
   !! # Description
   !! This module contains tools to validate the analityc derivatives of
   !! implmented excess Gibbs free energy models ([[GeModel]]). Also, allows to
   !! evaluate the consistency tests described in Thermodynamic Models:
   !! Fundamentals & Computational Aspects 2 ed. by Michelsen and Mollerup
   !! Chapter 5 section 4.
   !!
   !! Available tools:
   !!
   !! - [[numeric_ge_derivatives]]: From an instantiated [[GeModel]] evaluate
   !! all the excess Gibbs free energy derivatives from the central finite
   !! difference method.
   !!
   !! - [[ge_consistency]]: From an instantiated GeModel evaluate all the
   !! Michelsen and Mollerup consistency tests
   !!
   !! # References
   !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
   !! Fundamentals & computational aspects (2. ed). Tie-Line Publications.
   !!
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel

   implicit none
contains
   subroutine ge_consistency(model, n, t, eq58, eq59, eq60, eq61)
      !! # ge_consistency
      !! \(G^E\) models consistency tests
      !!
      !! # Description
      !! Evaluate the \(G^E\) models consistency tests described in
      !! Thermodynamic Models: Fundamentals & Computational Aspects 2 ed. by
      !! Michelsen and Mollerup (MM) Chapter 5 section 4. The "eq" are
      !! evaluations of the left hand side of the following expressions:
      !!
      !! Equation 58
      !!
      !!  \[
      !!   \sum_i^{NC} n_i \text{ln} \gamma_i - \frac{G^E}{RT} = 0
      !!  \]
      !!
      !! Equation 59
      !!
      !!  \[
      !!   \text{ln} \gamma_i - \frac{1}{RT}
      !!   \frac{\partial G^E}{\partial n_i} = 0
      !!  \]
      !!
      !! Equation 60
      !!
      !!  \[
      !!   \frac{\partial \text{ln} \gamma_i}{\partial n_j} -
      !!   \frac{\partial \text{ln} \gamma_j}{\partial n_i} = 0
      !!  \]
      !!
      !! Equation 61
      !!
      !!  \[
      !!   \sum_i^{NC} n_i
      !!   \frac{\partial \text{ln} \gamma_i}{\partial n_j} = 0
      !!  \]
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos, only: pr
      !!  use yaeos, only: Groups, setup_unifac, UNIFAC
      !!  use yaeos__consistency_gemodel, only: ge_consistency
      !!
      !!  type(UNIFAC) :: model
      !!
      !!  integer, parameter :: nc = 4, ng = 4
      !!
      !!  type(Groups) :: molecules(nc)
      !!
      !!  real(pr) :: n(nc), T
      !!  real(pr) :: dt, dn
      !!
      !!  real(pr) :: eq58, eq59(nc), eq60(nc,nc), eq61(nc)
      !!
      !!  T = 303.15
      !!  n = [400.0, 100.0, 300.0, 200.0]
      !!
      !!  ! Hexane [CH3, CH2]
      !!  molecules(1)%groups_ids = [1, 2]
      !!  molecules(1)%number_of_groups = [2, 4]
      !!
      !!  ! Ethanol [CH3, CH2, OH]
      !!  molecules(2)%groups_ids = [1, 2, 14]
      !!  molecules(2)%number_of_groups = [1, 1, 1]
      !!
      !!  ! Toluene [ACH, ACCH3]
      !!  molecules(3)%groups_ids = [9, 11]
      !!  molecules(3)%number_of_groups = [5, 1]
      !!
      !!  ! Cyclohexane [CH2]
      !!  molecules(4)%groups_ids = [2]
      !!  molecules(4)%number_of_groups = [6]
      !!
      !!  model = setup_unifac(molecules)
      !!
      !!  ! ====================================================================
      !!  ! Consistency tests
      !!  ! --------------------------------------------------------------------
      !!  call ge_consistency(model, n, t, eq58, eq59, eq60, eq61)
      !! ```
      !!
      !! # References
      !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
      !! Fundamentals & computational aspects (2. ed). Tie-Line Publications.
      !!
      class(GeModel), intent(in) :: model
      !! \(G^E\) model
      real(pr), intent(in) :: n(:)
      !! Moles number vector
      real(pr), intent(in) :: t
      !! Temperature [K]
      real(pr), optional, intent(out) :: eq58
      !! MM Eq. 58
      real(pr), optional, intent(out) :: eq59(size(n))
      !! MM Eq. 59
      real(pr), optional, intent(out) :: eq60(size(n),size(n))
      !! MM Eq. 60
      real(pr), optional, intent(out) :: eq61(size(n))
      !! MM Eq. 61

      real(pr) :: Ge, Gen(size(n)), Gen2(size(n), size(n))
      real(pr) :: ln_gammas(size(n))

      integer i, j

      call model%excess_gibbs(n, t, Ge=Ge, Gen=Gen, Gen2=Gen2)
      call model%ln_activity_coefficient(n, t, ln_gammas)

      ! ========================================================================
      ! Equation 58
      ! ------------------------------------------------------------------------
      if (present(eq58)) then
         eq58 = sum(n * ln_gammas) - Ge / R / T
      end if

      ! ========================================================================
      ! Equation 59
      ! ------------------------------------------------------------------------
      if (present(eq59)) then
         eq59 = Gen / R / T - ln_gammas
      end if

      ! ========================================================================
      ! Equation 60
      ! ------------------------------------------------------------------------
      if (present(eq60)) then
         eq60 = 0.0_pr
         do i=1,size(n)
            do j=1,size(n)
               eq60(i,j) = Gen2(i,j) / R / T - Gen2(j,i) / R / T
            end do
         end do
      end if

      ! ========================================================================
      ! Equation 61
      ! ------------------------------------------------------------------------
      if (present(eq61)) then
         eq61 = 0.0_pr
         do j=1,size(n)
            eq61(j) = sum(n * Gen2(:, j) / R / T)
         end do
      end if
   end subroutine ge_consistency

   subroutine numeric_ge_derivatives(&
      model, n, t, d_n, d_t, Ge, GeT, Gen, GeT2, GeTn, Gen2 &
      )
      !! # numeric_ge_derivatives
      !! Numeric \(G^E\) model derivatives
      !!
      !! # Description
      !! Tool to facilitate the development of new [[GeModel]] by testing
      !! the implementation of analytic derivatives.
      !!
      !! # Examples
      !!
      !! ```fortran
      !! use yaeos, only: Groups, setup_unifac, UNIFAC
      !! use yaeos__consistency_gemodel, only: numeric_ge_derivatives
      !!
      !! type(UNIFAC) :: model
      !!
      !! integer, parameter :: nc = 4, ng = 4
      !!
      !! type(Groups) :: molecules(nc)
      !!
      !! real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      !! real(pr) :: Ge_n, Gen_n(nc), GeT_n, GeT2_n, GeTn_n(nc), Gen2_n(nc, nc)
      !! real(pr) :: ln_gammas(nc)
      !!
      !! real(pr) :: n(nc), T
      !! real(pr) :: dt, dn
      !!
      !! T = 303.15
      !! n = [400.0, 100.0, 300.0, 200.0] ! always test with sum(n) > 1
      !!
      !! dt = 0.1_pr
      !! dn = 0.1_pr
      !!
      !! ! Hexane [CH3, CH2]
      !! molecules(1)%groups_ids = [1, 2]
      !! molecules(1)%number_of_groups = [2, 4]
      !!
      !! ! Ethanol [CH3, CH2, OH]
      !! molecules(2)%groups_ids = [1, 2, 14]
      !! molecules(2)%number_of_groups = [1, 1, 1]
      !!
      !! ! Toluene [ACH, ACCH3]
      !! molecules(3)%groups_ids = [9, 11]
      !! molecules(3)%number_of_groups = [5, 1]
      !!
      !! ! Cyclohexane [CH2]
      !! molecules(4)%groups_ids = [2]
      !! molecules(4)%number_of_groups = [6]
      !!
      !! model = setup_unifac(molecules)
      !!
      !! ! =====================================================================
      !! ! Call analytic derivatives
      !! ! ---------------------------------------------------------------------
      !! call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !!
      !! ! =====================================================================
      !! ! Call numeric derivatives
      !! ! ---------------------------------------------------------------------
      !! call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeT=GeT_n)
      !! call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen=Gen_n)
      !! call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeT2=GeT2_n)
      !! call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeTn=GeTn_n)
      !! call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen2=Gen2_n)
      !! ```
      !!
      class(GeModel), intent(in) :: model
      !! \(G^E\) model
      real(pr), intent(in) :: n(:)
      !! Moles number vector
      real(pr), intent(in) :: t
      !! Temperature [K]
      real(pr), intent(in) :: d_n
      !! Moles finite difference step
      real(pr), intent(in) :: d_t
      !! Temperature finite difference step
      real(pr), intent(out) :: Ge
      !! Residual Helmoltz energy
      real(pr), optional, intent(out) :: GeT
      !! \(\frac{dGe}{dT}\)
      real(pr), optional, intent(out) :: Gen(size(n))
      !! \(\frac{dGe}{dn_i}\)
      real(pr), optional, intent(out) :: GeT2
      !! \(\frac{d^2Ge}{dT^2}\)
      real(pr), optional, intent(out) :: GeTn(size(n))
      !! \(\frac{d^2Ge}{dTdn_i}\)
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))
      !! \(\frac{d^2Ge}{dn_{ij}}\)

      ! Auxiliary
      real(pr) :: Ge_aux1, Ge_aux2, Ge_aux3, Ge_aux4
      real(pr) :: dn_aux1(size(n)), dn_aux2(size(n))
      integer :: i, j

      ! ========================================================================
      ! Ar valuations
      ! ------------------------------------------------------------------------
      ! on point valuation
      call model%excess_gibbs(n, t, Ge=Ge)

      ! ========================================================================
      ! Central numeric derivatives
      ! ------------------------------------------------------------------------
      ! Temperature
      if (present(GeT) .or. present(GeT2)) then
         call model%excess_gibbs(n, t + d_t, Ge=Ge_aux1)
         call model%excess_gibbs(n, t - d_t, Ge=Ge_aux2)

         if (present(GeT)) GeT = (Ge_aux1 - Ge_aux2) / (2 * d_t)
         if (present(GeT2)) GeT2 = (Ge_aux1 - 2 * Ge + Ge_aux2) / d_t**2
      end if

      ! Mole first derivatives
      if (present(Gen)) then
         Gen = 0.0_pr

         do i = 1, size(n), 1
            dn_aux1 = 0.0_pr
            dn_aux1(i) = d_n

            call model%excess_gibbs(n + dn_aux1, t, Ge=Ge_aux1)
            call model%excess_gibbs(n - dn_aux1, t, Ge=Ge_aux2)

            Gen(i) = (Ge_aux1 - Ge_aux2) / (2 * d_n)
         end do
      end if

      ! ========================================================================
      ! Central cross derivatives
      ! ------------------------------------------------------------------------
      ! Temperature - Mole
      if (present(GeTn)) then
         GeTn = 0.0_pr

         do i = 1, size(n), 1
            dn_aux1 = 0.0_pr
            dn_aux1(i) = d_n

            call model%excess_gibbs(n + dn_aux1, t + d_t, Ge=Ge_aux1)
            call model%excess_gibbs(n + dn_aux1, t - d_t, Ge=Ge_aux2)
            call model%excess_gibbs(n - dn_aux1, t + d_t, Ge=Ge_aux3)
            call model%excess_gibbs(n - dn_aux1, t - d_t, Ge=Ge_aux4)

            GeTn(i) = &
               (Ge_aux1 - Ge_aux2 - Ge_aux3 + Ge_aux4) / (4 * d_t * d_n)
         end do
      end if

      ! Mole second derivatives
      if (present(Gen2)) then
         Gen2 = 0.0_pr

         do i = 1, size(n), 1
            do j = 1, size(n), 1
               if (i .eq. j) then
                  dn_aux1 = 0.0_pr
                  dn_aux1(i) = d_n

                  call model%excess_gibbs(n + dn_aux1, t, Ge=Ge_aux1)
                  call model%excess_gibbs(n - dn_aux1, t, Ge=Ge_aux2)

                  Gen2(i, j) = (Ge_aux1 - 2 * Ge + Ge_aux2) / d_n**2
               else
                  dn_aux1 = 0.0_pr
                  dn_aux2 = 0.0_pr

                  dn_aux1(i) = d_n
                  dn_aux2(j) = d_n

                  call model%excess_gibbs(n + dn_aux1 + dn_aux2, t, Ge=Ge_aux1)
                  call model%excess_gibbs(n + dn_aux1 - dn_aux2, t, Ge=Ge_aux2)
                  call model%excess_gibbs(n - dn_aux1 + dn_aux2, t, Ge=Ge_aux3)
                  call model%excess_gibbs(n - dn_aux1 - dn_aux2, t, Ge=Ge_aux4)

                  Gen2(i, j) = &
                     (Ge_aux1 - Ge_aux2 - Ge_aux3 + Ge_aux4) / (4 * d_n**2)
               end if
            end do
         end do
      end if
   end subroutine numeric_ge_derivatives
end module yaeos__consistency_gemodel

