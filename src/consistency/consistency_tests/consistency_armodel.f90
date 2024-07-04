module yaeos__consistency_armodel
   !! # yaeos__consistency_armodel
   !! Consistency checks of Helmholtz free energy models ([[ArModel]]).
   !!
   !! # Description
   !! This module contains tools to validate the analityc derivatives of
   !! implmented Helmholtz free energy models ([[ArModel]]). Also, allows to
   !! evaluate the consistency tests described in Thermodynamic Models:
   !! Fundamentals & Computational Aspects 2 ed. by Michelsen and Mollerup
   !! Chapter 2 section 3.
   !!
   !! Available tools:
   !!
   !! - [[numeric_ar_derivatives]]: From an instantiated [[ArModel]] evaluate
   !! all the Helmholtz free energy derivatives from the central finite
   !! difference method.
   !!
   !! - [[ar_consistency]]: From an instantiated [[ArModel]] evaluate all the
   !! Michelsen and Mollerup consistency tests.
   !!
   !! # References
   !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
   !! Fundamentals & computational aspects (2. ed). Tie-Line Publications.
   !!
   use yaeos__constants, only: pr, R
   use yaeos__models_ar, only: ArModel
   use yaeos__thermoprops, only: enthalpy_residual_vt, gibbs_residual_vt
   use yaeos__thermoprops, only: fugacity_vt, pressure

   implicit none
contains
   subroutine ar_consistency(&
      eos, n, v, t, eq31, eq33, eq34, eq36, eq37 &
      )
      !! # ar_consistency
      !! \(A^r\) models consistency tests.
      !!
      !! # Description
      !! The evaluated equations are taken from Fundamentals & Computational
      !! Aspects 2 ed. by Michelsen and Mollerup Chapter 2 section 3. The
      !! "eq" are evaluations of the left hand side of the following
      !! expressions:
      !!
      !! Equation 31:
      !!
      !! \[\sum_i n_i ln \hat{\phi}_i - \frac{G^r(T,P,n)}{RT} = 0\]
      !!
      !! Equation 33:
      !!
      !! \[
      !!    \left(\frac{\partial ln \hat{\phi}_i}{\partial n_j} \right)_{T,P}
      !!    - \left(\frac{\partial ln \hat{\phi}_j}{\partial n_i} \right)_{T,P}
      !!    = 0
      !! \]
      !!
      !! Equation 34:
      !!
      !! \[
      !!    \sum_i n_i
      !!    \left(\frac{\partial ln \hat{\phi}_i}{\partial n_j} \right)_{T,P}
      !!    = 0
      !! \]
      !!
      !! Equation 36:
      !!
      !! \[
      !!    \left(\frac{\partial}{\partial P}
      !!    \sum_i n_i ln \hat{\phi}_i \right)_{T,n} - \frac{(Z - 1)n}{P} = 0
      !! \]
      !!
      !! Equation 37:
      !!
      !! \[
      !!    \sum_i n_i \left(\frac{\partial ln \hat{\phi}_i}{\partial T}
      !!    \right)_{P,n} + \frac{H^r(T,P,n)}{RT^2} = 0
      !! \]
      !!
      !! The consistency test could be applied to any instantiated [[ArModel]]
      !! as shown in the following example.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos, only: pr, SoaveRedlichKwong, ArModel
      !!  use yaeos__consistency_armodel, only: ar_consistency
      !!
      !!  class(ArModel), allocatable :: model
      !!  real(pr) :: tc(4), pc(4), w(4)
      !!
      !!  real(pr) :: n(4), t, v
      !!
      !!  real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37
      !!
      !!  n = [1.5, 0.2, 0.7, 2.3]
      !!  tc = [190.564, 425.12, 300.11, 320.25]
      !!  pc = [45.99, 37.96, 39.23, 40.21]
      !!  w = [0.0115478, 0.200164, 0.3624, 0.298]
      !!
      !!  t = 600_pr
      !!  v = 0.5_pr
      !!
      !!  model = SoaveRedlichKwong(tc, pc, w)
      !!
      !!  call ar_consistency(&
      !!     model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
      !!     )
      !! ```
      !! All `eqXX` variables should be close to zero.
      !!
      !! # References
      !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
      !! Fundamentals & computational aspects (2. ed). Tie-Line Publications.
      !!
      class(ArModel), intent(in) :: eos !! Equation of state
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), optional, intent(out) :: eq31 !! MM Eq. 31
      ! TODO real(pr), optional, intent(out) :: eq32
      real(pr), optional, intent(out) :: eq33(size(n), size(n)) !! MM Eq. 33
      real(pr), optional, intent(out) :: eq34(size(n)) !! MM Eq. 34
      real(pr), optional, intent(out) :: eq36 !! MM Eq. 36
      real(pr), optional, intent(out) :: eq37 !! MM Eq. 37

      integer i, j

      ! ========================================================================
      ! Previous calculations
      ! ------------------------------------------------------------------------
      real(pr) :: Grp, Grv, Hrv, p, dpdn(size(n)), ntot, z
      real(pr) :: lnphi(size(n)), lnphip(size(n)), dlnPhidP(size(n))
      real(pr) :: dlnPhidT(size(n)), dlnPhidn(size(n), size(n))

      call pressure(eos, n, v, t, p, dpdn=dpdn)

      call gibbs_residual_vt(eos, n, v, t, Grv)

      call enthalpy_residual_vt(eos, n, v, t, Hr=Hrv)

      call fugacity_vt(&
         eos, n, v, t, lnphip=lnphip, &
         dlnPhidP=dlnPhidP, dlnPhidT=dlnPhidT, dlnPhidn=dlnPhidn &
         )

      ntot = sum(n)

      lnphi(:) = lnphip(:) - log(p)

      z = p * v / ntot / R / t

      Grp = Grv - ntot * R * t * log(z)

      ! ========================================================================
      ! Equation 31
      ! ------------------------------------------------------------------------
      if (present(eq31)) eq31 = sum(n(:) * lnphi(:)) - Grp / (R * t)

      ! ========================================================================
      ! Equation 32
      ! ------------------------------------------------------------------------
      ! TODO

      ! ========================================================================
      ! Equation 33
      ! ------------------------------------------------------------------------
      if (present(eq33)) then
         do i = 1, size(n), 1
            do j = 1, size(n), 1
               eq33(i, j) = dlnPhidn(i,j) - dlnPhidn(j,i)
            end do
         end do
      end if

      ! ========================================================================
      ! Equation 34
      ! ------------------------------------------------------------------------
      if (present(eq34)) then
         eq34 = 0.0_pr

         do j = 1, size(n), 1
            do i = 1, size(n), 1
               eq34(j) = eq34(j) + n(i) * dlnPhidn(i,j)
            end do
         end do
      end if

      ! ========================================================================
      ! Equation 36
      ! ------------------------------------------------------------------------
      if (present(eq36)) eq36 = sum(n(:) * dlnPhidP(:)) - (z - 1) * ntot / p

      ! ========================================================================
      ! Equation 37
      ! ------------------------------------------------------------------------
      if (present(eq37)) then
         eq37 = sum(n(:) * dlnPhidT(:)) + Hrv / (R * t**2)
      end if
   end subroutine ar_consistency

   subroutine numeric_ar_derivatives(&
      eos, n, v, t, d_n, d_v, d_t, &
      Ar, ArV, ArT, Arn, ArV2, ArT2, ArTV, ArVn, ArTn, Arn2 &
      )
      !! # numeric_ar_derivatives
      !! Evaluate the Helmholtz derivatives with central finite difference.
      !!
      !! # Description
      !! Tool to facilitate the development of new [[ArModel]] by testing
      !! the implementation of analytic derivatives.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos, only: pr, SoaveRedlichKwong, ArModel
      !!  use yaeos__consistency_armodel, only: numeric_ar_derivatives
      !!
      !!  class(ArModel), allocatable :: model
      !!  real(pr) :: tc(4), pc(4), w(4)
      !!
      !!  real(pr) :: n(4), t, v
      !!
      !!  real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      !!  real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      !!  real(pr) :: Arn2_num(size(n), size(n))
      !!
      !!  n = [1.5, 0.2, 0.7, 2.3]
      !!  tc = [190.564, 425.12, 300.11, 320.25]
      !!  pc = [45.99, 37.96, 39.23, 40.21]
      !!  w = [0.0115478, 0.200164, 0.3624, 0.298]
      !!
      !!  t = 600_pr
      !!  v = 0.5_pr
      !!
      !!  model = SoaveRedlichKwong(tc, pc, w)
      !!
      !!  call numeric_ar_derivatives(&
      !!     model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
      !!     Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
      !!     ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
      !!     Arn2=Arn2_num &
      !!     )
      !! ```
      !!
      class(ArModel), intent(in) :: eos !! Equation of state
      real(pr), intent(in) :: n(:) !! Moles number vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: d_n !! Moles finite difference step
      real(pr), intent(in) :: d_t !! Temperature finite difference step
      real(pr), intent(in) :: d_v !! Volume finite difference step
      real(pr), intent(out) :: Ar !! Residual Helmoltz energy
      real(pr), optional, intent(out) :: ArV !! \(\frac{dAr}{dV}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
      real(pr), optional, intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
      real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
      real(pr), optional, intent(out) :: ArTV !! \(\frac{d^2Ar}{dTdV}\)
      real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVdn_i}\)
      real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTdn_i}\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\frac{d^2Ar}{dn_{ij}}\)

      ! Auxiliary
      real(pr) :: Ar_aux1, Ar_aux2, Ar_aux3, Ar_aux4
      real(pr) :: dn_aux1(size(n)), dn_aux2(size(n))
      integer :: i, j

      ! ========================================================================
      ! Ar valuations
      ! ------------------------------------------------------------------------
      ! on point valuation
      call eos%residual_helmholtz(n, v, t, Ar=Ar)

      ! ========================================================================
      ! Central numeric derivatives
      ! ------------------------------------------------------------------------
      ! Volume
      if (present(ArV) .or. present(ArV2)) then
         call eos%residual_helmholtz(n, v + d_v, t, Ar=Ar_aux1)
         call eos%residual_helmholtz(n, v - d_v, t, Ar=Ar_aux2)

         if (present(ArV)) ArV = (Ar_aux1 - Ar_aux2) / (2 * d_v)
         if (present(ArV2)) ArV2 = (Ar_aux1 - 2 * Ar + Ar_aux2) / d_v**2
      end if

      ! Temperature
      if (present(ArT) .or. present(ArT2)) then
         call eos%residual_helmholtz(n, v, t + d_t, Ar=Ar_aux1)
         call eos%residual_helmholtz(n, v, t - d_t, Ar=Ar_aux2)

         if (present(ArT)) ArT = (Ar_aux1 - Ar_aux2) / (2 * d_t)
         if (present(ArT2)) ArT2 = (Ar_aux1 - 2 * Ar + Ar_aux2) / d_t**2
      end if

      ! Mole first derivatives
      if (present(Arn)) then
         Arn = 0.0_pr

         do i = 1, size(n), 1
            dn_aux1 = 0.0_pr
            dn_aux1(i) = d_n

            call eos%residual_helmholtz(n + dn_aux1, v, t, Ar=Ar_aux1)
            call eos%residual_helmholtz(n - dn_aux1, v, t, Ar=Ar_aux2)

            Arn(i) = (Ar_aux1 - Ar_aux2) / (2 * d_n)
         end do
      end if

      ! ========================================================================
      ! Central cross derivatives
      ! ------------------------------------------------------------------------
      ! Temperature - Volume
      if (present(ArTV)) then
         call eos%residual_helmholtz(n, v + d_v, t + d_t, Ar=Ar_aux1)
         call eos%residual_helmholtz(n, v + d_v, t - d_t, Ar=Ar_aux2)
         call eos%residual_helmholtz(n, v - d_v, t + d_t, Ar=Ar_aux3)
         call eos%residual_helmholtz(n, v - d_v, t - d_t, Ar=Ar_aux4)

         ArTV = (Ar_aux1 - Ar_aux2 - Ar_aux3 + Ar_aux4) / (4 * d_t * d_v)
      end if

      ! Temperature - Mole
      if (present(ArTn)) then
         ArTn = 0.0_pr

         do i = 1, size(n), 1
            dn_aux1 = 0.0_pr
            dn_aux1(i) = d_n

            call eos%residual_helmholtz(n + dn_aux1, v, t + d_t, Ar=Ar_aux1)
            call eos%residual_helmholtz(n + dn_aux1, v, t - d_t, Ar=Ar_aux2)
            call eos%residual_helmholtz(n - dn_aux1, v, t + d_t, Ar=Ar_aux3)
            call eos%residual_helmholtz(n - dn_aux1, v, t - d_t, Ar=Ar_aux4)

            ArTn(i) = &
               (Ar_aux1 - Ar_aux2 - Ar_aux3 + Ar_aux4) / (4 * d_t * d_n)
         end do
      end if

      ! Volume - Mole
      if (present(ArVn)) then
         ArVn = 0.0_pr

         do i = 1, size(n), 1
            dn_aux1 = 0.0_pr
            dn_aux1(i) = d_n

            call eos%residual_helmholtz(n + dn_aux1, v + d_v, t, Ar=Ar_aux1)
            call eos%residual_helmholtz(n + dn_aux1, v - d_v, t, Ar=Ar_aux2)
            call eos%residual_helmholtz(n - dn_aux1, v + d_v, t, Ar=Ar_aux3)
            call eos%residual_helmholtz(n - dn_aux1, v - d_v, t, Ar=Ar_aux4)

            ArVn(i) = &
               (Ar_aux1 - Ar_aux2 - Ar_aux3 + Ar_aux4) / (4 * d_v * d_n)
         end do
      end if

      ! Mole second derivatives
      if (present(Arn2)) then
         Arn2 = 0.0_pr

         do i = 1, size(n), 1
            do j = 1, size(n), 1
               if (i .eq. j) then
                  dn_aux1 = 0.0_pr
                  dn_aux1(i) = d_n

                  call eos%residual_helmholtz(n + dn_aux1, v, t, Ar=Ar_aux1)
                  call eos%residual_helmholtz(n - dn_aux1, v, t, Ar=Ar_aux2)

                  Arn2(i, j) = (Ar_aux1 - 2 * Ar + Ar_aux2) / d_n**2
               else
                  dn_aux1 = 0.0_pr
                  dn_aux2 = 0.0_pr

                  dn_aux1(i) = d_n
                  dn_aux2(j) = d_n

                  call eos%residual_helmholtz(&
                     n + dn_aux1 + dn_aux2, v, t, Ar=Ar_aux1 &
                     )
                  call eos%residual_helmholtz(&
                     n + dn_aux1 - dn_aux2, v, t, Ar=Ar_aux2 &
                     )
                  call eos%residual_helmholtz(&
                     n - dn_aux1 + dn_aux2, v, t, Ar=Ar_aux3 &
                     )
                  call eos%residual_helmholtz(&
                     n - dn_aux1 - dn_aux2, v, t, Ar=Ar_aux4 &
                     )

                  Arn2(i, j) = &
                     (Ar_aux1 - Ar_aux2 - Ar_aux3 + Ar_aux4) / (4 * d_n**2)
               end if
            end do
         end do
      end if
   end subroutine numeric_ar_derivatives
end module yaeos__consistency_armodel
