module yaeos_consistency_armodel
   !! Consistency checks of Helmholtz free energy models.

   !! This module contains tools to validate the analityc derivatives of !
   !! implmented Helmholtz free energy models (ArModel).
   !!
   !! Available tools:
   !!
   !! - numeric_ar_derivatives: From an instantiated ArModel evaluate all the
   !! Helmholtz free energy derivatives from the central finite difference
   !! method.
   use yaeos_constants, only: pr, R
   use yaeos_models_ar, only: ArModel

   implicit none
contains
   subroutine numeric_ar_derivatives(&
      eos, n, v, t, d_n, d_v, d_t, &
      Ar, ArV, ArT, Arn, ArV2, ArT2, ArTV, ArVn, ArTn, Arn2 &
      )
      !! Evaluate the Helmholtz derivatives with central finite difference.
      !!
      class(ArModel), intent(in) :: eos !! Model
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
      call eos%residual_helmholtz(n, v + d_v, t, Ar=Ar_aux1)
      call eos%residual_helmholtz(n, v - d_v, t, Ar=Ar_aux2)

      if (present(ArV)) ArV = (Ar_aux1 - Ar_aux2) / (2 * d_v)
      if (present(ArV2)) ArV2 = (Ar_aux1 - 2 * Ar + Ar_aux2) / d_v**2

      ! Temperature
      call eos%residual_helmholtz(n, v, t + d_t, Ar=Ar_aux1)
      call eos%residual_helmholtz(n, v, t - d_t, Ar=Ar_aux2)

      if (present(ArT)) ArT = (Ar_aux1 - Ar_aux2) / (2 * d_t)
      if (present(ArT2)) ArT2 = (Ar_aux1 - 2 * Ar + Ar_aux2) / d_t**2

      ! Mole first derivatives
      if (present(Arn)) Arn = 0.0_pr

      do i = 1, size(n), 1
         dn_aux1 = 0.0_pr
         dn_aux1(i) = d_n

         call eos%residual_helmholtz(n + dn_aux1, v, t, Ar=Ar_aux1)
         call eos%residual_helmholtz(n - dn_aux1, v, t, Ar=Ar_aux2)

         if (present(Arn)) Arn(i) = (Ar_aux1 - Ar_aux2) / (2 * d_n)
      end do


      ! ========================================================================
      ! Central cross derivatives
      ! ------------------------------------------------------------------------
      ! Temperature - Volume
      call eos%residual_helmholtz(n, v + d_v, t + d_t, Ar=Ar_aux1)
      call eos%residual_helmholtz(n, v + d_v, t - d_t, Ar=Ar_aux2)
      call eos%residual_helmholtz(n, v - d_v, t + d_t, Ar=Ar_aux3)
      call eos%residual_helmholtz(n, v - d_v, t - d_t, Ar=Ar_aux4)

      if (present(ArTV)) ArTV = &
         (Ar_aux1 - Ar_aux2 - Ar_aux3 + Ar_aux4) / (4 * d_t * d_v)

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

end module yaeos_consistency_armodel
