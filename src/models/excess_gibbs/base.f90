module yaeos__models_ge_base
   !! Base module for excess Gibbs energy models in YAEOS
   !! This module provides implementations of excess Gibbs energy models
   !! using only Fortran native types.
   use yaeos__constants, only: pr, R
contains

   ! ==========================================================================
   ! NRTL Model, modified by Huron and vidal
   ! --------------------------------------------------------------------------

   subroutine nrtl_hv_ge(&
      n, T, &
      b, &
      alpha, &
      tau, dtaudt, dtaudt2, &
      Ge, Gen, GeT, GeT2, GeTn, Gen2)
      real(pr), intent(in) :: n(:) !! Number of moles vector.
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: b(:) !! \(b\) parameter
      real(pr), intent(in) :: alpha(:, :) !! \(alpha\) matrix
      real(pr), intent(in) :: tau(:, :) !! \(\tau\) matrix
      real(pr), intent(in) :: dtaudt(:, :)
      !! First derivative of \(\tau\) with respect to temperature
      real(pr), intent(in) :: dtaudt2(:, :)
      !! Second derivative of \(\tau\) with respect to temperature
      real(pr), optional, intent(out) :: Ge
      !! Excess Gibbs energy
      real(pr), optional, intent(out) :: Gen(:)
      !! Excess Gibbs energy derivative with repect to number of moles
      real(pr), optional, intent(out) :: GeT
      !! Excess Gibbs energy derivative with respect to temperature
      real(pr), optional, intent(out) :: GeT2
      !! Excess Gibbs energy second derivative with respect to temperature
      real(pr), optional, intent(out) :: GeTn(:)
      !! Excess Gibbs energy derivative with respect to temperature and number of moles
      real(pr), optional, intent(out) :: Gen2(:, :)
      !! Excess Gibbs energy second derivative with respect to number of moles

      real(pr) :: E(size(n), size(n)), U, D
      real(pr) :: dEdT(size(n), size(n)), dEdT2(size(n), size(n))

      real(pr) :: Dn(size(n))

      real(pr) :: xi(size(n), size(n)), theta(size(n), size(n)), omega(size(n), size(n)), eta(size(n), size(n))
      real(pr) :: xiT(size(n)), etaT(size(n), size(n)), thetaT(size(n)), omegaT(size(n), size(n))
      real(pr) :: xiTT(size(n)), thetaTT(size(n))

      real(pr) :: denom

      integer :: i, j, k, l, nc

      logical :: p_ge, p_get, p_get2, p_gent, p_gen2, p_gen
      real(pr) :: aux_Ge, aux_Gen(size(n)), aux_GeT, aux_GeTn(size(n))

      integer :: m

      nc = size(n)
      p_ge = present(Ge)
      p_get = present(GeT)
      p_gen = present(Gen)
      p_get2 = present(GeT2)
      p_gent = present(GeTn)
      p_gen2 = present(Gen2)

      E = exp(-alpha * tau)
      dEdT = - alpha * dtaudt * E
      dEdT2 = -alpha * (dtaudt2 * E + dtaudt * dEdT)

      do i=1,nc
         xi(:, i)     = E(:, i) * b * tau(:, i) * n
         eta(:, i)    = E(:, i) * b * tau(:, i)
         theta(:, i)  = E(:, i) * b * n
         omega(:, i)  = E(:, i) * b

         xiT(i)    = sum(b * n * (dEdT(:, i) * tau(:, i) + E(:, i)*dtaudt(:, i)))
         etaT(:, i)   = b * (dEdT(:, i) * tau(:, i) + E(:, i) * dtaudt(:, i))
         thetaT(i) = sum(dEdT(:, i) * b * n)
         omegaT(:, i) = dEdT(:, i) * b

         xiTT(i) = sum(&
            b * n * &
            (dEdT2(:, i) * tau(:, i) &
            + 2 * dEdT(:, i) * dtaudt(:, i) &
            + E(:, i) * dtaudt2(:, i)))
         thetaTT(i) = sum(b * n * (dEdT2(:, i)))
      end do

      do i=1,nc
         Dn(i) = sum(theta(:, i))
      end do

      ! ==============================================================
      ! Calculate the excess Gibbs energy
      ! It is also needed to calculate the derivatives wrt T
      ! --------------------------------------------------------------
      if (p_ge .or. p_get .or. p_get2) then
         aux_Ge = 0
         do i=1,nc
            aux_Ge = aux_Ge + n(i) * sum(xi(:, i) / sum(theta(:, i)))
         end do
      end if


      ! ==============================================================
      ! Derivatives wrt number of moles
      ! --------------------------------------------------------------
      if (p_gen .or. p_gent) then
         aux_Gen = 0.0
         do i=1,nc
            aux_Gen(i) = sum(xi(:, i)) / sum(theta(:, i))
            do k=1,nc
               aux_Gen(i) = aux_Gen(i) + n(k) * (&
                  eta(i, k)/sum(theta(:, k)) &
                  - omega(i, k) * sum(xi(:, k))/sum(theta(:,k))**2 &
                  )
            end do
         end do
      end if

      if (p_gen2) then
         Gen2 = 0.0
         do i=1,nc
            do j=1,nc
               Gen2(i, j) = &
                  - omega(j, i) * sum(xi(:, i)) / sum(theta(:, i))**2 &
                  - omega(i, j) * sum(xi(:, j)) / sum(theta(:, j))**2 &
                  + eta(j, i) / sum(theta(:, i)) &
                  + eta(i, j) / sum(theta(:, j))

               do k=1,nc
                  denom = sum(theta(:, k))
                  Gen2(i, j) = Gen2(i, j) + &
                     2 * n(k) * omega(i, k) * omega(j, k) * sum(xi(:, k)) / denom**3 &
                     - n(k) * omega(i, k) * eta(j, k) / denom**2 &
                     - n(k) * omega(j, k) * eta(i, k) / denom**2
               end do
            end do
         end do
      end if

      if (p_genT) then
         do i=1,nc
            aux_GeTn(i) = xiT(i)/Dn(i) - sum(xi(:, i)) * thetaT(i)/Dn(i)**2
            do k=1,nc
               aux_GeTn(i) = aux_GeTn(i) + n(k) * (&
                  etaT(i, k)/Dn(k) - eta(i, k) * thetaT(k) / Dn(k)**2 &
                  - (omega(i, k) * xiT(k) + omegaT(i, k) * sum(xi(:, k))) / Dn(k)**2 &
                  + 2 * thetaT(k) * omega(i, k) * sum(xi(:, k)) / Dn(k)**3 &
                  )
            end do
         end do
      end if

      ! ==============================================================
      ! Derivatives wrt temperature
      ! --------------------------------------------------------------
      if (p_get .or. p_get2) then
         aux_GeT = 0
         do i=1,nc
            aux_GeT = aux_GeT + n(i) * (xiT(i)/sum(theta(:, i)) - sum(xi(:, i))*(thetaT(i))/Dn(i)**2)
         end do
      end if

      if (p_get2) then
         GeT2 = 0
         do i=1,nc
            GeT2 = GeT2 + n(i) * (&
               xiTT(i)/Dn(i) - 2*xiT(i)*thetaT(i)/Dn(i)**2 &
               - sum(xi(:, i))*thetaTT(i)/Dn(i)**2 + 2*sum(xi(:, i))*thetaT(i)**2/Dn(i)**3)
         end do
      end if

      if (present(Ge))   Ge   = aux_Ge * R * T
      if (present(Gen))  Gen  = aux_Gen * R * T
      if (present(GeTn)) GeTn = (R*T*aux_Gen)/T + R * T * aux_GeTn
      if (present(Gen2)) Gen2 = Gen2 * R * T
      if(present(GeT))   GeT  = aux_GeT * R * T + (aux_Ge * R * T) / T
      if(present(Get2))  GeT2 = -(aux_Ge*R*T)/T**2 + (aux_GeT * R * T + (aux_Ge*R*T) / T)/T + R * aux_GeT + R*T * GeT2
   end subroutine nrtl_hv_ge

   elemental subroutine nrtl_hv_tdep(T, gij, tau, dtaudt, dtaudt2)
      !! Temperature dependent parameters for NRTL model
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: gij !! Interaction parameters
      real(pr), intent(out) :: tau !! \(\tau\) matrix
      real(pr), intent(out) :: dtaudt
      !! First derivative of \(\tau\) with respect to temperature
      real(pr), intent(out) :: dtaudt2
      !! Second derivative of \(\tau\) with respect to temperature

      tau = gij/(R*T)
      dtaudt = -gij/(R*T**2)
      dtaudt2 = 2*gij/(R*T**3)
   end subroutine nrtl_hv_tdep

   elemental subroutine nrtl_hv_tdep_linear(T, A, B, tau, dtaudt, dtaudt2)
      !! Temperature dependent parameters for NRTL model
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: A !! Interaction parameters
      real(pr), intent(in) :: B !! Interaction parameters
      real(pr), intent(out) :: tau !! \(\tau\) matrix
      real(pr), intent(out) :: dtaudt
      !! First derivative of \(\tau\) with respect to temperature
      real(pr), intent(out) :: dtaudt2
      !! Second derivative of \(\tau\) with respect to temperature

      tau = A + B/T
      dtaudt = -B/(T**2)
      dtaudt2 = 2*B/(T**3)
   end subroutine nrtl_hv_tdep_linear
end module yaeos__models_ge_base
