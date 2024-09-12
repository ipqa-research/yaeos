module yaeos__models_ge_uniquac
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   implicit none

   type, extends(GeModel) :: UNIQUAC
      real(pr) :: z = 10.0_pr
      !! Model coordination number
      real(pr), allocatable :: qs(:)
      !! Molecule's relative areas \(Q_i\)
      real(pr), allocatable :: rs(:)
      !! Molecule's relative volumes \(R_i\)
      real(pr), allocatable :: aij(:,:)
      !! Interaction parameters matrix \(a_{ij}\)
      real(pr), allocatable :: bij(:,:)
      !! Interaction parameters matrix \(b_{ij}\)
      real(pr), allocatable :: cij(:,:)
      !! Interaction parameters matrix \(c_{ij}\)
      real(pr), allocatable :: dij(:,:)
      !! Interaction parameters matrix \(d_{ij}\)
      real(pr), allocatable :: eij(:,:)
      !! Interaction parameters matrix \(e_{ij}\)
   contains
      procedure :: excess_gibbs
      procedure :: taus
   end type UNIQUAC

contains
   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      class(UNIQUAC), intent(in) :: self
      !! UNIQUAC model
      real(pr), intent(in) :: n(:)
      !! Moles vector [mol]
      real(pr), intent(in) :: T
      !! Temperature [K]
      real(pr), optional, intent(out) :: Ge
      !! Excess Gibbs energy
      real(pr), optional, intent(out) :: GeT
      !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2
      !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n))
      !! \(\frac{dG^E}{dn}\)
      real(pr), optional, intent(out) :: GeTn(size(n))
      !! \(\frac{d^2G^E}{dTdn}\)
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))
      !! \(\frac{d^2G^E}{dn^2}\)

      real(pr) :: theta_i(size(n)), dtheta_i_dn(size(n), size(n))
      real(pr) :: d2theta_i_dn(size(n), size(n))
      real(pr) :: phi_i(size(n)), dphi_i_dn(size(n), size(n))
      real(pr) :: d2phi_i_dn(size(n), size(n))

      integer :: i, j

      logical :: dt, dt2, dn, dtn, dn2

      real(pr) :: Ge_comb, Ge_res
      real(pr) :: tau_ij(size(n), size(n)), dtau_ij(size(n), size(n))
      real(pr) :: d2tau_ij(size(n), size(n))

      ! Auxiliars
      integer :: nc
      real(pr) :: n_tot
      real(pr) :: xi(size(n))
      real(pr) :: sum_niqi, sum_niri

      real(pr) :: log_phi_i_xi(size(n))
      real(pr) :: sum_theta_j_tau_ji(size(n))
      real(pr) :: sum_theta_j_dtau_ji(size(n))
      real(pr) :: sum_theta_j_d2tau_ji(size(n))

      ! Auxiliars for compositionals derivatives
      real(pr) :: dxi_dnj(size(n), size(n))
      real(pr) :: Gen_comb(size(n)), Gen_res(size(n))
      real(pr) :: Gen_aux(size(n))

      real(pr) :: Ge_aux, GeT_aux, GeT2_aux, diff_aux(size(n))

      ! =======================================================================
      ! Logical variables
      ! -----------------------------------------------------------------------
      dt = present(GeT)
      dt2 = present(GeT2)
      dn = present(Gen)
      dtn = present(GeTn)
      dn2 = present(Gen2)

      ! =======================================================================
      ! Auxiliars
      ! -----------------------------------------------------------------------
      nc = size(n)

      n_tot = sum(n)

      xi = n / n_tot

      sum_niqi = sum(n * self%qs)
      sum_niri = sum(n * self%rs)

      ! =======================================================================
      ! tau call (temperature dependence term)
      ! -----------------------------------------------------------------------
      if (dt .and. .not. dt2) call self%taus(T, tau_ij, dtau_ij)
      if (.not. dt .and. dt2) call self%taus(T, tau_ij, tauT2=d2tau_ij)
      if (dt .and. dt2) call self%taus(T, tau_ij, dtau_ij, d2tau_ij)
      if (.not. dt .and. .not. dt2) call self%taus(T, tau_ij)


      ! =======================================================================
      ! theta_i
      ! -----------------------------------------------------------------------
      theta_i = n * self%qs / sum_niqi

      if (dn .or. dtn) then
         dtheta_i_dn = 0
         do concurrent(i=1:nc, j=1:nc)
            if (i == j) then
               dtheta_i_dn(i,i) = &
                  (self%qs(i) * sum_niqi - n(i) * self%qs(i)**2) / sum_niqi**2
            else
               dtheta_i_dn(i,j) = -n(i) * self%qs(i) * self%qs(j) / sum_niqi**2
            end if
         end do
      end if

      if (dn2) then
         d2theta_i_dn = 0
         do concurrent(i=1:nc, j=1:nc)
            if (i == j) then
               d2theta_i_dn(i,i) = &
                  (&
                  -2.0_pr * (self%qs(i) * sum_niqi &
                  - n(i) * self%qs(i)**2) * self%qs(i) * sum_niqi &
                  ) &
                  / sum_niqi**4
            else
               d2theta_i_dn(i,j) = &
                  (&
                  - self%qs(i) * self%qs(j) * sum_niqi**2 &
                  + 2.0_pr * n(i) * self%qs(i)**2 * self%qs(j) * sum_niqi &
                  ) &
                  / sum_niqi**4
            end if
         end do
      end if

      ! =======================================================================
      ! phi_i
      ! -----------------------------------------------------------------------
      phi_i = n * self%rs / sum_niri

      if (dn .or. dtn) then
         dphi_i_dn = 0
         do concurrent(i=1:nc, j=1:nc)
            if (i == j) then
               dphi_i_dn(i,i) = &
                  (self%rs(i) * sum_niri - n(i) * self%rs(i)**2) / sum_niri**2
            else
               dphi_i_dn(i,j) = -n(i) * self%rs(i) * self%rs(j) / sum_niri**2
            end if
         end do
      end if

      if (dn2) then
         d2phi_i_dn = 0
         do concurrent(i=1:nc, j=1:nc)
            if (i == j) then
               d2phi_i_dn(i,i) = &
                  (&
                  -2.0_pr * (self%rs(i) * sum_niri &
                  - n(i) * self%rs(i)**2) * self%rs(i) * sum_niri &
                  ) &
                  / sum_niri**4
            else
               d2phi_i_dn(i,j) = &
                  (&
                  - self%rs(i) * self%rs(j) * sum_niri**2 &
                  + 2.0_pr * n(i) * self%rs(i)**2 * self%rs(j) * sum_niri &
                  ) &
                  / sum_niri**4
            end if
         end do
      end if

      ! =======================================================================
      ! Ge
      ! -----------------------------------------------------------------------
      ! Combinatorial term
      log_phi_i_xi = log(phi_i * xi)

      Ge_comb = R * T * ( &
         sum(n * log_phi_i_xi) &
         + self%z / 2.0_pr * sum(n * self%qs * (log(theta_i) - log(phi_i))) &
         )

      ! Residual term
      sum_theta_j_tau_ji = 0.0_pr

      do i=1,nc
         sum_theta_j_tau_ji(i) = sum(theta_i * tau_ij(i,:))
      end do

      Ge_res = - R * T * sum(n * self%qs * log(sum_theta_j_tau_ji))

      Ge_aux = Ge_comb + Ge_res

      ! =======================================================================
      ! Derivatives
      ! -----------------------------------------------------------------------
      ! Compositional derivarives
      ! dn
      if (dn) then
         ! Mole fraction derivatives
         dxi_dnj = 0.0_pr

         do concurrent (i=1:nc, j=1:nc)
            if (i == j) then
               dxi_dnj(i,j) = (n_tot - n(i)) / n_tot**2
            else
               dxi_dnj(i,j) = -n(i) / n_tot**2
            end if
         end do

         ! Combinatorial term
         do i=1,nc
            Gen_comb(i) = (&
               log(phi_i(i) * n_tot / n(i)) &
               + sum(n * (dphi_i_dn(i,:) / phi_i - dxi_dnj(i,:) * n_tot / n)) &
               + self%z / 2.0_pr * ( &
               self%qs(i) * (log(theta_i(i)) - log(phi_i(i))) &
               + sum(self%qs * n * (dtheta_i_dn(i,:) / theta_i(i) - dphi_i_dn(i,:) / phi_i)) &
               ) &
            )
         end do

         ! Residual term
         do i = 1,nc
            Gen_res(i) = - (&
               sum(self%qs * log(sum_theta_j_tau_ji)) &
               
            )
         end do

         Gen_aux = R * T * (Gen_comb + Gen_res)
      end if

      ! Temperature derivatives
      if (dt .or. dt2) then
         sum_theta_j_dtau_ji = 0.0_pr

         do i=1,nc
            sum_theta_j_dtau_ji(i) = sum(theta_i * dtau_ij(i,:))
         end do

         GeT_aux = ( &
            Ge_aux / T &
            -R * T * sum( &
            self%qs * n / n_tot * sum_theta_j_dtau_ji / sum_theta_j_tau_ji &
            ) &
            )
      end if

      if (dt2) then
         sum_theta_j_d2tau_ji = 0.0_pr

         do i=1,nc
            sum_theta_j_d2tau_ji(i) = sum(theta_i * d2tau_ij(i,:))
         end do

         diff_aux = (&
            sum_theta_j_d2tau_ji / sum_theta_j_tau_ji &
            - (sum_theta_j_dtau_ji / sum_theta_j_tau_ji)**2 &
            )

         GeT2_aux = -R * ( &
            T * sum(self%qs * n / n_tot * diff_aux) &
            + 2.0_pr * sum(self%qs * n / n_tot * sum_theta_j_dtau_ji / sum_theta_j_tau_ji) &
            )
      end if
      ! =======================================================================
      ! Excess Gibbs energy returns
      ! -----------------------------------------------------------------------
      if (present(Ge)) Ge = Ge_aux
      if (dt) GeT = GeT_aux
      if (dt2) GeT2 = GeT2_aux
      if (dn) Gen = Gen_aux
      if (dtn) GeTn = 0.0_pr
      if (dn2) Gen2 = 0.0_pr
   end subroutine excess_gibbs

   subroutine taus(self, T, tau, tauT, tauT2)
      class(UNIQUAC), intent(in) :: self
      !! UNIQUAC model
      real(pr), intent(in) :: T
      !! Temperature [K]
      real(pr), optional, intent(out) :: tau(size(self%qs), size(self%qs))
      !! UNIQUAC temperature dependence term
      real(pr), optional, intent(out) :: tauT(size(self%qs), size(self%qs))
      !! \(\frac{d\tau_{ij}}{dT}\)
      real(pr), optional, intent(out) :: tauT2(size(self%qs), size(self%qs))
      !! \(\frac{d^2\tau_{ij}}{dT^2}\)

      ! aux
      real(pr) :: tau_aux(size(self%qs), size(self%qs))
      real(pr) :: u(size(self%qs), size(self%qs))
      real(pr) :: du(size(self%qs), size(self%qs))
      real(pr) :: d2u(size(self%qs), size(self%qs))

      ! Logical
      logical :: tt, dt, dt2

      tt = present(tau)
      dt = present(tauT)
      dt2 = present(tauT2)

      ! temperature function
      u = self%aij + self%bij/T + self%cij*log(T) + self%dij*T + self%eij*T**2

      ! tau_ij
      tau_aux = exp(u)

      ! dT
      if (dt .or. dt2) then
         du = -self%bij / T**2 + self%cij / T + self%dij + 2.0_pr * self%eij * T
      end if

      ! d2T
      if (dt2) then
         d2u = 2.0_pr * self%bij / T**3 - self%cij / T**2 + 2.0_pr * self%eij
      end if


      if (tt) tau = tau_aux
      if (dt) tauT = tau_aux * du
      if (dt2) tauT2 = tau_aux * du**2 + tau_aux * d2u
   end subroutine taus

   type(UNIQUAC) function setup_uniquac(qs, rs, aij, bij, cij, dij, eij)
      real(pr), intent(in) :: qs(:)
      real(pr), intent(in) :: rs(:)
      real(pr), optional, intent(in) :: aij(:,:)
      real(pr), optional, intent(in) :: bij(:,:)
      real(pr), optional, intent(in) :: cij(:,:)
      real(pr), optional, intent(in) :: dij(:,:)
      real(pr), optional, intent(in) :: eij(:,:)

      ! aij
      if (present(aij)) then
         setup_uniquac%aij = aij
      else
         setup_uniquac%aij = 0.0_pr
      end if

      ! bij
      if (present(bij)) then
         setup_uniquac%bij = bij
      else
         setup_uniquac%bij = 0.0_pr
      end if

      ! cij
      if (present(cij)) then
         setup_uniquac%cij = cij
      else
         setup_uniquac%cij = 0.0_pr
      end if

      ! dij
      if (present(dij)) then
         setup_uniquac%dij = dij
      else
         setup_uniquac%dij = 0.0_pr
      end if

      ! eij
      if (present(eij)) then
         setup_uniquac%eij = eij
      else
         setup_uniquac%eij = 0.0_pr
      end if

      setup_uniquac%qs = qs
      setup_uniquac%rs = rs
   end function setup_uniquac
end module yaeos__models_ge_uniquac
