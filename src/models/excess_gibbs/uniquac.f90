module yaeos__models_ge_uniquac
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   implicit none

   type, extends(GeModel) :: UNIQUAC
      real(pr) :: z = 10.0
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

      real(pr) :: thetak(size(n))
      real(pr) :: dthetak_dni(size(n), size(n))
      real(pr) :: d2thetak_dnidnj(size(n), size(n), size(n))

      real(pr) :: phik(size(n))
      real(pr) :: dphik_dn(size(n), size(n))
      real(pr) :: d2phik_dnidnj(size(n), size(n), size(n))

      real(pr) :: tau(size(n), size(n))
      real(pr) :: dtau(size(n), size(n))
      real(pr) :: d2tau(size(n), size(n))

      real(pr) :: Ge_comb, Ge_res

      integer :: i, j, k, l

      logical :: dt, dt2, dn, dtn, dn2

      ! Auxiliars
      integer :: nc
      real(pr) :: n_tot
      real(pr) :: xk(size(n))
      real(pr) :: r_i, q_i, r_j, q_j, r_k, q_k
      real(pr) :: sum_nq, sum_nr
      real(pr) :: Ge_aux

      real(pr) :: sum_thetal_tau_lk(size(n))
      real(pr) :: sum_theta_j_dtau_ji(size(n))
      real(pr) :: sum_theta_j_d2tau_ji(size(n))

      ! Auxiliars for compositionals derivatives
      real(pr) :: dxi_dnj(size(n), size(n))
      real(pr) :: d2xk_dnidnj(size(n), size(n), size(n))
      real(pr) :: Gen_comb(size(n)), Gen_res(size(n))
      real(pr) :: Gen_aux(size(n))
      real(pr) :: Gen2_aux(size(n), size(n))
      real(pr) :: Gen2_comb_term1(size(n), size(n))
      real(pr) :: Gen2_comb_term2(size(n), size(n))
      real(pr) :: Gen2_res(size(n), size(n))
      real(pr) :: sum_dthetak_dni_taulk(size(n), size(n))
      real(pr) :: sum_d2thetak_dnidnj_taulk(size(n))

      ! Temperature derivatives
      real(pr) :: GeT_aux, GeT2_aux, diff_aux(size(n))

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

      xk = n / n_tot

      sum_nq = sum(n * self%qs)
      sum_nr = sum(n * self%rs)

      ! =======================================================================
      ! tau call (temperature dependence term)
      ! -----------------------------------------------------------------------
      if (dt .and. .not. dt2) call self%taus(T, tau, dtau)
      if (.not. dt .and. dt2) call self%taus(T, tau, tauT2=d2tau)
      if (dt .and. dt2) call self%taus(T, tau, dtau, d2tau)
      if (.not. dt .and. .not. dt2) call self%taus(T, tau)


      ! =======================================================================
      ! theta_k
      ! -----------------------------------------------------------------------
      thetak = n * self%qs / sum_nq

      if (dn .or. dtn) then
         dthetak_dni = 0
         do concurrent(k=1:nc, i=1:nc)
            if (i == k) then
               dthetak_dni(i,i) = &
                  (self%qs(i) * sum_nq - n(i) * self%qs(i)**2) / sum_nq**2
            else
               dthetak_dni(k,i) = -n(k) * self%qs(i) * self%qs(k) / sum_nq**2
            end if
         end do
      end if

      if (dn2) then
         d2thetak_dnidnj = 0
         do concurrent(k=1:nc, i=1:nc, j=1:nc)
            if (i==k .and. j==k) then
               q_i = self%qs(i)

               d2thetak_dnidnj(k,i,j) = (&
                  2.0_pr * (q_i**3 * n(i) - q_i**2 * sum_nq) / sum_nq**3 &
                  )
            else if (i==k) then
               q_i = self%qs(i)
               q_j = self%qs(j)

               d2thetak_dnidnj(k,i,j) = (&
                  (2.0_pr * n(i) * q_i**2 * q_j - q_i*q_j*sum_nq) / sum_nq**3 &
                  )
            else if (j==k) then
               q_i = self%qs(i)
               q_j = self%qs(j)

               d2thetak_dnidnj(k,i,j) = (&
                  (2.0_pr * n(j) * q_j**2 * q_i - q_i*q_j*sum_nq) / sum_nq**3 &
                  )
            else
               q_i = self%qs(i)
               q_j = self%qs(j)
               q_k = self%qs(k)

               d2thetak_dnidnj(k,i,j) = (&
                  2.0_pr * n(k) * q_k * q_i * q_j / sum_nq**3 &
                  )
            end if
         end do
      end if

      ! =======================================================================
      ! phi_k
      ! -----------------------------------------------------------------------
      phik = n * self%rs / sum_nr

      if (dn .or. dtn) then
         dphik_dn = 0
         do concurrent(k=1:nc, i=1:nc)
            if (i == k) then
               dphik_dn(i,i) = &
                  (-n(i) * self%rs(i)**2 + self%rs(i) * sum_nr) / sum_nr**2
            else
               dphik_dn(k,i) = -n(i) * self%rs(i) * self%rs(k) / sum_nr**2
            end if
         end do
      end if

      if (dn2) then
         d2phik_dnidnj = 0
         do concurrent(k=1:nc, i=1:nc, j=1:nc)
            if (i==k .and. j==k) then
               r_i = self%rs(i)

               d2phik_dnidnj(k,i,j) = (&
                  2.0_pr * (r_i**3 * n(i) - r_i**2 * sum_nr) / sum_nr**3 &
                  )
            else if (i==k) then
               r_i = self%rs(i)
               r_j = self%rs(j)

               d2phik_dnidnj(k,i,j) = (&
                  (2.0_pr * n(i) * r_i**2 * r_j - r_j*r_i*sum_nr) / sum_nr**3 &
                  )
            else if (j==k) then
               r_i = self%rs(i)
               r_j = self%rs(j)

               d2phik_dnidnj(k,i,j) = (&
                  (2.0_pr * n(j) * r_j**2 * r_i - r_j*r_i*sum_nr) / sum_nr**3 &
                  )
            else
               r_i = self%rs(i)
               r_j = self%rs(j)
               r_k = self%rs(k)

               d2phik_dnidnj(k,i,j) = (&
                  2.0_pr * n(k) * r_k * r_i * r_j / sum_nr**3 &
                  )
            end if
         end do
      end if

      ! =======================================================================
      ! Ge
      ! -----------------------------------------------------------------------
      ! Combinatorial term
      Ge_comb = R * T * ( &
         sum(n * log(phik / xk)) &
         + self%z / 2.0_pr * sum(n * self%qs * (log(thetak) - log(phik))) &
         )

      ! Residual term
      sum_thetal_tau_lk = 0.0_pr

      do k=1,nc
         sum_thetal_tau_lk(k) = sum(thetak * tau(:,k))
      end do

      Ge_res = - R * T * sum(n * self%qs * log(sum_thetal_tau_lk))

      Ge_aux = Ge_comb + Ge_res

      ! =======================================================================
      ! Ge Derivatives
      ! -----------------------------------------------------------------------
      ! Compositional derivarives
      ! dn
      if (dn) then
         dxi_dnj = 0

         do concurrent(i=1:nc, j=1:nc)
            if (i == j) then
               dxi_dnj(i,i) = (n_tot - n(i)) / n_tot**2
            else
               dxi_dnj(i,j) = -n(i) / n_tot**2
            end if
         end do

         ! Combinatorial term
         do i=1,nc
            Gen_comb(i) = ( &
               log(phi_i(i) / xi(i)) + sum(n * (dphi_i_dn(:, i) / phi_i - dxi_dnj(:, i) / xi)) &
               + self%z / 2.0_pr * (self%qs(i) * log(theta_i(i) / phi_i(i)) + &
               sum(self%qs * n * (dtheta_i_dn(:,i) / theta_i - dphi_i_dn(:,i) / phi_i))) &
               )
         end do

         ! Residual term
         do i=1,nc
            Gen_res(i) = -(&
               self%qs(i) * log(sum_thetal_tau_lk(i)) + &
               sum(self%qs * n * sum(dtheta_i_dn(i,:) * tau_ij(i,:)) / sum_thetal_tau_lk(i)) &
               )
         end do

         Gen_aux = R * T * (Gen_comb + Gen_res)
      end if

      ! dn2
      if (dn2) then
         ! God help me
         d2xk_dnidnj = 0

         do concurrent (k=1:nc, i=1:nc, j=1:nc)
            if (i==k .and. j==k) then
               d2xk_dnidnj(k,i,j) = -2 * (n_tot - n(i)) / n_tot**3
            else if (i==k) then
               d2xk_dnidnj(k,i,j) = (2 * n(i) - n_tot) / n_tot**3
            else if (j==k) then
               d2xk_dnidnj(k,i,j) = (2 * n(j) - n_tot) / n_tot**3
            else
               d2xk_dnidnj(k,i,j) = 2 * n(k) / n_tot**3
            end if
         end do


         sum_d2thetak_dnidnj_taulk = 0.0_pr

         do k=1,nc
            do concurrent(i=1:nc, j=1:nc)
               sum_d2thetak_dnidnj_taulk(k) = (&
                  sum_d2thetak_dnidnj_taulk(k) + sum(d2thetak_dnidnj(:,i,j) * tau_ij(k,:)) &
                  )
            end do
         end do


         sum_dthetak_dni_taulk = 0.0_pr

         do concurrent(k=1:nc, i=1:nc)
            sum_dthetak_dni_taulk(i,k) = sum(dtheta_i_dn(i,:) * tau_ij(k,:))
         end do

         do concurrent(i=1:nc, j=1:nc)
            ! Combinatorial term
            ! first term
            Gen2_comb_term1(i,j) = (&
               dphi_i_dn(i,j) / phi_i(i) - dxi_dnj(i,j) / xi(i) &
               + dphi_i_dn(j,i) / phi_i(j) - dxi_dnj(j,i) / xi(j) &
               + sum(&
               n * &
               (d2phik_dnidnj(:,i,j) / phi_i - dphi_i_dn(:,i) * dphi_i_dn(:,j) / phi_i**2 &
               - d2xk_dnidnj(:,i,j) / xi - dxi_dnj(:,i) * dxi_dnj(:,j) / xi**2) &
               ) &
               )

            ! second term
            Gen2_comb_term2(i,j) = (&
               self%z / 2.0_pr * ( &
               self%qs(i) * (dtheta_i_dn(i,j) / theta_i(i) - dphi_i_dn(i,j) / phi_i(i)) &
               + self%qs(j) * (dtheta_i_dn(j,i) / theta_i(j) - dphi_i_dn(j,i) / phi_i(j)) &
               + sum(&
               self%qs * n * (&
               d2thetak_dnidnj(:,i,j) / theta_i - dtheta_i_dn(:,i) * dtheta_i_dn(:,j) / theta_i**2 &
               - d2phik_dnidnj(:,i,j) / phi_i - dphi_i_dn(:,i) * dphi_i_dn(:,j) / phi_i**2) &
               ) &
               ) &
               )

            ! Residual term
            Gen2_res(i,j) = -(&
               self%qs(i) * sum(dtheta_i_dn(:,j) * tau_ij(i,:)) / sum_thetal_tau_lk(i) &
               + self%qs(j) * sum(dtheta_i_dn(:,i) * tau_ij(j,:)) / sum_thetal_tau_lk(j) &
               + sum(&
               self%qs * n * (&
               sum_d2thetak_dnidnj_taulk / sum_thetal_tau_lk &
               - sum_dthetak_dni_taulk(i, :) * sum_dthetak_dni_taulk(j, :) / sum_thetal_tau_lk**2 &
               ) &
               ) &
               )

            Gen2_aux(i,j) = R * T * (Gen2_comb_term1(i,j) + Gen2_comb_term2(i,j) + Gen2_res(i,j))
         end do

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
            self%qs * n / n_tot * sum_theta_j_dtau_ji / sum_thetal_tau_lk &
            ) &
            )
      end if

      if (dt2) then
         sum_theta_j_d2tau_ji = 0.0_pr

         do i=1,nc
            sum_theta_j_d2tau_ji(i) = sum(theta_i * d2tau_ij(i,:))
         end do

         diff_aux = (&
            sum_theta_j_d2tau_ji / sum_thetal_tau_lk &
            - (sum_theta_j_dtau_ji / sum_thetal_tau_lk)**2 &
            )

         GeT2_aux = -R * ( &
            T * sum(self%qs * n / n_tot * diff_aux) &
            + 2.0_pr * sum(self%qs * n / n_tot * sum_theta_j_dtau_ji / sum_thetal_tau_lk) &
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
      if (dn2) Gen2 = Gen2_aux
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
