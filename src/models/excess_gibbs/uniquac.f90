module yaeos__models_ge_uniquac
   !! # UNIQUAC module
   !! UNIQUAC (**uni**versal **qua**si**c**hemical) Excess Gibbs free energy
   !! model.
   !!
   !! ## References
   !! 1. Maurer, G., & Prausnitz, J. M. (1978). On the derivation and extension
   !! of the UNIQUAC equation. Fluid Phase Equilibria, 2(2), 91-99.
   !! 2. Gmehling, Jurgen, Barbel Kolbe, Michael Kleiber, and Jurgen Rarey.
   !! Chemical Thermodynamics for Process Simulation. 1st edition. Weinheim:
   !! Wiley-VCH, 2012.
   !! 3. Caleb Bell and Contributors (2016-2024). Thermo: Chemical properties
   !! component of Chemical Engineering Design Library (ChEDL)
   !! https://github.com/CalebBell/thermo.
   !!
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   use yaeos__math, only: derivative_dxk_dni, derivative_d2xk_dnidnj
   implicit none

   type, extends(GeModel) :: UNIQUAC
      !! # UNIQUAC model
      !!
      !! # Description
      !! UNIQUAC (**uni**versal **qua**si**c**hemical) Excess Gibbs free energy
      !! model.
      !!
      !! \[
      !! \frac{G^E}{RT} = \sum_k n_k \ln\frac{\phi_k}{x_k}
      !! + \frac{z}{2}\sum_k q_k n_k \ln\frac{\theta_k}{\phi_k}
      !! - \sum_k q_k n_k \ln\left(\sum_l \theta_l \tau_{kl} \right)
      !! \]
      !!
      !! With:
      !!
      !! \[ x_k = \frac{n_k}{\sum_l n_l} \]
      !!
      !! \[ \phi_k = \frac{r_k n_k}{\sum_l r_l n_l} \]
      !!
      !! \[ \theta_k = \frac{q_k n_k}{\sum_l q_l n_l} \]
      !!
      !! \[ \tau_{lk} = \exp \left[\frac{-\Delta U_{lk}}{R T} \right] \]
      !!
      !! \[
      !! \frac{-\Delta U_{lk}}{R T} = a_{lk}+\frac{b_{lk}}{T}+c_{lk}\ln T +
      !! d_{lk}T + e_{lk}{T^2}
      !! \]
      !!
      !! # Example
      !!
      !! ```fortran
      !! use yaeos, only: pr, setup_uniquac, UNIQUAC
      !!
      !! integer, parameter :: nc = 3
      !!
      !! real(pr) :: rs(nc), qs(nc)
      !! real(pr) :: b(nc, nc)
      !! real(pr) :: n(nc)
      !!
      !! real(pr) :: ln_gammas(nc), T
      !!
      !! type(UNIQUAC) :: model
      !!
      !! rs = [0.92_pr, 2.1055_pr, 3.1878_pr]
      !! qs = [1.4_pr, 1.972_pr, 2.4_pr]
      !!
      !! T = 298.15_pr
      !!
      !! ! Calculate bij from DUij. We need -DU/R to get bij
      !! b(1, :) = [0.0_pr, -526.02_pr, -309.64_pr]
      !! b(2, :) = [318.06_pr, 0.0_pr, 91.532_pr]
      !! b(3, :) = [-1325.1_pr, -302.57_pr, 0.0_pr]
      !!
      !! model = setup_uniquac(qs, rs, bij=b)
      !!
      !! n = [0.8_pr, 0.1_pr, 0.2_pr]
      !!
      !! call model%ln_activity_coefficient(n, T, ln_gammas)
      !!
      !! print *, exp(ln_gammas) ! [8.856, 0.860, 1.425]
      !!
      !! ```
      !!
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
      !! Calculate the excess Gibbs free energy and its derivatives of the
      !! UNIQUAC model.
      !!
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

      ! Main terms
      real(pr) :: thetak(size(n))
      real(pr) :: dthetak_dni(size(n), size(n))
      real(pr) :: d2thetak_dnidnj(size(n), size(n), size(n))

      real(pr) :: tau(size(n), size(n))
      real(pr) :: dtau(size(n), size(n))
      real(pr) :: d2tau(size(n), size(n))

      ! Indexes and lofical for optional arguments
      integer :: i, j, k

      logical :: dt, dt2, dn, dtn, dn2, dt_or_dtn

      ! Auxiliars
      integer :: nc
      real(pr) :: n_tot
      real(pr) :: r_i, q_i, r_j, q_j, r_k, q_k
      real(pr) :: sum_nq, sum_nr
      real(pr) :: Ge_comb, Ge_res
      real(pr) :: Ge_aux

      ! Auxiliars for temperature derivatives
      real(pr) :: sum_thetal_tau_lk(size(n))
      real(pr) :: sum_theta_l_dtau_lk(size(n))
      real(pr) :: sum_theta_l_d2tau_lk(size(n))

      real(pr) :: GeT_aux, GeT2_aux, diff_aux(size(n))

      ! Auxiliars for compositionals derivatives
      real(pr) :: Gen_comb(size(n)), Gen_res(size(n))
      real(pr) :: Gen_aux(size(n))

      real(pr) :: dln_nk_dni(size(n), size(n))

      ! Auxiliars for second compositionals derivatives
      ! terms of the math expression
      real(pr) :: trm1, trm2, trm3, trm4, trm5, trm6

      real(pr) :: sum_d2theta_tau_lk(size(n))
      real(pr) :: Gen2_aux(size(n), size(n))

      ! cross derivative auxiliars
      real(pr) :: sum_dtheta_l_dtau_lk(size(n), size(n))
      real(pr) :: sum_dtheta_l_tau_lk(size(n), size(n))
      real(pr) :: GeTn_aux(size(n))

      ! =======================================================================
      ! Logical variables
      ! -----------------------------------------------------------------------
      dt = present(GeT)
      dt2 = present(GeT2)
      dn = present(Gen)
      dtn = present(GeTn)
      dn2 = present(Gen2)

      dt_or_dtn = dt .or. dtn

      ! =======================================================================
      ! Auxiliars
      ! -----------------------------------------------------------------------
      nc = size(n)

      n_tot = sum(n)

      sum_nq = sum(n * self%qs)
      sum_nr = sum(n * self%rs)
      ! =======================================================================
      ! tau call (temperature dependence term)
      ! -----------------------------------------------------------------------
      if (dt_or_dtn .and. .not. dt2) call self%taus(T, tau, dtau)
      if (.not. dt_or_dtn .and. dt2) call self%taus(T, tau, dtau, d2tau)
      if (dt_or_dtn .and. dt2) call self%taus(T, tau, dtau, d2tau)
      if (.not. dt_or_dtn .and. .not. dt2) call self%taus(T, tau)

      ! =======================================================================
      ! theta_k
      ! -----------------------------------------------------------------------
      thetak = n * self%qs / sum_nq

      if (dn .or. dn2 .or. dtn) then
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
      ! Ge
      ! -----------------------------------------------------------------------
      ! Combinatorial term
      Ge_comb = ( &
         sum(n * log(n_tot * self%rs / sum_nr)) &
         + self%z / 2.0_pr * sum(n * self%qs * log(self%qs * sum_nr / self%rs / sum_nq)) &
         )

      ! Residual term
      sum_thetal_tau_lk = 0.0_pr

      do k=1,nc
         sum_thetal_tau_lk(k) = sum(thetak * tau(:,k))
      end do

      Ge_res = -sum(n * self%qs * log(sum_thetal_tau_lk))

      Ge_aux = R * T * (Ge_comb + Ge_res)

      ! =======================================================================
      ! Ge Derivatives
      ! -----------------------------------------------------------------------
      if (dn .or. dtn .or. dn2) then
         do concurrent (k=1:nc, i=1:nc)
            sum_dtheta_l_tau_lk(i,k) = sum(dthetak_dni(:,i) * tau(:,k))
         end do
      end if

      ! Compositional derivarives
      ! dn
      if (dn .or. dtn) then
         do i=1,nc
            ! Combinatorial term
            Gen_comb(i) = (&
               log(n_tot * self%rs(i) / sum_nr) &
               + sum(n * (sum_nr - n_tot * self%rs(i)) / n_tot / sum_nr) &
               + self%z / 2.0_pr * self%qs(i) * log(self%qs(i) * sum_nr / self%rs(i) / sum_nq) &
               + self%z / 2.0_pr * sum(n * self%qs * (self%rs(i) * sum_nq - self%qs(i) * sum_nr) / sum_nq / sum_nr) &
               )

            ! Residual term
            Gen_res(i) = (&
               -self%qs(i) * log(sum_thetal_tau_lk(i)) &
               -sum(n * self%qs * sum_dtheta_l_tau_lk(i,:) / sum_thetal_tau_lk) &
               )
         end do

         Gen_aux = R * T * (Gen_comb + Gen_res)
      end if

      ! dn2
      if (dn2) then
         do concurrent(i=1:nc, j=1:nc)
            trm1 = (sum_nr - n_tot * self%rs(j)) / n_tot / sum_nr

            trm2 = (&
               (sum_nr - n_tot * self%rs(i)) / n_tot / sum_nr &
               + sum(n * (self%rs(i) * self%rs(j) / sum_nr**2 - 1.0_pr / n_tot**2)) &
               )

            trm3 = (&
               self%z / 2.0_pr * self%qs(i) * (self%rs(j) * sum_nq - &
               self%qs(j) * sum_nr) / sum_nq / sum_nr &
               )

            trm4 = (&
               self%z / 2.0_pr * self%qs(j) * (self%rs(i) * sum_nq - &
               self%qs(i) * sum_nr) / sum_nq / sum_nr &
               + self%z / 2.0_pr * sum(n * self%qs * &
               (self%qs(i) * self%qs(j) / sum_nq**2 - &
               self%rs(i) * self%rs(j) / sum_nr**2))&
               )

            trm5 = -self%qs(i) * (sum(dthetak_dni(:,j) * tau(:,i)) / sum_thetal_tau_lk(i))

            do k=1,nc
               sum_d2theta_tau_lk(k) = sum(d2thetak_dnidnj(:,i,j) * tau(:,k))
            end do

            trm6 = (&
               -self%qs(j) * (sum(dthetak_dni(:,i) * tau(:,j)) / sum_thetal_tau_lk(j)) &
               -sum(self%qs * n * (&
               sum_d2theta_tau_lk * sum_thetal_tau_lk &
               - sum_dtheta_l_tau_lk(i,:) * sum_dtheta_l_tau_lk(j,:) &
               ) / sum_thetal_tau_lk**2) &
               )

            Gen2_aux(i,j) = R * T * (trm1 + trm2 + trm3 + trm4 + trm5 + trm6)
         end do
      end if

      ! Temperature derivatives
      if (dt .or. dt2 .or. dtn) then
         sum_theta_l_dtau_lk = 0.0_pr

         do k=1,nc
            sum_theta_l_dtau_lk(k) = sum(thetak * dtau(:,k))
         end do
      end if

      if (dt) then
         GeT_aux = ( &
            Ge_aux/T - R*T*sum(n * self%qs * sum_theta_l_dtau_lk / sum_thetal_tau_lk)&
            )
      end if

      if (dt2) then
         sum_theta_l_d2tau_lk = 0.0_pr

         do k=1,nc
            sum_theta_l_d2tau_lk(k) = sum(thetak * d2tau(:,k))
         end do

         diff_aux = (&
            sum_theta_l_d2tau_lk / sum_thetal_tau_lk &
            - (sum_theta_l_dtau_lk / sum_thetal_tau_lk)**2 &
            )

         GeT2_aux = -R * ( &
            T * sum(self%qs * n * diff_aux) &
            + 2.0_pr*sum(self%qs * n * sum_theta_l_dtau_lk / sum_thetal_tau_lk)&
            )
      end if

      ! Cross derivative Tn
      if (dtn) then
         do concurrent (k=1:nc, i=1:nc)
            sum_dtheta_l_dtau_lk(i,k) = sum(dthetak_dni(:,i) * dtau(:,k))
         end do
      end if

      if (dtn) then
         do i=1,nc
            GeTn_aux(i) = ( &
               1.0_pr / T  * Gen_aux(i) &
               -R * T * (&
               self%qs(i) * sum_theta_l_dtau_lk(i) / sum_thetal_tau_lk(i) &
               + sum(n * self%qs * (&
               sum_dtheta_l_dtau_lk(i,:) * sum_thetal_tau_lk &
               - sum_theta_l_dtau_lk * sum_dtheta_l_tau_lk(i,:)) &
               / sum_thetal_tau_lk**2) &
               ) &
               )
         end do
      end if
      ! =======================================================================
      ! Excess Gibbs energy returns
      ! -----------------------------------------------------------------------
      if (present(Ge)) Ge = Ge_aux
      if (dt) GeT = GeT_aux
      if (dt2) GeT2 = GeT2_aux
      if (dn) Gen = Gen_aux
      if (dtn) GeTn = GeTn_aux
      if (dn2) Gen2 = Gen2_aux
   end subroutine excess_gibbs

   subroutine taus(self, T, tau, tauT, tauT2)
      !! Calculate the temperature dependence term of the UNIQUAC model.
      !!
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
         du = -self%bij / T**2 + self%cij / T + self%dij + 2.0_pr * T * self%eij
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
      !! Instantiate a UNIQUAC model.
      !!
      !! Non provided interaction parameters are set to zero matrices.
      !!
      real(pr), intent(in) :: qs(:)
      !! Molecule's relative volumes \(Q_i\)
      real(pr), intent(in) :: rs(size(qs))
      !! Molecule's relative areas \(R_i\)
      real(pr), optional, intent(in) :: aij(size(qs),size(qs))
      !! Interaction parameters matrix \(a_{ij}\), zero matrix if no provided.
      real(pr), optional, intent(in) :: bij(size(qs),size(qs))
      !! Interaction parameters matrix \(b_{ij}\), zero matrix if no provided.
      real(pr), optional, intent(in) :: cij(size(qs),size(qs))
      !! Interaction parameters matrix \(c_{ij}\), zero matrix if no provided.
      real(pr), optional, intent(in) :: dij(size(qs),size(qs))
      !! Interaction parameters matrix \(d_{ij}\), zero matrix if no provided.
      real(pr), optional, intent(in) :: eij(size(qs),size(qs))
      !! Interaction parameters matrix \(e_{ij}\), zero matrix if no provided.

      ! aij
      if (present(aij)) then
         setup_uniquac%aij = aij
      else
         allocate(setup_uniquac%aij(size(qs),size(qs)))
         setup_uniquac%aij = 0.0_pr
      end if

      ! bij
      if (present(bij)) then
         setup_uniquac%bij = bij
      else
         allocate(setup_uniquac%bij(size(qs),size(qs)))
         setup_uniquac%bij = 0.0_pr
      end if

      ! cij
      if (present(cij)) then
         setup_uniquac%cij = cij
      else
         allocate(setup_uniquac%cij(size(qs),size(qs)))
         setup_uniquac%cij = 0.0_pr
      end if

      ! dij
      if (present(dij)) then
         setup_uniquac%dij = dij
      else
         allocate(setup_uniquac%dij(size(qs),size(qs)))
         setup_uniquac%dij = 0.0_pr
      end if

      ! eij
      if (present(eij)) then
         setup_uniquac%eij = eij
      else
         allocate(setup_uniquac%eij(size(qs),size(qs)))
         setup_uniquac%eij = 0.0_pr
      end if

      setup_uniquac%qs = qs
      setup_uniquac%rs = rs
   end function setup_uniquac
end module yaeos__models_ge_uniquac
