module yaeos__models_ar_saft_pcsaft
   !! PC-SAFT Implementation (Gross & Sadowski, 2001)
   !! Approach: Hard Chain + Dispersion (Placeholder)

   use yaeos__constants, only: pr, R
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   use hyperdual_mod

   implicit none

   private

   public :: PcSaft, init_pcsaft

   ! =========================================================================
   ! PC-SAFT UNIVERSAL CONSTANTS (Gross & Sadowski, 2001, Table A1)
   ! =========================================================================
   ! A_COEFFS(k, i): k = power of eta (0..6), i = type (0:m=1, 1:m=inf, 2:corr)
   real(pr), dimension(0:2, 0:6) :: A_COEFFS = reshape([ &
      0.9105631445_pr, -0.3084016918_pr, -0.0906148351_pr, &
      0.6361281449_pr, 0.1860531159_pr, 0.4527842806_pr,  &
      2.6861347891_pr, -2.5030047259_pr, 0.5962700728_pr, &
      -26.547362491_pr, 21.419793629_pr, -1.7241829131_pr,&
      97.759208784_pr, -65.255885330_pr, -4.1302112531_pr,&
      -159.59154087_pr, 83.318680481_pr, 13.776631870_pr, &
      91.297774084_pr, -33.746922930_pr, -8.6728470368_pr ], [3, 7])

   real(pr), dimension(0:2, 0:6) :: B_COEFFS = reshape([ &
      0.7240946941_pr, -0.5755498075_pr, 0.0976883116_pr, &
      2.2382791861_pr, 0.6995095521_pr, -0.2557574982_pr, &
      -4.0025849485_pr, 3.8925673390_pr, -9.1558561530_pr, &
      -21.003576815_pr, -17.215471648_pr, 20.642075974_pr, &
      26.855641363_pr, 192.67226447_pr, -38.804430052_pr, &
      206.55133841_pr, -161.82646165_pr, 93.626774077_pr, &
      -355.60235612_pr, -165.20769346_pr, -29.666905585_pr], [3, 7])


   type, extends(ArModelAdiff) :: PcSaft
      !! # `PcSaft`
      !! PC-SAFT Equation of State Model
      real(pr), allocatable :: m(:)         !! Number of segments
      real(pr), allocatable :: sigma(:)     !! Segment diameter [Angstrom]
      real(pr), allocatable :: epsilon_k(:) !! Energy / k_B [K]
      real(pr), allocatable :: kij(:,:)     !! Binary interaction parameters (optional)
      real(pr), allocatable :: eps_ab(:)    !! Association energy [K]
      real(pr), allocatable :: kap_ab(:)    !! Association volume [A^3]
      integer, allocatable :: na_sites(:)  !! Number of aceptors sites
      integer, allocatable :: nb_sites(:)  !! Number of donor sites
   contains
      procedure :: Ar => Ar_impl
      procedure :: get_v0 => get_v0_impl
   end type PcSaft

   ! Module private constants
   real(pr), parameter :: PI = 3.14159265359_pr
   real(pr), parameter :: N_AVO = 6.02214076e23_pr

   ! Critical conversion factor:
   ! Zeta must be dimensionless.
   ! V comes in Liters. Sigma in Angstroms.
   ! 1 L = 10^27 A^3.
   ! rho_num [1/A^3] = (n [mol] * N_AVO) / (V [L] * 10^27)
   ! Factor = 0.602214...
   real(pr), parameter :: UNITS_FACTOR = 0.000602214086_pr

contains

   type(PcSaft) function init_pcsaft(m, sigma, epsilon_k, kij, eps_ab, kap_ab, na_sites, nb_sites) result(model)
      use yaeos__equilibria_critical, only: get_critical_constants
      real(pr), intent(in) :: m(:)
      real(pr), intent(in) :: sigma(:)
      real(pr), intent(in) :: epsilon_k(:)
      real(pr), intent(in), optional :: kij(:,:)
      real(pr), intent(in), optional :: eps_ab(:)
      real(pr), intent(in), optional :: kap_ab(:)
      integer, intent(in), optional :: na_sites(:)
      integer, intent(in), optional :: nb_sites(:)

      integer :: nc
      model%m = m
      model%sigma = sigma
      model%epsilon_k = epsilon_k
      if (present(kij)) then
         model%kij = kij
      end if

      if (present(eps_ab) .and. present(kap_ab) .and. present(na_sites) .and. present(nb_sites)) then
         model%eps_ab = eps_ab
         model%kap_ab = kap_ab
         model%na_sites = na_sites
         model%nb_sites = nb_sites
      end if

      !nc = size(m)
      !allocate(model%components%Tc(nc))
      !allocate(model%components%Pc(nc))
      !allocate(model%components%w(nc))
      !call get_critical_constants(model)
   end function init_pcsaft


   ! ====================================================================
   ! Main Function: Residual Helmholtz Energy
   ! ====================================================================
   function Ar_impl(self, n, V, T) result(ar_total)
      class(PcSaft) :: self
      type(hyperdual), intent(in) :: n(:)  ! Moles [mol]
      type(hyperdual), intent(in) :: V     ! Volume [L]
      type(hyperdual), intent(in) :: T     ! Temperature [K]
      type(hyperdual) :: ar_total          ! A_res Total [bar L]

      ! Local variables
      type(hyperdual) :: d(size(n))       ! T-dependent diameter
      type(hyperdual) :: zeta(0:3)        ! Density moments
      type(hyperdual) :: a_hs             ! A_hard_sphere / RT (dimensionless)
      type(hyperdual) :: a_chain          ! A_chain / RT (dimensionless)
      type(hyperdual) :: a_disp           ! A_disp / RT (dimensionless)
      type(hyperdual) :: a_assoc          ! A_assoc / RT (dimensionless)
      type(hyperdual) :: n_tot
      type(hyperdual) :: rho, eta, m_ave, x(size(n))

      integer :: nc, i

      nc = size(n)
      n_tot = sum(n)
      x = n / n_tot
      m_ave = sum(x * self%m)

      ! 1. Calculate Segment Diameter d(T) [Eq. A.4]
      ! IMPORTANT: 'd' is hyperdual because it depends on T.
      do i = 1, nc
         d(i) = self%sigma(i) * (1.0_pr - 0.12_pr * exp(-3.0_pr * self%epsilon_k(i) / T))
      end do

      ! 2. Calculate Density Moments (Zetas) [Eq. A.5]
      ! Passing n, V and d.
      call calculate_zetas(n, V, self%m, d, zeta)

      ! 3. Calculate Hard Sphere Contribution (Mixture) [Eq. A.6 - A.7]
      a_hs = calculate_hard_sphere(zeta, n_tot)

      ! 4. Calculate Chain Contribution [Eq. A.8 - A.10]
      ! Requires the radial distribution function (g_hs)
      a_chain = calculate_chain(x, n_tot, self%m, d, zeta)

      ! 5. Calculate Dispersion Contribution
      ! We pass self%kij if it exists, otherwise optional
      a_disp = 0.0_pr
      if (allocated(self%kij)) then
         a_disp = calculate_dispersion(n, V, T, zeta, self%m, self%epsilon_k, self%sigma, self%kij)
      else
         a_disp = calculate_dispersion(n, V, T, zeta, self%m, self%epsilon_k, self%sigma)
      end if

      ! 6. Calculate Association (Optional)
      if (allocated(self%eps_ab) .and. allocated(self%kap_ab) .and. allocated(self%na_sites)) then
         a_assoc = calculate_association(n, V, T, zeta, d, self%sigma, self%eps_ab, self%kap_ab, self%na_sites, self%nb_sites)
      else
         a_assoc = 0.0_pr
      end if

      print *, "Ar_hs: ", R * T%f0 * m_ave%f0 * a_hs%f0 * 100
      print *, "Ar_chain: ", R * T%f0 * a_chain%f0 * 100 / n_tot%f0
      print *, "Ar_disp: ", R * T%f0 * a_disp%f0 * 100 / n_tot%f0
      print *, "Ar_assoc: ", R * T%f0 * a_assoc%f0 * 100 / 1000 / n_tot%f0

      ! 7. Sum and convert to Energy units [bar * L]
      ar_total = (R * T) * ( n_tot * m_ave * a_hs + a_chain + a_disp + a_assoc)
   end function Ar_impl

   ! ====================================================================
   ! Auxiliary Routines (Hard Chain)
   ! ====================================================================
   subroutine calculate_zetas(n, V, m, d, zeta)
      type(hyperdual), intent(in) :: n(:), V, d(:)
      real(pr), intent(in) :: m(:)
      type(hyperdual), intent(out) :: zeta(0:3)

      integer :: k, i
      type(hyperdual) :: term
      type(hyperdual) :: rho
      type(hyperdual) :: x(size(n))

      rho = sum(n) * N_AVO / (V * 1.0e27_pr)  ! Number density [1/A^3]
      x = n / sum(n)

      zeta = 0.0_pr

      do k = 0, 3
         zeta(k) = rho * (PI / 6.0_pr) * sum(x * m * (d**k))
      end do
   end subroutine calculate_zetas

   function calculate_hard_sphere(zeta, n_tot) result(val)
      type(hyperdual), intent(in) :: zeta(0:3), n_tot
      type(hyperdual) :: val

      ! Boublik-Mansoori-Carnahan-Starling Equation (Eq. A.6 and A.7)
      ! A_hs/RT = (1/zeta0) * [ ... ]
      ! But careful: Gross-Sadowski Eq A.6 is on a molar basis.
      ! A_hs_total/RT = n_total * a_hs_molar

      ! We implement the expanded form to avoid division by zero if zeta0 is very small,
      ! although in PC-SAFT zeta0 is never zero if there is mass.

      type(hyperdual) :: term1, term2, term3, one_m_z3

      one_m_z3 = 1.0_pr - zeta(3)

      term1 = (3.0_pr * zeta(1) * zeta(2)) / one_m_z3
      term2 = (zeta(2)**3) / (zeta(3) * (one_m_z3**2))
      term3 = ((zeta(2)**3)/(zeta(3)**2) - zeta(0)) * log(one_m_z3)

      ! Correct formula for TOTAL Free Energy A_hs/RT:
      val = (1.0_pr / zeta(0)) * (term1 + term2 + term3)
   end function calculate_hard_sphere

   function calculate_chain(x, n_tot, m, d, zeta) result(val)
      type(hyperdual), intent(in) :: x(:), n_tot, d(:), zeta(0:3)
      real(pr), intent(in) :: m(:)
      type(hyperdual) :: val

      ! Eq. A.8: A_chain = - sum( x_i * (m_i - 1) * ln(g_ii_hs) )
      ! In total units: A_chain_total/RT = - sum( n_i * (m_i - 1) * ln(g_ii_hs) )

      integer :: i
      type(hyperdual) :: g_ii, one_m_z3, z2_term, z2_term2
      type(hyperdual) :: di_2

      val = 0.0_pr
      one_m_z3 = 1.0_pr - zeta(3)

      ! Pre-calculate common RDF terms (Eq. A.9)
      ! g_ij = 1/(1-z3) + (di*dj/(di+dj)) * 3z2/(1-z3)^2 + ...
      ! For Chain we only need g_ii (di=dj), so di*dj/(di+dj) = di/2

      do i = 1, size(x)
         ! Contact RDF g_ii (Eq. A.9 simplified for i=j)
         ! di*di / (di+di) = di / 2

         di_2 = d(i) / 2.0_pr

         g_ii = (1.0_pr / one_m_z3) + &
            (3.0_pr * di_2 * zeta(2)) / (one_m_z3**2) + &
            (2.0_pr * (di_2**2) * (zeta(2)**2)) / (one_m_z3**3)

         ! Summation
         val = val - x(i) * (m(i) - 1.0_pr) * log(g_ii)
      end do

      val = val * n_tot

   end function calculate_chain

   function calculate_dispersion(n, V, T, zeta, m, eps_k, sig, kij) result(val)
      type(hyperdual), intent(in) :: n(:), V, T, zeta(0:3)
      real(pr), intent(in) :: m(:), eps_k(:), sig(:)
      real(pr), intent(in), optional :: kij(:,:)
      type(hyperdual) :: val

      ! Local variables
      integer :: i, j, k, nc
      type(hyperdual) :: n_tot, rho, eta, one_m_eta
      type(hyperdual) :: m_ave, m2_es3, m2_e2s3, I1, I2, C1, x(size(n))
      type(hyperdual) :: term, a_k, b_k
      type(hyperdual) :: a1_term, a2_term ! <--- Separate the terms
      real(pr) :: eps_ij, sig_ij, kij_val

      nc = size(n)
      n_tot = sum(n)
      rho = n_tot * N_AVO / (V * 1.0e27_pr)  ! Number density [1/A^3]
      eta = zeta(3)
      one_m_eta = 1.0_pr - eta
      x = n / n_tot

      ! 1. Mixture Averages (Same as before)
      m_ave = 0.0_pr; m2_es3 = 0.0_pr; m2_e2s3 = 0.0_pr

      do i = 1, nc
         m_ave = m_ave + x(i)*m(i) ! Segment average
      end do

      do i = 1, nc
         do j = 1, nc
            sig_ij = 0.5_pr * (sig(i) + sig(j))
            kij_val = 0.0_pr; if (present(kij)) kij_val = kij(i,j)
            eps_ij = sqrt(eps_k(i) * eps_k(j)) * (1.0_pr - kij_val)

            term = x(i) * x(j) * m(i) * m(j) * (sig_ij**3)

            m2_es3 = m2_es3 + term * (eps_ij / T)
            m2_e2s3 = m2_e2s3 + term * (eps_ij / T)**2
         end do
      end do

      ! 2. Integrals I1 and I2 (Same as before)
      I1 = 0.0_pr; I2 = 0.0_pr
      do k = 0, 6
         a_k = A_COEFFS(0,k) + (m_ave - 1.0_pr)/m_ave * A_COEFFS(1,k) + &
            (m_ave - 1.0_pr)/m_ave * (m_ave - 2.0_pr)/m_ave * A_COEFFS(2,k)
         b_k = B_COEFFS(0,k) + (m_ave - 1.0_pr)/m_ave * B_COEFFS(1,k) + &
            (m_ave - 1.0_pr)/m_ave * (m_ave - 2.0_pr)/m_ave * B_COEFFS(2,k)

         I1 = I1 + a_k * (eta**k)
         I2 = I2 + b_k * (eta**k)
      end do

      ! 3. Compressibility C1 (Same as before)
      C1 = 1.0_pr + m_ave * (8.0_pr * eta - 2.0_pr * (eta**2)) / (one_m_eta**4)
      C1 = C1 + (1.0_pr - m_ave) * &
         (20.0_pr*eta - 27.0_pr*(eta**2) + 12.0_pr*(eta**3) - 2.0_pr*(eta**4)) / &
         ((1.0_pr - eta) * (2.0_pr - eta))**2
      C1 = 1._pr / C1

      ! ---------------------------------------------------------
      ! 4. Final Sum
      ! Eq. A.11 Gross & Sadowski (2001)
      ! a_res/RT = -2*pi*rho*I1*... - pi*rho*m_ave*C1*I2*...
      ! ---------------------------------------------------------

      ! 1st order term (Takes 2 * PI)
      a1_term = -2.0_pr * PI * rho * I1 * m2_es3

      ! 2nd order term (Takes 1 * PI and multiplies by m_ave)
      a2_term = -1.0_pr * PI * rho * m_ave * C1 * I2 * m2_e2s3

      ! Sum and make extensive
      val = (a1_term + a2_term) * n_tot

   end function calculate_dispersion

   function calculate_association(n, V, T, zeta, d, sigmas, eps_ab, kap_ab, na_sites, nb_sites) result(val)
      type(hyperdual), intent(in) :: n(:), V, T, zeta(0:3), d(:)
      real(pr), intent(in)        :: sigmas(:), eps_ab(:), kap_ab(:) ! Parametros puros (eps/kb en [K], kappa es adimensional)
      integer, intent(in)         :: na_sites(:), nb_sites(:)
      type(hyperdual)             :: val

      integer :: i, j, iter, n_comp
      type(hyperdual)             :: rho_num, total_n
      type(hyperdual), allocatable :: XA(:), XB(:), XA_old(:), XB_old(:)
      type(hyperdual), allocatable :: delta(:,:)
      type(hyperdual)             :: sum_A, sum_B, g_ij, d_ij, eps_mix, kap_mix, di_dj_term
      real(pr)                    :: error_f0
      
      real(pr), parameter :: ALPHA = 0.5_pr
      n_comp = size(n)

      allocate(XA(n_comp), XB(n_comp), XA_old(n_comp), XB_old(n_comp))
      allocate(delta(n_comp, n_comp))

      ! 1. Densidad numérica total [1/Angstrom^3]
      total_n = sum(n)
      rho_num = total_n * N_AVO / (V * 1.0e27_pr)

      ! 2. Cálculo de la Matriz Delta (Interacción Aceptor_i -> Donador_j)
      do i = 1, n_comp
         do j = 1, n_comp
            
            ! Si alguno de los dos componentes no tiene capacidad de asociarse, Delta = 0
            if ((na_sites(i) == 0 .and. nb_sites(i) == 0) .or. &
               (na_sites(j) == 0 .and. nb_sites(j) == 0)) then
               delta(i,j) = 0.0_pr
               cycle
            end if

            ! RDF de Esferas Duras a contacto (g_hs)
            d_ij = 0.5_pr * (d(i) + d(j))

            ! g_ij = 1.0_pr / (1.0_pr - zeta(3)) + &
            !       (3.0_pr * d(i) * d(j) / (d(i) + d(j)) * zeta(2)) / ((1.0_pr - zeta(3))**2) + &
            !       (2.0_pr * (d(i) * d(j) / (d(i) + d(j)))**2 * zeta(2)**2) / ((1.0_pr - zeta(3))**3)
            
            g_ij = 1.0_pr / (1.0_pr - zeta(3)) + &
                  (3.0_pr * d_ij * zeta(2)) / ((1.0_pr - zeta(3))**2) + &
                  (2.0_pr * (d_ij)**2 * zeta(2)**2) / ((1.0_pr - zeta(3))**3)

            ! Regla de mezclado Geométrica de Wolbach-Sandler para Energía
            eps_mix = 0.5_pr * (eps_ab(i) + eps_ab(j))

            ! Corrección geométrica para el Volumen de Asociación (Kappa)
            kap_mix = sqrt(kap_ab(i) * kap_ab(j)) * (sqrt(sigmas(i) * sigmas(j)) &
            / (0.5_pr * (sigmas(i) + sigmas(j))))**3

            ! Cálculo de Delta_ij [Angstroms^3]
            delta(i,j) = (d_ij**3) * g_ij * kap_mix * (exp(eps_mix / T) - 1.0_pr)
         end do
      end do

      ! 3. Solución acoplada de XA y XB (Sustitución Sucesiva con Amortiguación)
      ! Inicialización con un guess hiperdual limpio (derivadas en 0, real en 0.2)
      XA = 0.2_pr
      XB = 0.2_pr

      do iter = 1, 300
         XA_old = XA
         XB_old = XB

         do i = 1, n_comp
            if (na_sites(i) == 0 .and. nb_sites(i) == 0) cycle

            sum_A = 0.0_pr
            sum_B = 0.0_pr

            do j = 1, n_comp
               if (na_sites(j) == 0 .and. nb_sites(j) == 0) cycle
               
               ! Recordar: Fracción molar x_j = n(j)/total_n
               ! Para XA(i) necesito los sitios Donadores (nb) de j y sus XB(j)
               sum_A = sum_A + (n(j) / total_n) * rho_num * real(nb_sites(j), pr) * XB_old(j) * delta(i,j)
               
               ! Para XB(i) necesito los sitios Aceptores (na) de j y sus XA(j)
               ! Nota: delta(j,i) representa Aceptor_j interaccionando con Donador_i
               sum_B = sum_B + (n(j) / total_n) * rho_num * real(na_sites(j), pr) * XA_old(j) * delta(j,i)
            end do

            XA(i) = 1.0_pr / (1.0_pr + sum_A)
            XB(i) = 1.0_pr / (1.0_pr + sum_B)
         end do

         ! Criterio de convergencia basado ESTRICTAMENTE en la parte real %f0
         error_f0 = 0.0_pr
         do i = 1, n_comp
            error_f0 = max(error_f0, abs(XA(i)%f0 - XA_old(i)%f0), abs(XB(i)%f0 - XB_old(i)%f0))
         end do

         ! if (error_f0 < 1.0e-12_pr) exit

         ! Aplicar factor de amortiguación sobre todo el objeto hyperdual (mantiene álgebra de derivadas)
         XA = ALPHA * XA + (1.0_pr - ALPHA) * XA_old
         XB = ALPHA * XB + (1.0_pr - ALPHA) * XB_old
      end do

      ! 4. Cálculo del potencial de Helmholtz Residual Total (A_assoc / RT)
      val = 0.0_pr
      do i = 1, n_comp
         if (na_sites(i) == 0 .and. nb_sites(i) == 0) cycle
         
         val = val + n(i) * ( &
               real(na_sites(i), pr) * (log(XA(i)) - 0.5_pr * XA(i) + 0.5_pr) + &
               real(nb_sites(i), pr) * (log(XB(i)) - 0.5_pr * XB(i) + 0.5_pr) )
      end do

      deallocate(XA, XB, XA_old, XB_old, delta)

   end function calculate_association

   ! ---------------------------------------------------------
   ! Method get_v0: Volume lower limit (Covolume)
   ! Represents the physical volume occupied by the segments
   ! ---------------------------------------------------------
   function get_v0_impl(self, n, P, T) result(v0)
      class(PcSaft), intent(in) :: self
      real(pr), intent(in) :: n(:)  ! Moles of each component
      real(pr), intent(in) :: P     ! System pressure
      real(pr), intent(in) :: T     ! System temperature
      real(pr) :: v0

      integer :: i
      real(pr) :: sum_seg_vol

      ! v0 = (Pi/6) * Sum(n_i * m_i * sigma_i^3) * Factor_Conversion
      ! This makes zeta_3 = 1 if V = v0.

      sum_seg_vol = 0.0_pr
      do i = 1, size(n)
         sum_seg_vol = sum_seg_vol + n(i) * self%m(i) * (self%sigma(i)**3)
      end do

      ! Note: UNITS_FACTOR must be approx 0.000602214 to convert
      ! (mol * Ang^3) to Liters.
      ! Make sure it matches the one used in 'calculate_zetas'.
      v0 = (PI / 6.0_pr) * UNITS_FACTOR * sum_seg_vol

   end function get_v0_impl

end module yaeos__models_ar_saft_pcsaft
