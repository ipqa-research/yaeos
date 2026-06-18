module tapenade_model_template
   use yaeos__tapenade_ar_api, only: ArModelTapenade
   use yaeos__tapenade_interfaces
   use yaeos, only: R
   implicit none

   ! ---------------------------------------------------------------------------
   ! PC-SAFT UNIVERSAL CONSTANTS (Gross & Sadowski, 2001, Table A1)
   ! ---------------------------------------------------------------------------
   real(8), dimension(0:2, 0:6) :: A_COEFFS = reshape([ &
      0.9105631445, -0.3084016918, -0.0906148351, &
      0.6361281449 , 0.1860531159 , 0.4527842806,  &
      2.6861347891 , -2.5030047259, 0.5962700728, &
      -26.547362491, 21.419793629 ,-1.7241829131,&
      97.759208784 ,-65.255885330 ,-4.1302112531,&
      -159.59154087, 83.318680481 ,13.776631870, &
      91.297774084 ,-33.746922930 ,-8.6728470368 ], [3, 7])

   real(8), dimension(0:2, 0:6) :: B_COEFFS = reshape([ &
      0.7240946941, -0.5755498075, 0.0976883116, &
      2.2382791861, 0.6995095521, -0.2557574982, &
      -4.0025849485, 3.8925673390, -9.1558561530, &
      -21.003576815, -17.215471648, 20.642075974, &
      26.855641363, 192.67226447, -38.804430052, &
      206.55133841, -161.82646165, 93.626774077, &
      -355.60235612, -165.20769346, -29.666905585], [3, 7])

   !type, extends(ArModelTapenade) :: PCSAFT
   type :: PCSAFT
      real(8), allocatable :: m(:)
      !! Number of segments
      real(8), allocatable :: sigma(:)
      !! Segment diameter [Angstrom]
      real(8), allocatable :: epsilon_k(:)
      !! Energy / k_B [K]
      real(8), allocatable :: kij(:,:)
      !! Binary interaction parameters (optional)
      ! contains
      !  procedure :: ar
      !  procedure :: ar_d
      !  procedure :: ar_b
      !  procedure :: ar_d_b
      !  procedure :: ar_d_d
      !  procedure :: v0
   end type PCSAFT

   ! Module private constants
   real(8), parameter :: PI = 3.14159265359_8
   real(8), parameter :: N_AVO = 6.02214076e23_8

   ! Conversion factor:
   ! Zeta must be dimensionless. V comes in Liters. Sigma in Angstroms.
   ! 1 L = 10^27 A^3. rho_num [1/A^3] = (n [mol] * N_AVO) / (V [L] * 10^27)
   real(8), parameter :: UNITS_FACTOR = 0.000602214086_8

contains
   type(PCSAFT) function setup_model(m, sigma, epsilon_k, kij) result(model)
      use yaeos__equilibria_critical, only: get_critical_constants
      real(8), intent(in) :: m(:)
      !! Number of segments
      real(8), intent(in) :: sigma(:)
      !! Segment diameter [Angstrom]
      real(8), intent(in) :: epsilon_k(:)
      !! Energy / k_B [K]
      real(8), intent(in), optional :: kij(:,:)
      !! Binary interaction parameters (optional)

      integer :: nc
      nc = size(m)

      model%m = m
      model%sigma = sigma
      model%epsilon_k = epsilon_k
      
      if (present(kij)) then
         model%kij = kij
      end if

      allocate(model%components%Tc(nc))
      allocate(model%components%Pc(nc))
      allocate(model%components%w(nc))
      call get_critical_constants(model)
   end function setup_model

   subroutine ar(model, n, v, t, arval)
      type(PCSAFT), intent(in) :: model
      real(8), intent(in) :: n(:), v, t
      real(8), intent(out) :: arval

      ! Auxiliars
      integer :: i, j, k
      real(8) :: x(size(n)) ! Mole fractions
      real(8) :: n_tot ! Total moles
      integer :: nc ! Number of components
      real(8) :: m_average ! Average number of segments

      ! Internal variables declarations
      ! Diameter declarations
      real(8) :: d(size(n))
      ! Zetas declarations
      real(8) :: zetas(0:3)
      real(8) :: rho
      ! Ar hard-sphere contribution declarations
      real(8) :: term1_hs, term2_hs, term3_hs, one_m_z3_hs
      real(8) :: ar_hs
      ! Ar chain contribution declarations
      real(8) :: g_ii, one_m_z3_chain
      real(8) :: di_2
      real(8) :: ar_chain
      ! Ar dispersion contribution declarations
      real(8) :: rho_disp, eta, one_m_eta
      real(8) :: m_ave, m2_es3, m2_e2s3, I1, I2, C1
      real(8) :: term_disp, a_k, b_k
      real(8) :: a1_term_disp, a2_term_disp
      real(8) :: eps_ij, sig_ij, kij_val
      real(8) :: ar_dispersion
      ! Attributes declarations as variables
      real(8) :: m(size(n))
      real(8) :: sigma(size(n))
      real(8) :: epsilon_k(size(n))
      real(8) :: kij(size(n),size(n))


      m = model%m
      sigma = model%sigma
      epsilon_k = model%epsilon_k
      kij = model%kij

      ! Residual Helmholtz free energy calculation
      x = n / sum(n)  ! Mole fractions
      n_tot = sum(n)  ! Total moles
      nc = size(n)    ! Number of components

      ! Diameter calculation
      do i = 1, size(n)
         d(i) = sigma(i) * (1.0_8 - 0.12_8 * exp(-3.0_8 * epsilon_k(i) / t))
      end do
      
      ! Zetas calculation
      rho = sum(n) * N_AVO / (V * 1.0e27_8)  ! Number density [1/A^3]
      
      zetas = 0.0_8

      do k = 0, 3
         zetas(k) = rho * (PI / 6.0_8) * sum(x * m * (d**k))
      end do

      ! Ar hard-sphere contribution
      one_m_z3_hs = 1.0_8 - zeta(3)

      term1_hs = (3.0_8 * zeta(1) * zeta(2)) / one_m_z3_hs
      term2_hs = (zeta(2)**3) / (zeta(3) * (one_m_z3_hs**2))
      term3_hs = ((zeta(2)**3)/(zeta(3)**2) - zeta(0)) * log(one_m_z3_hs)

      ar_hs = (1.0_8 / zeta(0)) * (term1_hs + term2_hs + term3_hs) ! A_hs/RT

      ! Ar chain contribution
      ar_chain = 0.0_8
      one_m_z3_chain = 1.0_8 - zeta(3)

      do i = 1, size(x)
         di_2 = d(i) / 2.0_8

         g_ii = (1.0_8 / one_m_z3_chain) + &
            (3.0_8 * di_2 * zeta(2)) / (one_m_z3_chain**2) + &
            (2.0_8 * (di_2**2) * (zeta(2)**2)) / (one_m_z3_chain**3)

         ! Summation
         ar_chain = ar_chain - x(i) * (m(i) - 1.0_8) * log(g_ii)
      end do

      ! Ar dispersion contribution
      rho_disp = n_tot * N_AVO / (V * 1.0e27_pr)  ! Number density [1/A^3]
      eta = zeta(3)
      one_m_eta = 1.0_pr - eta

      ! 1. Mixture Averages (Same as before)
      m_ave = 0.0_pr; m2_es3 = 0.0_pr; m2_e2s3 = 0.0_pr

      do i = 1, nc
         m_ave = m_ave + x(i)*m(i) ! Segment average
      end do

      do i = 1, nc
         do j = 1, nc
            sig_ij = 0.5_pr * (sigma(i) + sigma(j))
            kij_val = 0.0_pr; if (present(kij)) kij_val = kij(i,j)
            eps_ij = sqrt(epsilon_k(i) * epsilon_k(j)) * (1.0_pr - kij_val)

            term_disp = x(i) * x(j) * m(i) * m(j) * (sig_ij**3)

            m2_es3 = m2_es3 + term_disp * (eps_ij / T)
            m2_e2s3 = m2_e2s3 + term_disp * (eps_ij / T)**2
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
      ! a_res/RT = -2*pi*rho_disp*I1*... - pi*rho_disp*m_ave*C1*I2*...
      ! ---------------------------------------------------------

      ! 1st order term (Takes 2 * PI)
      a1_term_disp = -2.0_pr * PI * rho_disp * I1 * m2_es3

      ! 2nd order term (Takes 1 * PI and multiplies by m_ave)
      a2_term_disp = -1.0_pr * PI * rho_disp * m_ave * C1 * I2 * m2_e2s3

      ! Sum and make extensive
      ar_dispersion = a1_term_disp + a2_term_disp   

      ! Final evaluation
      m_average = 0.0_8
      do i = 1, nc
         m_average = m_average + x(i)*m(i) ! Segment average
      end do

      arval = (R * T * n_tot) * (m_average * a_hs + a_chain + a_disp)
   end subroutine ar

   ! ---------------------------------------------------------
   ! Method get_v0: Volume lower limit (Covolume)
   ! Represents the physical volume occupied by the segments
   ! ---------------------------------------------------------
   function get_v0(self, n, P, T) result(v0)
      class(PCSAFT), intent(in) :: self
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
   end function get_v0
end module tapenade_model_template
