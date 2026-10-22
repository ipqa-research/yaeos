program test_flash_second_order
   use yaeos, only: pr, EquilibriumState, flash, flash_second_order
   use yaeos, only: ArModel, PengRobinson76, PengRobinson78
   use yaeos, only: GeModel, NRTL, UNIFAC, Groups, setup_unifac
   use yaeos, only: UNIQUAC, setup_uniquac
   use yaeos__equilibria_stability, only: min_tpd
   use testing_aux, only: assert, test_title
   use auxiliar_functions, only: allclose, rel_error

   implicit none

   print *, test_title("Second-order Gibbs flash")

   call test_pr76_vle()
   call test_pr78_ternary()
   call test_unifac_lle()
   call test_nrtl_lle()
   call test_ar_lle_liquid_root()
   call test_uniquac_ss_failure_point()

contains

   subroutine test_pr76_vle()
      !! Methane/n-butane with Peng-Robinson 76 (from `test_flash_pt`).
      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: ss, so
      real(pr) :: tc(nc), pc(nc), w(nc), z(nc), k0(nc)
      real(pr) :: T, P
      integer :: its

      tc = [190.564_pr, 425.12_pr]
      pc = [45.99_pr, 37.96_pr]
      w = [0.0115478_pr, 0.200164_pr]
      model = PengRobinson76(tc, pc, w)

      z = [0.4_pr, 0.6_pr]
      T = 294._pr
      P = 60._pr
      k0 = (pc/P)*exp(5.373_pr*(1 + w)*(1 - tc/T))

      ss = flash(model, z, T, p_spec=P, k0=k0, iters=its)
      so = flash_second_order(model, z, T, P, k0=k0, iters=its)

      call assert(ss%kind == "split", "PR76: successive substitutions")
      call assert(so%kind == "split", "PR76: second order converges")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "PR76: liquid composition")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "PR76: vapour composition")
      call assert(abs(so%beta - ss%beta) < 1e-5_pr, "PR76: vapour fraction")
      call assert(abs(so%Vx - ss%Vx) < 1e-5_pr, "PR76: liquid volume")
      call assert(abs(so%Vy - ss%Vy) < 1e-5_pr, "PR76: vapour volume")

      ! Wilson initialization is used when K0 is not provided.
      so = flash_second_order(model, z, T, P)
      call assert(so%kind == "split", "PR76: Wilson initialization")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "PR76: Wilson liquid")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "PR76: Wilson vapour")
   end subroutine test_pr76_vle

   subroutine test_pr78_ternary()
      !! Ternary PR78 mixture used in the multiphase flash tests.
      integer, parameter :: nc = 3

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: ss, so
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc, nc), z(nc)
      real(pr) :: T, P
      integer :: its

      tc = [190.564_pr, 304.1282_pr, 768.0_pr]
      pc = [45.992_pr, 73.773_pr, 10.7_pr]
      w = [0.01142_pr, 0.22394_pr, 0.8805_pr]

      kij = 0.0_pr
      kij(1, 2) = 0.1_pr
      kij(2, 1) = 0.1_pr
      kij(2, 3) = 0.2_pr
      kij(3, 2) = 0.2_pr

      model = PengRobinson78(tc, pc, w, kij=kij)

      z = [0.2_pr, 0.4_pr, 0.4_pr]
      z = z/sum(z)
      T = 260._pr
      P = 150._pr

      ss = flash(model, z, T, p_spec=P, iters=its)
      so = flash_second_order(model, z, T, P, iters=its)

      call assert(ss%kind == "split", "PR78: successive substitutions")
      call assert(so%kind == "split", "PR78: second order converges")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "PR78: liquid composition")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "PR78: vapour composition")
      call assert(abs(so%beta - ss%beta) < 1e-5_pr, "PR78: vapour fraction")
   end subroutine test_pr78_ternary

   subroutine test_unifac_lle()
      !! UNIFAC liquid-liquid split from `test_flash_ge`.
      integer, parameter :: nc = 3

      type(UNIFAC) :: model
      type(Groups) :: molecules(nc)
      type(EquilibriumState) :: ss, so
      real(pr) :: z(nc), w(nc), k0(nc)
      real(pr) :: T, P, mintpd
      integer :: its

      T = 250._pr
      P = 1._pr
      z = [0.2_pr, 0.5_pr, 0.3_pr]

      molecules(1)%groups_ids = [1, 42]
      molecules(1)%number_of_groups = [1, 1]
      molecules(2)%groups_ids = [1, 2]
      molecules(2)%number_of_groups = [2, 4]
      molecules(3)%groups_ids = [16]
      molecules(3)%number_of_groups = [1]
      model = setup_unifac(molecules)

      call min_tpd(model, z, P, T, mintpd, w)
      k0 = w/z

      ss = flash(model, z, T, k0=k0, iters=its)

      ! No K0: the second-order flash initializes itself with stability analysis.
      so = flash_second_order(model, z, T, P, iters=its)

      call assert(ss%kind == "split", "UNIFAC: successive substitutions")
      call assert(so%kind == "split", "UNIFAC: second order converges")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "UNIFAC: heavy composition")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "UNIFAC: light composition")
      call assert(abs(so%beta - ss%beta) < 1e-5_pr, "UNIFAC: phase fraction")
   end subroutine test_unifac_lle

   subroutine test_nrtl_lle()
      !! Strongly non-ideal binary NRTL mixture with two liquid phases.
      integer, parameter :: nc = 2

      type(NRTL) :: model
      type(EquilibriumState) :: ss, so
      real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc)
      real(pr) :: z(nc), w(nc), k0(nc), T, P, mintpd
      integer :: its

      a = 0.0_pr
      b = 0.0_pr
      c = 0.0_pr
      a(1, 2) = 5.0_pr
      a(2, 1) = 5.0_pr
      c(1, 2) = 0.3_pr
      c(2, 1) = 0.3_pr
      model = NRTL(a, b, c)

      z = [0.3_pr, 0.7_pr]
      T = 300._pr
      P = 1.0_pr

      call min_tpd(model, z, P, T, mintpd, w)
      k0 = w/z

      ss = flash(model, z, T, k0=k0, iters=its)
      so = flash_second_order(model, z, T, P, k0=k0, iters=its)

      call assert(ss%kind == "split", "NRTL: successive substitutions")
      call assert(so%kind == "split", "NRTL: second order converges")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "NRTL: heavy composition")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "NRTL: light composition")
      call assert(abs(so%beta - ss%beta) < 1e-5_pr, "NRTL: phase fraction")
   end subroutine test_nrtl_lle

   subroutine test_ar_lle_liquid_root()
      !! Cubic model with a liquid-liquid split. The light phase is forced to
      !! the liquid root, which is the purpose of the `light_root` argument.
      integer, parameter :: nc = 2

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: ss, so
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc, nc)
      real(pr) :: z(nc), wtrial(nc), k0(nc), T, P, mintpd
      integer :: its

      tc = [500.0_pr, 500.0_pr]
      pc = [30.0_pr, 30.0_pr]
      w = [0.3_pr, 0.3_pr]
      kij = 0.0_pr
      kij(1, 2) = 0.2_pr
      kij(2, 1) = 0.2_pr
      model = PengRobinson76(tc, pc, w, kij)

      z = [0.5_pr, 0.5_pr]
      T = 300._pr
      P = 1.0_pr

      call min_tpd(model, z, P, T, mintpd, wtrial)
      k0 = wtrial/z

      ss = flash(model, z, T, p_spec=P, k0=k0, iters=its)
      so = flash_second_order( &
         model, z, T, P, k0=k0, light_root="liquid", iters=its &
         )

      call assert(ss%kind == "split", "Ar LLE: successive substitutions")
      call assert(so%kind == "split", "Ar LLE: second order converges")
      call assert(maxval(abs(so%x - ss%x)) < 1e-5_pr, "Ar LLE: heavy composition")
      call assert(maxval(abs(so%y - ss%y)) < 1e-5_pr, "Ar LLE: light composition")
      call assert(abs(so%beta - ss%beta) < 1e-5_pr, "Ar LLE: phase fraction")
      call assert(so%Vx < 0.2_pr, "Ar LLE: heavy phase is liquid")
      call assert(so%Vy < 0.2_pr, "Ar LLE: light phase is liquid")
   end subroutine test_ar_lle_liquid_root

   subroutine test_uniquac_ss_failure_point()
      !! UNIQUAC point of the H2O2 feed trajectory where the
      !! successive-substitution flash does not converge. The comparison is
      !! against the Gibbs reference, not against `flash`.
      !!
      !! Component order: F, H2O2, PF, H2O, oil, ESBO, polyol.
      integer, parameter :: nc = 7

      type(UNIQUAC) :: model
      type(EquilibriumState) :: so, trivial
      real(pr) :: qs(nc), rs(nc)
      real(pr) :: aij(nc, nc), bij(nc, nc), cij(nc, nc)
      real(pr) :: dij(nc, nc), eij(nc, nc)
      real(pr) :: n(nc), z(nc), k0(nc), ones(nc)
      real(pr) :: x_ref(4), y_ref(4), beta_ref
      real(pr) :: T, P
      real(pr) :: ax(nc), ay(nc), ln_gamma_x(nc), ln_gamma_y(nc)
      integer :: its
      real(8) :: t1, t2

      n = [ &
         8.69085708798058754e+00_pr, 1.04599576389094473e+00_pr, &
         8.06134619967180610e-01_pr, 2.12080199832480858e+01_pr, &
         4.53403810029077592e-01_pr, 4.54221490480090750e-01_pr, &
         8.39249341818703722e-14_pr ]

      T = 333.15_pr
      P = 1.0_pr
      z = n/sum(n)

      qs = [1.532_pr, 2.4_pr, 2.388_pr, 1.4_pr, 31.347_pr, &
         30.492_pr, 41.232_pr]

      rs = [1.528_pr, 2.0_pr, 2.242_pr, 0.92_pr, 38.7157_pr, &
         38.8207_pr, 48.8112_pr]

      aij = 0.0_pr
      aij(1, 4) = 0.0650752585_pr
      aij(4, 1) = -1.82683869_pr

      bij = 0.0_pr

      bij(1, :) = [0.0_pr, 151.0_pr, 135.49896199_pr, &
         59.3246124_pr, -11095.62551048_pr, -3246.2424174_pr, &
         -995.317802_pr]
      bij(2, :) = [-199.0_pr, 0.0_pr, -101.29550998_pr, &
         -11.76497438_pr, -62.99691978_pr, -1940.97729389_pr, &
         -675.950034_pr]
      bij(3, :) = [52.89079805_pr, -3.03767903_pr, 0.0_pr, &
         107.87953158_pr, -1444.81779555_pr, -1113.00536155_pr, &
         -787.3298_pr]
      bij(4, :) = [947.601438_pr, 173.79237634_pr, &
         -210.13173272_pr, 0.0_pr, 158.23356021_pr, &
         190.74747368_pr, 180.790188_pr]
      bij(5, :) = [21.69729307_pr, -169.18769991_pr, &
         207.5615499_pr, -1252.35267809_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      bij(6, :) = [346.55860276_pr, 142.25776765_pr, &
         414.04790322_pr, -1308.66514598_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      bij(7, :) = [264.640828_pr, 37.832309_pr, 340.930152_pr, &
         -1362.762912_pr, 0.0_pr, 0.0_pr, 0.0_pr]

      cij = 0.0_pr
      dij = 0.0_pr
      eij = 0.0_pr

      model = setup_uniquac(qs, rs, aij, bij, cij, dij, eij)

      k0 = [1.04947543188346137e+00_pr, &
         5.59102656162917544e-01_pr, 2.05883043639840357e+00_pr, &
         2.63674340469321811e-01_pr, 4.31976399332579896e+07_pr, &
         1.79109736675408584e+01_pr, 4.05540135724048734e-02_pr]

      x_ref = [0.2656_pr, 0.0326_pr, 0.0236_pr, 0.6700_pr]
      y_ref = [0.2787_pr, 0.0182_pr, 0.0487_pr, 0.1767_pr]
      beta_ref = 0.0418376_pr

      ! The successive-substitution flash does not converge at this point,
      ! so the second-order result is checked against the Gibbs reference.
      so = flash_second_order(model, z, T, P, k0=k0, iters=its)
      
      call assert(so%kind == "split", "UNIQUAC SS-fail: split")
      call assert(abs(so%beta - beta_ref) < 1e-3_pr, &
         "UNIQUAC SS-fail: beta")
      call assert(maxval(abs(so%x(1:4) - x_ref)) < 1e-3_pr, &
         "UNIQUAC SS-fail: x")
      call assert(maxval(abs(so%y(1:4) - y_ref)) < 1e-3_pr, &
         "UNIQUAC SS-fail: y")
      call assert(maxval(abs(so%x - so%y)) > 1e-2_pr, &
         "UNIQUAC SS-fail: non-trivial")

      call model%ln_activity_coefficient(so%x, T, ln_gamma_x)
      call model%ln_activity_coefficient(so%y, T, ln_gamma_y)

      ax = exp(ln_gamma_x) * so%x
      ay = exp(ln_gamma_y) * so%y

      call assert(allclose(ax, ay, 1e-6_pr), "Iso activities K0")

      ! The stability-analysis initialization reaches the same split.
      so = flash_second_order(model, z, T, P, iters=its)
      
      call assert(so%kind == "split", "UNIQUAC SS-fail: min_tpd init")
      call assert(maxval(abs(so%x(1:4) - x_ref)) < 1e-3_pr, &
         "UNIQUAC SS-fail: x")
      call assert(maxval(abs(so%y(1:4) - y_ref)) < 1e-3_pr, &
         "UNIQUAC SS-fail: y")
      call assert(abs(so%beta - beta_ref) < 1e-3_pr, &
         "UNIQUAC SS-fail: min_tpd beta")
      
      call model%ln_activity_coefficient(so%x, T, ln_gamma_x)
      call model%ln_activity_coefficient(so%y, T, ln_gamma_y)

      ax = exp(ln_gamma_x) * so%x
      ay = exp(ln_gamma_y) * so%y

      call assert(allclose(ax, ay, 1e-6_pr), "Iso activities K0")

      ! K0 = 1 is exactly the trivial solution and must be rejected.
      ones = 1.0_pr
      trivial = flash_second_order(model, z, T, P, k0=ones, iters=its)
      call assert(trivial%kind == "failed", &
         "UNIQUAC SS-fail: trivial solution rejected")
   end subroutine test_uniquac_ss_failure_point
end program test_flash_second_order
