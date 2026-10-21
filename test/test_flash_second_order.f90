program test_flash_second_order
   use yaeos, only: pr, EquilibriumState, flash, flash_second_order
   use yaeos, only: ArModel, PengRobinson76, PengRobinson78
   use yaeos, only: GeModel, NRTL, UNIFAC, Groups, setup_unifac
   use yaeos__equilibria_stability, only: min_tpd
   use testing_aux, only: assert, test_title

   implicit none

   print *, test_title("Second-order Gibbs flash")

   call test_pr76_vle()
   call test_pr78_ternary()
   call test_unifac_lle()
   call test_nrtl_lle()
   call test_ar_lle_liquid_root()

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
      z = [0.2_pr, 0.7_pr, 0.1_pr]

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
end program test_flash_second_order
