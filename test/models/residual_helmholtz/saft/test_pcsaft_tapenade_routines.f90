program test_pcsaft_tapenade_routines
   !! Direct checks of the Tapenade-generated PC-SAFT routines that the
   !! `ArModelTapenade` adapter does not call: `AR_D` (first-order tangent) and
   !! `AR_D_D_D` (third-order nested tangent).
   !!
   !! Conventions (verified numerically): direction "d" seeds (n, V, T),
   !! directions "0" and "1" seed (V, T).
   !! * AR_D_D_D's `arval, arvald0, arvald, arvaldd` equal AR_D_D's outputs.
   !! * `arvald1, arvald0d, arvaldd0, arvalddd` are the derivatives along
   !!   direction 1 of `arval, arvald0, arvald, arvaldd`.
   use yaeos, only: pr
   use yaeos__models_ar_saft_pcsaft, only: PCSAFT, init_pcsaft, AR_D, AR_D_D, AR_D_D_D
   use testing_aux, only: assert, test_title
   implicit none
   type(PCSAFT) :: model
   integer :: icase

   write(*, *) test_title("PC-SAFT Tapenade tangent routines")

   ! Methane + n-hexane parameters (Gross & Sadowski 2001), states from dilute
   ! gas to dense liquid, with and without binary interaction.
   do icase = 1, 4
      select case (icase)
       case (1)
         model = init_pcsaft([1.0000_pr], [3.7039_pr], [150.03_pr])
         call check_state("pure, dilute", [1.0_pr], [0.2_pr], 20.0_pr, 250.0_pr)
       case (2)
         model = init_pcsaft([1.0000_pr, 3.0576_pr], [3.7039_pr, 3.7983_pr], [150.03_pr, 236.77_pr])
         call check_state("binary, kij = 0, dense", [0.4_pr, 0.6_pr], [0.3_pr, -0.2_pr], 0.12_pr, 300.0_pr)
       case (3)
         model = init_pcsaft([1.0000_pr, 3.0576_pr], [3.7039_pr, 3.7983_pr], [150.03_pr, 236.77_pr], &
            reshape([0.0_pr, 0.02_pr, 0.02_pr, 0.0_pr], [2, 2]))
         call check_state("binary, kij /= 0, moderate", [1.4_pr, 0.6_pr], [0.1_pr, 0.5_pr], 1.5_pr, 400.0_pr)
       case (4)
         model = init_pcsaft([1.0000_pr, 2.0020_pr, 3.0576_pr], [3.7039_pr, 3.6184_pr, 3.7983_pr], &
            [150.03_pr, 208.11_pr, 236.77_pr])
         call check_state("ternary, liquid-like", [0.2_pr, 0.3_pr, 0.5_pr], [0.05_pr, -0.1_pr, 0.2_pr], &
            0.09_pr, 320.0_pr)
      end select
   end do

contains

   subroutine check_state(label, n, nd, v, t)
      character(len=*), intent(in) :: label !! Case description for messages.
      real(pr), intent(in) :: n(:) !! Mole amounts [mol].
      real(pr), intent(in) :: nd(:) !! Composition seed of direction d.
      real(pr), intent(in) :: v !! Volume [L].
      real(pr), intent(in) :: t !! Temperature [K].

      real(pr), parameter :: h = 1.0e-5_pr !! Relative finite-difference step.
      real(pr), parameter :: tol = 1.0e-6_pr !! Relative tolerance for central differences.
      real(pr) :: vd !! Volume seed of direction d [L].
      real(pr) :: td !! Temperature seed of direction d [K].
      real(pr) :: vd0 !! Volume seed of direction 0 [L].
      real(pr) :: td0 !! Temperature seed of direction 0 [K].
      real(pr) :: vd1 !! Volume seed of direction 1 [L].
      real(pr) :: td1 !! Temperature seed of direction 1 [K].
      real(pr) :: zero_n(size(n)) !! Zero composition seed.
      real(pr) :: a3(8) !! AR_D_D_D outputs in argument order.
      real(pr) :: a2(4) !! AR_D_D outputs: arval, arvald0, arvald, arvaldd.
      real(pr) :: a2p(4) !! AR_D_D outputs at +h along direction 1.
      real(pr) :: a2m(4) !! AR_D_D outputs at -h along direction 1.
      real(pr) :: fd1(4) !! Central differences of a2 along direction 1.
      real(pr) :: ar !! AR_D primal value.
      real(pr) :: ard !! AR_D tangent along d.
      real(pr) :: arp !! Ar at +h along d.
      real(pr) :: arm !! Ar at -h along d.
      real(pr) :: dummy !! Unused tangent output.

      vd = 0.05_pr * v
      td = 0.01_pr * t
      vd0 = 0.1_pr * v
      td0 = -0.005_pr * t
      vd1 = -0.07_pr * v
      td1 = 0.01_pr * t
      zero_n = 0.0_pr

      ! AR_D against a central difference of Ar along d.
      call AR_D(model, n, nd, v, vd, t, td, ar, ard)
      call AR_D(model, n + h*nd, zero_n, v + h*vd, 0.0_pr, t + h*td, 0.0_pr, arp, dummy)
      call AR_D(model, n - h*nd, zero_n, v - h*vd, 0.0_pr, t - h*td, 0.0_pr, arm, dummy)
      call assert(close_to(ard, (arp - arm)/(2*h), tol), label // ": AR_D tangent vs finite difference")

      ! AR_D_D_D against AR_D_D and its derivative along direction 1.
      call AR_D_D_D(model, n, nd, v, vd1, vd0, vd, t, td1, td0, td, &
         a3(1), a3(2), a3(3), a3(4), a3(5), a3(6), a3(7), a3(8))
      call AR_D_D(model, n, nd, v, vd0, vd, t, td0, td, a2(1), a2(2), a2(3), a2(4))
      call assert(abs(a3(1) - ar) <= 1e-12_pr * abs(ar), label // ": AR_D_D_D primal equals AR_D primal")
      call assert(close_to(a3(3), a2(2), 1e-12_pr), label // ": AR_D_D_D arvald0 equals AR_D_D")
      call assert(close_to(a3(5), a2(3), 1e-12_pr), label // ": AR_D_D_D arvald equals AR_D_D")
      call assert(close_to(a3(7), a2(4), 1e-12_pr), label // ": AR_D_D_D arvaldd equals AR_D_D")

      call AR_D_D(model, n, nd, v + h*vd1, vd0, vd, t + h*td1, td0, td, a2p(1), a2p(2), a2p(3), a2p(4))
      call AR_D_D(model, n, nd, v - h*vd1, vd0, vd, t - h*td1, td0, td, a2m(1), a2m(2), a2m(3), a2m(4))
      fd1 = (a2p - a2m)/(2*h)
      call assert(close_to(a3(2), fd1(1), tol), label // ": arvald1 vs finite difference")
      call assert(close_to(a3(4), fd1(2), tol), label // ": arvald0d vs finite difference")
      call assert(close_to(a3(6), fd1(3), tol), label // ": arvaldd0 vs finite difference")
      call assert(close_to(a3(8), fd1(4), tol), label // ": arvalddd vs finite difference")
   end subroutine check_state

   logical function close_to(x, ref, rtol)
      !! Relative comparison with an absolute floor for values near zero.
      real(pr), intent(in) :: x !! Computed value.
      real(pr), intent(in) :: ref !! Reference value.
      real(pr), intent(in) :: rtol !! Relative tolerance.
      close_to = abs(x - ref) <= rtol * max(abs(ref), 1e-8_pr)
   end function close_to
end program test_pcsaft_tapenade_routines
