program test_boundaries_auxiliar
   !! Unit tests of the phase-boundary helper routines in
   !! `yaeos__equilibria_boundaries_auxiliar`, driven with synthetic
   !! continuation vectors \(X = [\ln K^1, \dots, \ln K^{np}, \beta, \ln P, \ln T]\)
   !! so that every branch is reached deterministically.
   use yaeos, only: pr
   use yaeos__equilibria_boundaries_auxiliar, only: &
      get_z, detect_critical, check_critical_jump, near_critical
   use testing_aux, only: assert, test_title
   implicit none

   write(*, *) test_title("Phase-boundary auxiliar helpers")

   call test_get_z()
   call test_detect_critical_crossing()
   call test_detect_critical_first_point()
   call test_detect_critical_binary_stop()
   call test_detect_critical_binary_crossing_stop()
   call test_detect_critical_step_back()
   call test_check_critical_jump()
   call test_near_critical()

contains

   subroutine test_get_z()
      real(pr) :: z0(3)
      real(pr) :: zi(3)
      real(pr) :: z(3)
      real(pr) :: dzda(3)

      z0 = [0.6_pr, 0.3_pr, 0.1_pr]
      zi = [0.0_pr, 0.2_pr, 0.8_pr]

      call get_z(0.25_pr, z0, zi, z, dzda)
      call assert(maxval(abs(z - (0.25_pr*zi + 0.75_pr*z0))) < 1e-15_pr, &
         "get_z: linear interpolation between z0 and z_inj")
      call assert(abs(sum(z) - 1) < 1e-14_pr, "get_z: normalized endpoints give a normalized z")
      call assert(maxval(abs(dzda - (zi - z0))) < 1e-15_pr, "get_z: dz/dalpha = z_inj - z0")

      call get_z(0.0_pr, z0, zi, z)
      call assert(maxval(abs(z - z0)) < 1e-15_pr, "get_z: alpha=0 returns z0 (dzda omitted)")
      call get_z(1.0_pr, z0, zi, z)
      call assert(maxval(abs(z - zi)) < 1e-15_pr, "get_z: alpha=1 returns z_inj")
   end subroutine test_get_z

   subroutine test_detect_critical_crossing()
      !! Ternary, one main phase. The extrapolated step changes the sign of
      !! every lnK, so a critical point must be detected, the phase kinds
      !! swapped and the point moved to the other side of the critical point.
      integer, parameter :: nc = 3
      integer, parameter :: np = 1
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: Xold(nc*np + np + 2)
      real(pr) :: dXdS(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      real(pr) :: dS
      real(pr) :: S
      integer :: ns
      logical :: found
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w
      real(pr) :: t

      Xold = [0.3_pr, -0.2_pr, 0.25_pr, 1.0_pr, log(50._pr), log(300._pr)]
      X = Xold
      dXdS = [-1.0_pr, 0.6_pr, -0.8_pr, 0.0_pr, 0.1_pr, 0.05_pr]
      dS = 0.4_pr
      ns = 5
      S = X(ns)
      kinds_x = "liquid"
      kind_w = "vapor"

      call detect_critical(nc, np, 2, kinds_x, kind_w, .false., Xold, &
         X, dXdS, ns, dS, S, found, Xc)

      call assert(found, "detect_critical: sign change of all lnK flags a critical point")
      call assert(kinds_x(1) == "vapor" .and. kind_w == "liquid", &
         "detect_critical: main and incipient phase kinds are swapped")

      ! The interpolation is made on the component with the largest change
      ! (component 1, from 0.3 to -0.1), at lnK_1 = 0
      t = 0.3_pr/0.4_pr
      call assert(abs(Xc(1)) < 1e-12_pr, "detect_critical: Xc has lnK=0 for the interpolating component")
      call assert(maxval(abs(Xc - (Xold + t*dXdS*dS))) < 1e-12_pr, &
         "detect_critical: Xc is the linear interpolation between Xold and Xnew")

      ! The new point is beyond the critical point, at least as far from
      ! it as the old one was (in lnK units) and in the direction of dXdS
      call assert(maxval(abs(X(1:nc))) >= maxval(abs(Xold(1:nc))), &
         "detect_critical: restarting point is as far from the CP as the old one")
      call assert(all(X(1:nc)*Xold(1:nc) < 0), "detect_critical: restarting point is on the other side")
      call assert(abs(S - X(ns)) < 1e-15_pr, "detect_critical: specification updated to X(ns)")
      call assert(abs(dS - 0.4_pr) < 1e-15_pr, "detect_critical: step size untouched for nc > 2")
   end subroutine test_detect_critical_crossing

   subroutine test_detect_critical_first_point()
      !! On the first point a sign change must not be considered a CP
      !! (liquid-liquid lines starting from a critical point).
      integer, parameter :: nc = 3
      integer, parameter :: np = 1
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: Xold(nc*np + np + 2)
      real(pr) :: dXdS(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      real(pr) :: dS
      real(pr) :: S
      integer :: ns
      logical :: found
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w

      Xold = [0.3_pr, -0.2_pr, 0.25_pr, 1.0_pr, log(50._pr), log(300._pr)]
      X = Xold
      dXdS = [-1.0_pr, 0.6_pr, -0.8_pr, 0.0_pr, 0.1_pr, 0.05_pr]
      dS = 0.4_pr
      ns = 5
      S = X(ns)
      kinds_x = "liquid"
      kind_w = "liquid"

      call detect_critical(nc, np, 1, kinds_x, kind_w, .false., Xold, &
         X, dXdS, ns, dS, S, found, Xc)
      call assert(.not. found, "detect_critical: no CP reported on the first point")
      call assert(maxval(abs(X - Xold)) < 1e-15_pr, "detect_critical: X untouched when no CP")
   end subroutine test_detect_critical_first_point

   subroutine test_detect_critical_binary_stop()
      !! Binary system sitting on the critical point (|lnK| < 1e-6) with
      !! binary_stop: the continuation must be stopped (dS = 0).
      integer, parameter :: nc = 2
      integer, parameter :: np = 1
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: Xold(nc*np + np + 2)
      real(pr) :: dXdS(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      real(pr) :: dS
      real(pr) :: S
      integer :: ns
      logical :: found
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w

      X = [1e-7_pr, -2e-7_pr, 1.0_pr, log(80._pr), 0.5_pr]
      Xold = X
      dXdS = [1.0_pr, -1.0_pr, 0.0_pr, 0.2_pr, 0.1_pr]
      dS = 0.01_pr
      ns = 4
      S = X(ns)
      kinds_x = "liquid"
      kind_w = "vapor"

      call detect_critical(nc, np, 5, kinds_x, kind_w, .true., Xold, &
         X, dXdS, ns, dS, S, found, Xc)
      call assert(dS == 0, "detect_critical: binary at the CP with binary_stop sets dS=0")
      call assert(.not. found, "detect_critical: binary stop does not report a crossing")
      call assert(kinds_x(1) == "liquid" .and. kind_w == "vapor", &
         "detect_critical: binary stop keeps the phase kinds")
   end subroutine test_detect_critical_binary_stop

   subroutine test_detect_critical_binary_crossing_stop()
      !! Binary crossing a CP with binary_stop: found and stopped.
      integer, parameter :: nc = 2
      integer, parameter :: np = 1
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: Xold(nc*np + np + 2)
      real(pr) :: dXdS(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      real(pr) :: dS
      real(pr) :: S
      integer :: ns
      logical :: found
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w

      X = [0.2_pr, -0.3_pr, 1.0_pr, log(80._pr), 0.5_pr]
      Xold = X
      dXdS = [-1.0_pr, 1.2_pr, 0.0_pr, 0.2_pr, 0.1_pr]
      dS = 0.5_pr
      ns = 4
      S = X(ns)
      kinds_x = "liquid"
      kind_w = "vapor"

      call detect_critical(nc, np, 3, kinds_x, kind_w, .true., Xold, &
         X, dXdS, ns, dS, S, found, Xc)
      call assert(found, "detect_critical: binary crossing is found")
      call assert(dS == 0, "detect_critical: binary crossing with binary_stop sets dS=0")
      call assert(abs(Xc(2)) < 1e-12_pr, "detect_critical: binary Xc interpolated at lnK_2 = 0")
   end subroutine test_detect_critical_binary_crossing_stop

   subroutine test_detect_critical_step_back()
      !! Near-critical point (all |lnK| below the limit). With a large step
      !! the point is pushed backwards (half steps) until it leaves the
      !! near-critical region; with a small step it is pushed forward.
      integer, parameter :: nc = 3
      integer, parameter :: np = 1
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: X0(nc*np + np + 2)
      real(pr) :: Xold(nc*np + np + 2)
      real(pr) :: dXdS(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      real(pr) :: dS
      real(pr) :: S
      integer :: ns
      logical :: found
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w
      real(pr) :: limit

      limit = 0.01_pr + nc*(0.1_pr - 0.01_pr)/18
      X0 = [0.01_pr, -0.01_pr, 0.005_pr, 1.0_pr, log(50._pr), log(300._pr)]
      dXdS = [-1.0_pr, 1.0_pr, -1.0_pr, 0.0_pr, 0.1_pr, 0.05_pr]

      ! Large step: one backward half-step is enough to leave the region
      X = X0
      Xold = X0
      dS = 0.2_pr
      ns = 1
      S = X(ns)
      kinds_x = "liquid"
      kind_w = "vapor"
      call detect_critical(nc, np, 4, kinds_x, kind_w, .false., Xold, &
         X, dXdS, ns, dS, S, found, Xc)
      call assert(found, "detect_critical: big step over a near-critical point is a crossing")
      call assert(kinds_x(1) == "vapor", "detect_critical: big step swaps kinds")

      ! Small step: moved forward until leaving the region, no crossing
      ! is extrapolated because the next step keeps the signs.
      X = X0
      Xold = X0
      dXdS = [1.0_pr, -1.0_pr, 1.0_pr, 0.0_pr, 0.1_pr, 0.05_pr]
      dS = 0.004_pr
      kinds_x = "liquid"
      kind_w = "vapor"
      call detect_critical(nc, np, 4, kinds_x, kind_w, .false., Xold, &
         X, dXdS, ns, dS, S, found, Xc)
      call assert(.not. found, "detect_critical: moving away from the CP is not a crossing")
      call assert(maxval(abs(X(1:nc))) >= min(limit, 0.05_pr), &
         "detect_critical: point pushed out of the near-critical region")
      call assert(X(1) > X0(1) .and. X(2) < X0(2), "detect_critical: small steps go forward along dXdS")
   end subroutine test_detect_critical_step_back

   subroutine test_check_critical_jump()
      integer, parameter :: nc = 3
      integer, parameter :: np = 2
      real(pr) :: X(nc*np + np + 2)
      real(pr) :: Xlast(nc*np + np + 2)
      real(pr) :: Xc(nc*np + np + 2)
      logical :: jumped
      character(len=14) :: kinds_x(np)
      character(len=14) :: kind_w

      ! Only the second main phase crosses its critical point
      Xlast = [0.5_pr, -0.4_pr, 0.3_pr, 0.1_pr, -0.3_pr, 0.05_pr, 0.6_pr, 0.4_pr, log(50._pr), log(300._pr)]
      X =     [0.6_pr, -0.5_pr, 0.2_pr, -0.1_pr, 0.1_pr, -0.05_pr, 0.6_pr, 0.4_pr, log(55._pr), log(305._pr)]
      kinds_x = ["liquid        ", "liquid        "]
      kind_w = "vapor"

      call check_critical_jump(nc, np, 1, kinds_x, kind_w, X, Xlast, Xc, .false., jumped)
      call assert(jumped, "check_critical_jump: crossing of the second phase lnK is detected")
      ! Largest change is component 2 of phase 2 (index 5): from -0.3 to 0.1
      call assert(abs(Xc(5)) < 1e-12_pr, "check_critical_jump: Xc interpolated at lnK=0 of the fastest component")
      call assert(abs(Xc(nc*np+np+1) - (log(50._pr) + 0.75_pr*(log(55._pr) - log(50._pr)))) < 1e-12_pr, &
         "check_critical_jump: lnP interpolated with the same weight")
      call assert(kinds_x(2) == "vapor" .and. kind_w == "liquid" .and. kinds_x(1) == "liquid", &
         "check_critical_jump: only the crossing phase swaps its kind with the incipient phase")

      ! Same crossing but the CP had already been handled: no kind swap
      kinds_x = ["liquid        ", "liquid        "]
      kind_w = "vapor"
      call check_critical_jump(nc, np, 1, kinds_x, kind_w, X, Xlast, Xc, .true., jumped)
      call assert(jumped, "check_critical_jump: jump reported even if already found")
      call assert(kinds_x(2) == "liquid" .and. kind_w == "vapor", &
         "check_critical_jump: kinds untouched when the CP was already found")

      ! No crossing
      call check_critical_jump(nc, np, 1, kinds_x, kind_w, Xlast, Xlast, Xc, .false., jumped)
      call assert(.not. jumped, "check_critical_jump: no jump without a sign change")
   end subroutine test_check_critical_jump

   subroutine test_near_critical()
      integer, parameter :: nc = 3
      integer, parameter :: np = 2
      real(pr) :: X(nc*np + np + 2)
      logical :: near
      integer :: l
      integer :: i
      integer :: j

      X = [1.0_pr, -2.0_pr, 0.5_pr, 0.05_pr, -0.08_pr, 0.02_pr, 0.5_pr, 0.5_pr, log(50._pr), log(300._pr)]
      call near_critical(nc, np, X, near, l, i, j)
      call assert(near, "near_critical: second phase with small lnK is near critical")
      call assert(l == 2, "near_critical: reports the near-critical phase index")
      call assert(i == 1 .and. j == 2, "near_critical: max/min lnK component indices")

      X(4:6) = [0.5_pr, -0.6_pr, 0.2_pr]
      call near_critical(nc, np, X, near, l, i, j)
      call assert(.not. near, "near_critical: far from critical for all phases")
   end subroutine test_near_critical
end program test_boundaries_auxiliar
