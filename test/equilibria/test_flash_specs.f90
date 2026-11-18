program test_flash_specs
   !! Two-phase flash at every supported specification (TP, TV, with and
   !! without initial K factors) for both `flash` and `flash_no_beta_limits`.
   !! Every converged result is checked for material balance, normalized
   !! phase compositions and isofugacity, and the TP and TV routes are
   !! cross-checked against each other. Single-phase states must be reported
   !! as `kind = "failed"` with `P = -1`.
   use yaeos, only: pr, ArModel, PengRobinson76, EquilibriumState, flash
   use yaeos__equilibria_flash, only: flash_no_beta_limits
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: model
   type(EquilibriumState) :: tp
   type(EquilibriumState) :: tv
   type(EquilibriumState) :: res
   type(EquilibriumState) :: tp10
   real(pr) :: z(nc)
   real(pr) :: T
   real(pr) :: P
   real(pr) :: V
   real(pr) :: V10
   real(pr) :: Pz
   integer :: iters

   write(*, *) test_title("Flash specifications and routes")

   ! Methane / n-butane
   model = PengRobinson76( &
      [190.564_pr, 425.12_pr], [45.99_pr, 37.96_pr], [0.0115478_pr, 0.200164_pr])
   z = [0.3_pr, 0.7_pr]
   T = 300.0_pr
   P = 40.0_pr

   ! ===========================================================================
   ! TP flash, reference solution
   ! ---------------------------------------------------------------------------
   tp = flash(model, z, T, p_spec=P, iters=iters)
   call assert(tp%kind == "split", "TP flash: two phases at 300 K, 40 bar")
   call check_solution(tp, "TP flash")
   call assert(abs(tp%P - P) < 1e-12_pr, "TP flash: pressure equals the specification")
   call assert(tp%iters == iters, "TP flash: iters argument equals the stored count")

   ! ===========================================================================
   ! TP flash with user K factors (the converged ones): immediate convergence
   ! ---------------------------------------------------------------------------
   res = flash(model, z, T, p_spec=P, k0=tp%y/tp%x)
   call assert(res%kind == "split", "TP flash with converged k0 is a split")
   call assert(abs(res%beta - tp%beta) < 1e-6_pr, "TP flash with k0: same beta as with Wilson K")
   call assert(res%iters <= 3, "TP flash started at the solution converges in a few iterations")

   ! ===========================================================================
   ! TV flash: specify the total volume of a TP solution and recover P.
   ! At 10 bar the homogeneous-fluid pressure at that volume is positive, so
   ! the Wilson K-factors are evaluated at it.
   ! ---------------------------------------------------------------------------
   tp10 = flash(model, z, T, p_spec=10.0_pr)
   call assert(tp10%kind == "split", "TP flash: two phases at 300 K, 10 bar")
   V10 = tp10%beta*tp10%Vy + (1 - tp10%beta)*tp10%Vx
   call model%pressure(z, V10, T, P=Pz)
   call assert(Pz > 0, "homogeneous pressure at the 10 bar two-phase volume is positive")

   tv = flash(model, z, T, v_spec=V10, iters=iters)
   call assert(tv%kind == "split", "TV flash: two phases at the 10 bar TP-solution volume")
   call check_solution(tv, "TV flash")
   call assert(tv%iters == iters, "TV flash: iters argument equals the stored count")
   call assert(abs(tv%P - 10.0_pr) < 1e-4_pr, "TV flash recovers the TP pressure")
   call assert(abs(tv%beta - tp10%beta) < 1e-5_pr, "TV flash recovers the TP vapor fraction")
   call assert(abs(tv%beta*tv%Vy + (1 - tv%beta)*tv%Vx - V10) < 1e-6_pr*V10, &
      "TV flash: phase volumes add up to the specified volume")

   ! TV flash with K factors provided. At the 40 bar volume the homogeneous
   ! pressure is negative, so the K factors must come from the user.
   V = tp%beta*tp%Vy + (1 - tp%beta)*tp%Vx
   res = flash(model, z, T, v_spec=V, k0=tp%y/tp%x)
   call assert(res%kind == "split", "TV flash with k0 is a split")
   call check_solution(res, "TV flash with k0")
   call assert(abs(res%P - P) < 1e-4_pr, "TV flash with k0 recovers the TP pressure")
   call assert(abs(res%beta - tp%beta) < 1e-5_pr, "TV flash with k0 recovers the TP vapor fraction")

   ! ===========================================================================
   ! Single-phase states: trivial solution flagged as failed
   ! ---------------------------------------------------------------------------
   ! Compressed liquid
   res = flash(model, z, T, p_spec=200.0_pr, iters=iters)
   call assert(res%kind == "failed", "TP flash of a compressed liquid is flagged as failed")
   call assert(res%P == -1, "failed flash returns P = -1")
   call assert(res%T == T, "failed flash keeps the temperature")
   call assert(res%iters == iters, "failed flash stores the iteration count")

   ! ===========================================================================
   ! flash_no_beta_limits
   ! ---------------------------------------------------------------------------
   ! TP
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, p_spec=P)
   call assert(res%kind == "split", "flash_no_beta_limits TP: split")
   call check_solution(res, "flash_no_beta_limits TP")
   call assert(abs(res%beta - tp%beta) < 1e-6_pr, "flash_no_beta_limits TP agrees with flash")

   ! TP with initial K factors
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, p_spec=P, k0=tp%y/tp%x)
   call assert(abs(res%beta - tp%beta) < 1e-6_pr, "flash_no_beta_limits TP with k0 agrees with flash")

   ! TV without K factors
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, v_spec=V10, iters=iters)
   call assert(res%kind == "split", "flash_no_beta_limits TV: split")
   call check_solution(res, "flash_no_beta_limits TV")
   call assert(abs(res%P - 10.0_pr) < 1e-4_pr, "flash_no_beta_limits TV recovers the TP pressure")
   call assert(abs(res%beta - tp10%beta) < 1e-5_pr, "flash_no_beta_limits TV recovers the TP beta")
   call assert(res%iters == iters, "flash_no_beta_limits TV: iters argument equals stored count")

   ! TV with K factors
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, v_spec=V, k0=tp%y/tp%x)
   call assert(res%kind == "split", "flash_no_beta_limits TV with k0: split")
   call check_solution(res, "flash_no_beta_limits TV with k0")
   call assert(abs(res%P - P) < 1e-4_pr, "flash_no_beta_limits TV with k0 recovers the TP pressure")

   ! Compressed liquid: trivial solution flagged as failed
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, p_spec=200.0_pr, iters=iters)
   call assert(res%kind == "failed", "flash_no_beta_limits of a compressed liquid is flagged as failed")
   call assert(res%P == -1, "flash_no_beta_limits failure returns P = -1")
   call assert(res%iters == iters, "flash_no_beta_limits failure stores the iteration count")

   ! Pressure and volume both specified: rejected, nothing computed
   res = flash_no_beta_limits(model, z, T, beta0=0.5_pr, v_spec=V, p_spec=P)
   call assert(.not. allocated(res%x), "flash_no_beta_limits rejects specifying both P and V")
   res = flash(model, z, T, v_spec=V, p_spec=P)
   call assert(.not. allocated(res%x), "flash rejects specifying both P and V")

   call test_ge_model()

contains

   subroutine test_ge_model()
      !! Activity-coefficient (GeModel) flash: only T can be specified and
      !! initial K factors are required. Both flash routines must agree and
      !! satisfy isoactivity.
      use yaeos, only: UNIFAC, Groups, setup_unifac, min_tpd
      integer, parameter :: ncg = 3
      type(UNIFAC) :: ge
      type(Groups) :: molecules(ncg)
      type(EquilibriumState) :: ref
      type(EquilibriumState) :: nbl
      real(pr) :: n(ncg)
      real(pr) :: wt(ncg)
      real(pr) :: k0(ncg)
      real(pr) :: lng_x(ncg)
      real(pr) :: lng_y(ncg)
      real(pr) :: tm_min
      real(pr) :: Tge

      Tge = 250.0_pr
      n = [0.2_pr, 0.7_pr, 0.1_pr]
      molecules(1)%groups_ids = [1, 42]
      molecules(1)%number_of_groups = [1, 1]
      molecules(2)%groups_ids = [1, 2]
      molecules(2)%number_of_groups = [2, 4]
      molecules(3)%groups_ids = [16]
      molecules(3)%number_of_groups = [1]
      ge = setup_unifac(molecules)

      call min_tpd(ge, n, P=1.0_pr, T=Tge, mintpd=tm_min, w=wt)
      call assert(tm_min < 0, "GeModel: feed is unstable (negative tm)")
      k0 = wt/n

      ref = flash(ge, n, Tge, k0=k0)
      call assert(ref%kind == "split", "GeModel flash: liquid-liquid split")

      ! Started from the converged K factors (from the min_tpd trial phase
      ! flash_no_beta_limits stalls at beta ~ 5e-7 with sum(y) ~ 2, since
      ! solve_rr clamps beta to the K-dependent limits).
      nbl = flash_no_beta_limits(ge, n, Tge, beta0=0.5_pr, k0=ref%y/ref%x)
      call assert(nbl%kind == "split", "GeModel flash_no_beta_limits: liquid-liquid split")
      call assert(abs(nbl%beta - ref%beta) < 1e-6_pr, "GeModel: both flash routines give the same beta")
      call assert(maxval(abs(nbl%x - ref%x)) < 1e-6_pr, "GeModel: both flash routines give the same x")
      call assert(maxval(abs(nbl%beta*nbl%y + (1 - nbl%beta)*nbl%x - n)) < 1e-10_pr, &
         "GeModel flash_no_beta_limits: material balance")
      call ge%ln_activity_coefficient(nbl%x, Tge, lngamma=lng_x)
      call ge%ln_activity_coefficient(nbl%y, Tge, lngamma=lng_y)
      call assert(maxval(abs(log(nbl%x) + lng_x - log(nbl%y) - lng_y)) < 1e-5_pr, &
         "GeModel flash_no_beta_limits: isoactivity")

      ! Pressure or volume specifications are rejected for GeModels
      nbl = flash(ge, n, Tge, p_spec=1.0_pr, k0=k0)
      call assert(.not. allocated(nbl%x), "GeModel flash rejects a pressure specification")
      nbl = flash_no_beta_limits(ge, n, Tge, beta0=0.5_pr, v_spec=1.0_pr, k0=k0)
      call assert(.not. allocated(nbl%x), "GeModel flash_no_beta_limits rejects a volume specification")
   end subroutine test_ge_model

   subroutine check_solution(fr, name)
      !! Material balance, normalization and isofugacity of a flash result.
      type(EquilibriumState), intent(in) :: fr
      character(len=*), intent(in) :: name
      real(pr) :: lnphi_x(nc)
      real(pr) :: lnphi_y(nc)
      real(pr) :: Vx
      real(pr) :: Vy

      call assert(fr%beta > 0 .and. fr%beta < 1, name // ": 0 < beta < 1")
      call assert(abs(sum(fr%x) - 1) < 1e-10_pr, name // ": liquid fractions sum to one")
      call assert(abs(sum(fr%y) - 1) < 1e-10_pr, name // ": vapor fractions sum to one")
      call assert(maxval(abs(fr%beta*fr%y + (1 - fr%beta)*fr%x - z)) < 1e-10_pr, &
         name // ": material balance")
      call assert(maxval(abs(fr%y - fr%x)) > 1e-2_pr, name // ": non-trivial split")

      call model%lnphi_pt(fr%x, fr%P, fr%T, V=Vx, root_type="liquid", lnPhi=lnphi_x)
      call model%lnphi_pt(fr%y, fr%P, fr%T, V=Vy, root_type="vapor", lnPhi=lnphi_y)
      call assert(maxval(abs((log(fr%x) + lnphi_x) - (log(fr%y) + lnphi_y))) < 1e-5_pr, &
         name // ": isofugacity between the phases")
      call assert(abs(Vx - fr%Vx) < 1e-6_pr*Vx, name // ": stored liquid volume is the liquid root")
      call assert(abs(Vy - fr%Vy) < 1e-6_pr*Vy, name // ": stored vapor volume is the vapor root")
   end subroutine check_solution
end program test_flash_specs
