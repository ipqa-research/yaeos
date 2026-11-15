program test_armodel_generic_volume
   !! The generic `ArModel%volume` Newton solver (used by models that do not
   !! override it, here a hyperdual-autodiff PR76) must give the same roots
   !! as the Michelsen solver of the equivalent `CubicEoS`, and select the
   !! stable root. Also checks the `lnPhiP` output of `lnphi_vt`.
   use yaeos, only: pr, R, ArModel
   use fixtures_models, only: binary_PR76_hd, binary_PR76
   use testing_aux, only: assert, test_title
   implicit none
   integer, parameter :: nc = 2
   class(ArModel), allocatable :: hd
   class(ArModel), allocatable :: cub
   real(pr) :: z(nc)

   write(*, *) test_title("Generic ArModel volume solver")

   hd = binary_PR76_hd()
   cub = binary_PR76()
   z = [0.3_pr, 0.7_pr]

   ! Liquid roots, from low to very high pressure
   call compare_root(200.0_pr, 1.0_pr, "liquid")
   call compare_root(200.0_pr, 100.0_pr, "liquid")
   call compare_root(200.0_pr, 1000.0_pr, "liquid")

   ! Vapor root at low pressure
   call compare_root(200.0_pr, 1.0_pr, "vapor")

   ! Stable root: vapor at 200 K and 1 bar, liquid at 150 K and 2 bar (where
   ! a vapor root also exists).
   call compare_root(200.0_pr, 1.0_pr, "stable")
   call compare_root(150.0_pr, 2.0_pr, "stable")
   call check_stable_choice(200.0_pr, 1.0_pr, "vapor")
   call check_stable_choice(150.0_pr, 2.0_pr, "liquid")

   call test_lnphip()

contains

   subroutine compare_root(T, P, root)
      real(pr), intent(in) :: T
      real(pr), intent(in) :: P
      character(len=*), intent(in) :: root
      real(pr) :: V_hd
      real(pr) :: V_cub
      real(pr) :: Pcalc
      real(pr) :: dPdV
      character(len=60) :: name

      write(name, '(A, " root at ", F6.1, " K, ", F7.1, " bar")') root, T, P
      call hd%volume(z, P, T, V_hd, root)
      call cub%volume(z, P, T, V_cub, root)
      call assert(V_hd > 0, trim(name) // ": generic solver converges")
      call assert(abs(V_hd - V_cub) < 1e-8_pr*V_cub, trim(name) // ": same volume as the Michelsen solver")
      call hd%pressure(z, V_hd, T, P=Pcalc, dPdV=dPdV)
      call assert(abs(Pcalc - P) < 1e-6_pr*P, trim(name) // ": P(V) equals the specified pressure")
      call assert(dPdV < 0, trim(name) // ": mechanically stable root")
   end subroutine compare_root

   subroutine check_stable_choice(T, P, expected)
      !! The stable root is the one with the lower residual Gibbs energy.
      real(pr), intent(in) :: T
      real(pr), intent(in) :: P
      character(len=*), intent(in) :: expected
      real(pr) :: V_stable
      real(pr) :: V_root
      real(pr) :: V_liq
      real(pr) :: V_vap
      real(pr) :: G_liq
      real(pr) :: G_vap

      call hd%volume(z, P, T, V_stable, "stable")
      call hd%volume(z, P, T, V_root, expected)
      call assert(abs(V_stable - V_root) < 1e-10_pr*V_root, "stable root is the " // expected // " root")

      call hd%volume(z, P, T, V_liq, "liquid")
      call hd%volume(z, P, T, V_vap, "vapor")
      call assert(V_vap > 10*V_liq, "two distinct roots exist")
      call hd%gibbs_residual_vt(z, V_liq, T, Gr=G_liq)
      call hd%gibbs_residual_vt(z, V_vap, T, Gr=G_vap)
      ! Gr(T, V) differs from Gr(T, P) by n R T ln(Z), add it to compare
      G_liq = G_liq - sum(z)*R*T*log(P*V_liq/(sum(z)*R*T))
      G_vap = G_vap - sum(z)*R*T*log(P*V_vap/(sum(z)*R*T))
      if (expected == "liquid") then
         call assert(G_liq < G_vap, "liquid is stable: lower Gibbs energy")
      else
         call assert(G_vap < G_liq, "vapor is stable: lower Gibbs energy")
      end if
   end subroutine check_stable_choice

   subroutine test_lnphip()
      !! `lnPhiP` is \(\ln(\phi_i P)\), both in the fast path (only `lnPhi`
      !! or `lnPhiP` requested) and with derivatives requested.
      real(pr) :: V
      real(pr) :: T
      real(pr) :: P
      real(pr) :: lnphi(nc)
      real(pr) :: lnphip(nc)
      real(pr) :: lnphip_full(nc)
      real(pr) :: dlnphidn(nc, nc)

      T = 200.0_pr
      call cub%volume(z, 1.0_pr, T, V, "vapor")

      call cub%lnphi_vt(z, V, T, P=P, lnPhi=lnphi)
      call cub%lnphi_vt(z, V, T, lnPhiP=lnphip)
      call assert(maxval(abs(lnphip - (lnphi + log(P)))) < 1e-12_pr, "lnPhiP fast path equals lnPhi + ln P")

      call cub%lnphi_vt(z, V, T, lnPhiP=lnphip_full, dlnPhidn=dlnphidn)
      call assert(maxval(abs(lnphip_full - lnphip)) < 1e-12_pr, "lnPhiP with derivatives equals the fast path")
   end subroutine test_lnphip
end program test_armodel_generic_volume
