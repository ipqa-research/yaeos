program test_cmr_cubic_eos
   !! Cubic mixing rule (`CMR`) used inside a `CubicEoS`.
   !!
   !! - Limiting cases: zero interactions give linear \(B\); a pure component
   !!   gives \(D = n^2 a_i (1 - k_{iii})\) and \(B = n b_i (1 - l_{iii})\).
   !! - Homogeneity (Euler) of \(D\) (degree 2) and \(B\) (degree 1).
   !! - Analytical derivatives of \(D\) and \(B\) against finite differences.
   !! - The \(\Delta_1\) mixing (RKPR rule) for constant \(\delta_1\).
   !! - Residual Helmholtz derivatives of the full model against finite differences.
   use yaeos, only: pr, CubicEoS, PengRobinson76, CMR
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(CubicEoS) :: base_model
   type(CubicEoS) :: model
   type(CubicEoS) :: scaled
   type(CMR) :: mixrule
   type(CMR) :: zero_rule

   real(pr) :: Tc(nc)
   real(pr) :: Pc(nc)
   real(pr) :: w(nc)
   real(pr) :: k(nc, nc, nc)
   real(pr) :: l(nc, nc, nc)
   real(pr) :: n(nc)
   real(pr) :: T
   real(pr) :: V

   write(*, *) test_title("CMR mixing rule inside a CubicEoS")

   Tc = [190.564_pr, 369.83_pr, 617.7_pr]
   Pc = [45.99_pr, 42.48_pr, 21.1_pr]
   w = [0.0115478_pr, 0.1523_pr, 0.4884_pr]
   n = [0.3_pr, 0.45_pr, 0.25_pr]
   T = 330._pr
   V = 0.6_pr

   base_model = PengRobinson76(Tc, Pc, w)
   call symmetric_tensors(k, l)

   mixrule = CMR(k=k, l=l)
   zero_rule = CMR(k=0*k, l=0*l)

   call test_zero_interactions()
   call test_pure_component_limit()
   call test_d_and_b_derivatives()
   call test_d1_mixing()
   call test_full_model_derivatives()

contains

   subroutine symmetric_tensors(k, l)
      !! Fully symmetric \(k_{ijk}\), \(l_{ijk}\) built from distinct values.
      real(pr), intent(out) :: k(nc, nc, nc)
      real(pr), intent(out) :: l(nc, nc, nc)
      integer :: i
      integer :: j
      integer :: m
      integer :: s(3)
      do i = 1, nc
         do j = 1, nc
            do m = 1, nc
               s = [i, j, m]
               k(i, j, m) = 0.01_pr*(sum(s) + product(s)) - 0.05_pr
               l(i, j, m) = 0.004_pr*(sum(s**2)) - 0.02_pr
            end do
         end do
      end do
   end subroutine symmetric_tensors

   subroutine test_zero_interactions()
      real(pr) :: B
      real(pr) :: dBi(nc)
      real(pr) :: dBij(nc, nc)

      call zero_rule%Bmix(n, base_model%b, B, dBi, dBij)
      call assert(abs(B - sum(n*base_model%b)) < 1e-13_pr, "CMR Bmix: l=0 gives the linear B = sum(n_i b_i)")
      call assert(maxval(abs(dBi - base_model%b)) < 1e-13_pr, "CMR Bmix: l=0 gives dB/dn_i = b_i")
      call assert(maxval(abs(dBij)) < 1e-12_pr, "CMR Bmix: l=0 gives zero second derivatives")
   end subroutine test_zero_interactions

   subroutine test_pure_component_limit()
      !! Pure component 2: the CMR with interactions must behave as a pure
      !! fluid with \(a_c (1-k_{222})\) and \(b (1-l_{222})\).
      real(pr) :: npure(nc)
      real(pr) :: Ar_cmr
      real(pr) :: Ar_ref
      real(pr) :: ArV_cmr
      real(pr) :: ArV_ref
      real(pr) :: ArT_cmr
      real(pr) :: ArT_ref

      npure = [0._pr, 1.3_pr, 0._pr]

      model = base_model
      call model%set_mixrule(mixrule)

      scaled = base_model
      scaled%ac(2) = base_model%ac(2)*(1 - k(2, 2, 2))
      scaled%b(2) = base_model%b(2)*(1 - l(2, 2, 2))

      call model%residual_helmholtz(npure, V, T, Ar=Ar_cmr, ArV=ArV_cmr, ArT=ArT_cmr)
      call scaled%residual_helmholtz(npure, V, T, Ar=Ar_ref, ArV=ArV_ref, ArT=ArT_ref)

      call assert(abs(Ar_cmr - Ar_ref) < 1e-10_pr*abs(Ar_ref), &
         "CMR pure limit: Ar equals the pure fluid with a(1-k_iii), b(1-l_iii)")
      call assert(abs(ArV_cmr - ArV_ref) < 1e-10_pr*abs(ArV_ref), "CMR pure limit: same ArV")
      call assert(abs(ArT_cmr - ArT_ref) < 1e-10_pr*abs(ArT_ref), "CMR pure limit: same ArT")
   end subroutine test_pure_component_limit

   subroutine test_d_and_b_derivatives()
      real(pr) :: ai(nc)
      real(pr) :: daidt(nc)
      real(pr) :: daidt2(nc)
      real(pr) :: D
      real(pr) :: dDdV
      real(pr) :: dDdT
      real(pr) :: dDdV2
      real(pr) :: dDdT2
      real(pr) :: dDi(nc)
      real(pr) :: dDdTV
      real(pr) :: dDidV(nc)
      real(pr) :: dDidT(nc)
      real(pr) :: dDij(nc, nc)
      real(pr) :: B
      real(pr) :: dBi(nc)
      real(pr) :: dBij(nc, nc)

      real(pr) :: Dp
      real(pr) :: Dm
      real(pr) :: dDi_p(nc)
      real(pr) :: dDi_m(nc)
      real(pr) :: dDT_p
      real(pr) :: dDT_m
      real(pr) :: Bp
      real(pr) :: Bm
      real(pr) :: dBi_p(nc)
      real(pr) :: dBi_m(nc)
      real(pr) :: dn(nc)
      real(pr) :: eps
      real(pr) :: dT
      real(pr) :: num
      integer :: i
      integer :: j

      call base_model%calculate_attractive_parameters(T, ai, daidt, daidt2)
      call mixrule%Dmix(n, V, T, ai, daidt, daidt2, D, dDdV, dDdT, dDdV2, dDdT2, &
         dDi, dDdTV, dDidV, dDidT, dDij)
      call mixrule%Bmix(n, base_model%b, B, dBi, dBij)

      ! Euler homogeneity
      call assert(abs(sum(n*dDi) - 2*D) < 1e-10_pr*abs(D), "CMR Dmix: sum(n_i dD/dn_i) = 2D")
      call assert(abs(sum(n*dBi) - B) < 1e-12_pr*abs(B), "CMR Bmix: sum(n_i dB/dn_i) = B")
      call assert(maxval(abs(matmul(dDij, n) - dDi)) < 1e-9_pr*maxval(abs(dDi)), &
         "CMR Dmix: sum_j n_j d2D/dn_i dn_j = dD/dn_i")
      call assert(maxval(abs(matmul(dBij, n))) < 1e-12_pr, "CMR Bmix: sum_j n_j d2B/dn_i dn_j = 0")
      ! dDidV is not checked: CMR_Dmix leaves it undefined (it is unused by
      ! the non-DDLC cubic Ar routine).
      call assert(abs(dDdV) == 0 .and. abs(dDdV2) == 0 .and. abs(dDdTV) == 0, &
         "CMR Dmix: no volume dependence")

      ! Composition derivatives by central differences
      eps = 1e-5_pr
      do i = 1, nc
         dn = 0
         dn(i) = eps
         call mixrule%Dmix(n + dn, V, T, ai, daidt, daidt2, Dp, dDdV, dDT_p, dDdV2, dDdT2, &
            dDi_p, dDdTV, dDidV, dDidT, dDij)
         call mixrule%Dmix(n - dn, V, T, ai, daidt, daidt2, Dm, dDdV, dDT_m, dDdV2, dDdT2, &
            dDi_m, dDdTV, dDidV, dDidT, dDij)
         num = (Dp - Dm)/(2*eps)
         call assert(abs(num - dDi(i)) < 1e-7_pr*abs(dDi(i)), "CMR Dmix: dD/dn_i vs finite differences")

         call mixrule%Bmix(n + dn, base_model%b, Bp, dBi_p, dBij)
         call mixrule%Bmix(n - dn, base_model%b, Bm, dBi_m, dBij)
         num = (Bp - Bm)/(2*eps)
         call assert(abs(num - dBi(i)) < 1e-8_pr*abs(dBi(i)), "CMR Bmix: dB/dn_i vs finite differences")
      end do

      ! Second composition derivatives
      call mixrule%Dmix(n, V, T, ai, daidt, daidt2, D, dDdV, dDdT, dDdV2, dDdT2, &
         dDi, dDdTV, dDidV, dDidT, dDij)
      call mixrule%Bmix(n, base_model%b, B, dBi, dBij)
      do j = 1, nc
         dn = 0
         dn(j) = eps
         call mixrule%Dmix(n + dn, V, T, ai, daidt, daidt2, Dp, dDdV, dDT_p, dDdV2, dDdT2, &
            dDi_p, dDdTV, dDidV, dDidT, dDij)
         call mixrule%Dmix(n - dn, V, T, ai, daidt, daidt2, Dm, dDdV, dDT_m, dDdV2, dDdT2, &
            dDi_m, dDdTV, dDidV, dDidT, dDij)
         call mixrule%Bmix(n + dn, base_model%b, Bp, dBi_p, dBij)
         call mixrule%Bmix(n - dn, base_model%b, Bm, dBi_m, dBij)
         call mixrule%Dmix(n, V, T, ai, daidt, daidt2, D, dDdV, dDdT, dDdV2, dDdT2, &
            dDi, dDdTV, dDidV, dDidT, dDij)
         call mixrule%Bmix(n, base_model%b, B, dBi, dBij)
         call assert(maxval(abs((dDi_p - dDi_m)/(2*eps) - dDij(:, j))) < 1e-6_pr*maxval(abs(dDij)), &
            "CMR Dmix: d2D/dn_i dn_j vs finite differences")
         call assert(maxval(abs((dBi_p - dBi_m)/(2*eps) - dBij(:, j))) < 1e-6_pr*maxval(abs(dBi)), &
            "CMR Bmix: d2B/dn_i dn_j vs finite differences")
         call assert(abs((dDT_p - dDT_m)/(2*eps) - dDidT(j)) < 1e-6_pr*abs(dDidT(j)), &
            "CMR Dmix: d2D/dT dn_i vs finite differences")
      end do

      ! Temperature derivatives
      dT = 1e-3_pr
      call base_model%calculate_attractive_parameters(T + dT, ai, daidt, daidt2)
      call mixrule%Dmix(n, V, T + dT, ai, daidt, daidt2, Dp, dDdV, dDT_p, dDdV2, dDdT2, &
         dDi_p, dDdTV, dDidV, dDidT, dDij)
      call base_model%calculate_attractive_parameters(T - dT, ai, daidt, daidt2)
      call mixrule%Dmix(n, V, T - dT, ai, daidt, daidt2, Dm, dDdV, dDT_m, dDdV2, dDdT2, &
         dDi_m, dDdTV, dDidV, dDidT, dDij)
      call base_model%calculate_attractive_parameters(T, ai, daidt, daidt2)
      call mixrule%Dmix(n, V, T, ai, daidt, daidt2, D, dDdV, dDdT, dDdV2, dDdT2, &
         dDi, dDdTV, dDidV, dDidT, dDij)
      call assert(abs((Dp - Dm)/(2*dT) - dDdT) < 1e-6_pr*abs(dDdT), "CMR Dmix: dD/dT vs finite differences")
      call assert(abs((dDT_p - dDT_m)/(2*dT) - dDdT2) < 1e-5_pr*abs(dDdT2), &
         "CMR Dmix: d2D/dT2 vs finite differences")
   end subroutine test_d_and_b_derivatives

   subroutine test_d1_mixing()
      real(pr) :: D1
      real(pr) :: dD1i(nc)
      real(pr) :: dD1ij(nc, nc)

      call mixrule%D1mix(n, base_model%del1, D1, dD1i, dD1ij)
      call assert(abs(D1 - (1 + sqrt(2._pr))) < 1e-13_pr, "CMR D1mix: constant delta1 is kept")
      call assert(maxval(abs(dD1i)) < 1e-13_pr, "CMR D1mix: constant delta1 has zero dD1/dn")
      call assert(maxval(abs(dD1ij)) < 1e-13_pr, "CMR D1mix: constant delta1 has zero d2D1/dn2")
   end subroutine test_d1_mixing

   subroutine test_full_model_derivatives()
      real(pr) :: Ar
      real(pr) :: ArV
      real(pr) :: ArT
      real(pr) :: Arn(nc)
      real(pr) :: Arp
      real(pr) :: Arm
      real(pr) :: dn(nc)
      real(pr) :: eps
      real(pr) :: Ar_qmr
      integer :: i

      model = base_model
      call model%set_mixrule(mixrule)

      call model%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn)
      call base_model%residual_helmholtz(n, V, T, Ar=Ar_qmr)
      call assert(abs(Ar - Ar_qmr) > 1e-4_pr*abs(Ar_qmr), "CMR model: interactions change Ar")

      eps = 1e-6_pr
      do i = 1, nc
         dn = 0
         dn(i) = eps
         call model%residual_helmholtz(n + dn, V, T, Ar=Arp)
         call model%residual_helmholtz(n - dn, V, T, Ar=Arm)
         call assert(abs((Arp - Arm)/(2*eps) - Arn(i)) < 1e-6_pr*abs(Arn(i)), &
            "CMR model: Arn vs finite differences")
      end do

      call model%residual_helmholtz(n, V + 1e-6_pr, T, Ar=Arp)
      call model%residual_helmholtz(n, V - 1e-6_pr, T, Ar=Arm)
      call assert(abs((Arp - Arm)/2e-6_pr - ArV) < 1e-6_pr*abs(ArV), "CMR model: ArV vs finite differences")

      call model%residual_helmholtz(n, V, T + 1e-3_pr, Ar=Arp)
      call model%residual_helmholtz(n, V, T - 1e-3_pr, Ar=Arm)
      call assert(abs((Arp - Arm)/2e-3_pr - ArT) < 1e-6_pr*abs(ArT), "CMR model: ArT vs finite differences")
   end subroutine test_full_model_derivatives
end program test_cmr_cubic_eos
