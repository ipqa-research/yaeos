program test_generic_cubic_setters
   !! `CubicEoS%set_delta1` and `CubicEoS%set_alpha`.
   !!
   !! A Peng-Robinson model whose \(\delta_1\) is set to 1 must be the same
   !! model as an SRK model whose \(a_c\), \(b\) and alpha function are replaced
   !! by the Peng-Robinson ones. Replacing the alpha function must change the
   !! attractive parameters exactly as the new Soave \(k\) dictates.
   use yaeos, only: pr, CubicEoS, PengRobinson76, SoaveRedlichKwong
   use yaeos__models_ar_cubic_alphas, only: AlphaSoave
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(CubicEoS) :: pr_model
   type(CubicEoS) :: model
   type(CubicEoS) :: srk
   type(AlphaSoave) :: alpha

   real(pr) :: Tc(nc)
   real(pr) :: Pc(nc)
   real(pr) :: w(nc)
   real(pr) :: n(nc)
   real(pr) :: V
   real(pr) :: T
   real(pr) :: Ar_a
   real(pr) :: Ar_b
   real(pr) :: ArV_a
   real(pr) :: ArV_b
   real(pr) :: ArT_a
   real(pr) :: ArT_b
   real(pr) :: Arn_a(nc)
   real(pr) :: Arn_b(nc)
   real(pr) :: a(nc)
   real(pr) :: dadt(nc)
   real(pr) :: dadt2(nc)
   real(pr) :: a_ref(nc)
   real(pr) :: k_new(nc)
   real(pr) :: Tr(nc)
   real(pr) :: sq(nc)

   write(*, *) test_title("CubicEoS set_delta1 and set_alpha")

   Tc = [190.564_pr, 369.83_pr, 617.7_pr]
   Pc = [45.99_pr, 42.48_pr, 21.1_pr]
   w = [0.0115478_pr, 0.1523_pr, 0.4884_pr]
   n = [0.4_pr, 0.35_pr, 0.25_pr]
   V = 0.5_pr
   T = 350._pr

   pr_model = PengRobinson76(Tc, Pc, w)

   ! ==========================================================================
   ! set_delta1 with the model's own delta1 is a no-op
   ! --------------------------------------------------------------------------
   model = pr_model
   call model%set_delta1(pr_model%del1)
   call assert(maxval(abs(model%del2 - pr_model%del2)) < 1e-14_pr, &
      "set_delta1: PR delta1 gives back the PR delta2 = 1 - sqrt(2)")
   call assert(maxval(abs(model%del2 - (1 - sqrt(2._pr)))) < 1e-14_pr, &
      "set_delta1: delta2 = (1-delta1)/(1+delta1)")

   ! ==========================================================================
   ! PR with delta1=1 is the SRK volume function with PR parameters
   ! --------------------------------------------------------------------------
   call model%set_delta1([1._pr, 1._pr, 1._pr])
   call assert(all(model%del1 == 1._pr), "set_delta1: delta1 stored")
   call assert(maxval(abs(model%del2)) < 1e-15_pr, "set_delta1: delta1=1 gives delta2=0")

   srk = SoaveRedlichKwong(Tc, Pc, w)
   srk%ac = pr_model%ac
   srk%b = pr_model%b
   call srk%set_alpha(pr_model%alpha)

   call model%residual_helmholtz(n, V, T, Ar=Ar_a, ArV=ArV_a, ArT=ArT_a, Arn=Arn_a)
   call srk%residual_helmholtz(n, V, T, Ar=Ar_b, ArV=ArV_b, ArT=ArT_b, Arn=Arn_b)
   call assert(abs(Ar_a - Ar_b) < 1e-10_pr*abs(Ar_b), "set_delta1/set_alpha: PR(delta1=1) Ar == SRK(PR params) Ar")
   call assert(abs(ArV_a - ArV_b) < 1e-10_pr*abs(ArV_b), "set_delta1/set_alpha: same ArV")
   call assert(abs(ArT_a - ArT_b) < 1e-10_pr*abs(ArT_b), "set_delta1/set_alpha: same ArT")
   call assert(maxval(abs(Arn_a - Arn_b)) < 1e-10_pr*maxval(abs(Arn_b)), "set_delta1/set_alpha: same Arn")

   ! And it must differ from the true PR model
   call pr_model%residual_helmholtz(n, V, T, Ar=Ar_b)
   call assert(abs(Ar_a - Ar_b) > 1e-3_pr*abs(Ar_b), "set_delta1: changing delta1 changes Ar")

   ! ==========================================================================
   ! set_alpha replaces the alpha function
   ! --------------------------------------------------------------------------
   k_new = [0.2_pr, 0.5_pr, 0.9_pr]
   alpha%k = k_new
   model = pr_model
   call model%set_alpha(alpha)
   call model%calculate_attractive_parameters(T, a, dadt, dadt2)

   Tr = T/Tc
   sq = sqrt(Tr)
   a_ref = pr_model%ac*(1 + k_new*(1 - sq))**2
   call assert(maxval(abs(a - a_ref)/a_ref) < 1e-12_pr, "set_alpha: a_i = ac_i*(1+k(1-sqrt(Tr)))^2 with the new k")
   ! da/dT = ac * 2*(1+k(1-sqrt Tr)) * (-k/(2 sqrt Tr)) / Tc
   call assert(maxval(abs(dadt - pr_model%ac*(1 + k_new*(1 - sq))*(-k_new/sq)/Tc)/abs(dadt)) < 1e-12_pr, &
      "set_alpha: analytic da/dT of the new alpha")

   ! At T = Tc_i the alpha function is one, a_i = ac_i, whatever k is
   call model%calculate_attractive_parameters(Tc(2), a, dadt, dadt2)
   call assert(abs(a(2) - pr_model%ac(2)) < 1e-12_pr*pr_model%ac(2), "set_alpha: alpha(Tr=1) = 1")

   ! Setting the alpha twice keeps the last one (deallocation path)
   alpha%k = 2*k_new
   call model%set_alpha(alpha)
   call model%calculate_attractive_parameters(T, a, dadt, dadt2)
   a_ref = pr_model%ac*(1 + 2*k_new*(1 - sq))**2
   call assert(maxval(abs(a - a_ref)/a_ref) < 1e-12_pr, "set_alpha: second call replaces the alpha function")
end program test_generic_cubic_setters
