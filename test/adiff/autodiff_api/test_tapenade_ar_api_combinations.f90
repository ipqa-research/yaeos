program test_tapenade_ar_api_combinations
   !! `ArModelTapenade%residual_helmholtz` with subsets of the optional
   !! outputs. The routine reuses byproducts of each automatic
   !! differentiation sweep; requesting \(\partial^2 A_r/\partial T\partial n\)
   !! without \(\partial^2 A_r/\partial n^2\) or \(\partial^2 A_r/\partial V\partial n\)
   !! must still produce the same values as the analytical Peng-Robinson model.
   use yaeos, only: pr, CubicEoS, PengRobinson76
   use yaeos__models_ar_cubic_alphas, only: AlphaSoave
   use autodiff_tapenade_pr76, only: TPR76, setup_model
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 2
   type(TPR76) :: tape
   type(CubicEoS) :: analytic
   type(AlphaSoave) :: alpha

   real(pr) :: Tc(nc)
   real(pr) :: Pc(nc)
   real(pr) :: w(nc)
   real(pr) :: kij(nc, nc)
   real(pr) :: lij(nc, nc)
   real(pr) :: n(nc)
   real(pr) :: V
   real(pr) :: T

   real(pr) :: Ar
   real(pr) :: ArV
   real(pr) :: ArT
   real(pr) :: ArTV
   real(pr) :: ArV2
   real(pr) :: ArT2
   real(pr) :: Arn(nc)
   real(pr) :: ArVn(nc)
   real(pr) :: ArTn(nc)
   real(pr) :: Arn2(nc, nc)

   real(pr) :: t_ArV
   real(pr) :: t_ArT
   real(pr) :: t_ArTV
   real(pr) :: t_ArT2
   real(pr) :: t_ArV2
   real(pr) :: t_Arn(nc)
   real(pr) :: t_ArTn(nc)
   real(pr) :: t_ArVn(nc)

   real(pr) :: f_ArV
   real(pr) :: f_ArT
   real(pr) :: f_ArTV
   real(pr) :: f_ArT2
   real(pr) :: f_Arn(nc)
   real(pr) :: f_ArTn(nc)
   real(pr) :: f_ArVn(nc)
   real(pr) :: f_ArV2
   real(pr) :: f_Arn2(nc, nc)

   real(pr), parameter :: rtol = 1e-6_pr
   real(pr), parameter :: self_tol = 1e-12_pr

   write(*, *) test_title("Tapenade Ar API: partial requests of derivatives")

   Tc = [190.564_pr, 425.12_pr]
   Pc = [45.99_pr, 37.96_pr]
   w = [0.0115478_pr, 0.200164_pr]
   kij = reshape([0._pr, 0.1_pr, 0.1_pr, 0._pr], [nc, nc])
   lij = reshape([0._pr, -0.02_pr, -0.02_pr, 0._pr], [nc, nc])

   tape = setup_model(Tc, Pc, w, kij, lij)
   analytic = PengRobinson76(Tc, Pc, w, kij, lij)
   ! The fixture uses default-real literals for some constants: copy its
   ! parameters so that both models are as close as possible (the remaining
   ! single-precision constants inside the fixture limit the agreement to
   ! ~1e-8, hence the tolerance).
   analytic%ac = tape%ac
   analytic%b = tape%b
   analytic%del1 = [tape%del1, tape%del1]
   analytic%del2 = [tape%del2, tape%del2]
   alpha%k = tape%k
   call analytic%set_alpha(alpha)

   n = [0.3_pr, 0.7_pr]
   V = 0.25_pr
   T = 250._pr

   call analytic%residual_helmholtz(n, V, T, Ar=Ar, ArV=ArV, ArT=ArT, ArTV=ArTV, &
      ArV2=ArV2, ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

   ! --------------------------------------------------------------------------
   ! ArTn sweep provides ArT2, ArTV and the gradient
   ! --------------------------------------------------------------------------
   call tape%residual_helmholtz(n, V, T, ArT=t_ArT, ArV=t_ArV, Arn=t_Arn, &
      ArTV=t_ArTV, ArT2=t_ArT2, ArTn=t_ArTn)
   call assert(maxval(abs(t_ArTn - ArTn)) < rtol*maxval(abs(ArTn)), "ArTn sweep: ArTn matches analytic PR76")
   call assert(abs(t_ArT2 - ArT2) < rtol*abs(ArT2), "ArTn sweep: ArT2 byproduct matches analytic PR76")
   call assert(abs(t_ArTV - ArTV) < rtol*abs(ArTV), "ArTn sweep: ArTV byproduct matches analytic PR76")
   call assert(abs(t_ArT - ArT) < rtol*abs(ArT), "ArTn sweep: ArT byproduct matches analytic PR76")
   call assert(abs(t_ArV - ArV) < rtol*abs(ArV), "ArTn sweep: ArV byproduct matches analytic PR76")
   call assert(maxval(abs(t_Arn - Arn)) < rtol*maxval(abs(Arn)), "ArTn sweep: Arn byproduct matches analytic PR76")

   ! Same values as a request of every derivative (Arn2 sweep path)
   call tape%residual_helmholtz(n, V, T, ArV=f_ArV, ArT=f_ArT, ArTV=f_ArTV, &
      ArV2=f_ArV2, ArT2=f_ArT2, Arn=f_Arn, ArVn=f_ArVn, ArTn=f_ArTn, Arn2=f_Arn2)
   call assert(maxval(abs(t_ArTn - f_ArTn)) < self_tol*maxval(abs(f_ArTn)), &
      "ArTn sweep: ArTn equals the full-request value")
   call assert(abs(t_ArT2 - f_ArT2) < self_tol*abs(f_ArT2), "ArTn sweep: ArT2 equals the full-request value")
   call assert(abs(t_ArTV - f_ArTV) < self_tol*abs(f_ArTV), "ArTn sweep: ArTV equals the full-request value")
   call assert(abs(t_ArV - f_ArV) < self_tol*abs(f_ArV), "ArTn sweep: ArV equals the full-request value")
   call assert(abs(t_ArT - f_ArT) < self_tol*abs(f_ArT), "ArTn sweep: ArT equals the full-request value")
   call assert(maxval(abs(t_Arn - f_Arn)) < self_tol*maxval(abs(f_Arn)), &
      "ArTn sweep: Arn equals the full-request value")

   ! Only the gradient wrt n (no ArT/ArV) together with ArTn
   t_Arn = 0
   call tape%residual_helmholtz(n, V, T, Arn=t_Arn, ArTn=t_ArTn)
   call assert(maxval(abs(t_Arn - Arn)) < rtol*maxval(abs(Arn)), "ArTn sweep: Arn alone matches")
   call assert(maxval(abs(t_ArTn - ArTn)) < rtol*maxval(abs(ArTn)), "ArTn sweep: ArTn alone matches")

   ! --------------------------------------------------------------------------
   ! ArVn sweep provides ArV2 and ArTV; ArT2 is left to the ArTn sweep
   ! --------------------------------------------------------------------------
   call tape%residual_helmholtz(n, V, T, ArVn=t_ArVn, ArV2=t_ArV2, ArTV=t_ArTV, &
      ArT2=t_ArT2, ArTn=t_ArTn)
   call assert(maxval(abs(t_ArVn - ArVn)) < rtol*maxval(abs(ArVn)), "ArVn+ArTn sweeps: ArVn matches")
   call assert(abs(t_ArV2 - ArV2) < rtol*abs(ArV2), "ArVn+ArTn sweeps: ArV2 matches")
   call assert(abs(t_ArTV - ArTV) < rtol*abs(ArTV), "ArVn+ArTn sweeps: ArTV matches")
   call assert(abs(t_ArT2 - ArT2) < rtol*abs(ArT2), "ArVn+ArTn sweeps: ArT2 matches")
   call assert(maxval(abs(t_ArTn - ArTn)) < rtol*maxval(abs(ArTn)), "ArVn+ArTn sweeps: ArTn matches")
end program test_tapenade_ar_api_combinations
