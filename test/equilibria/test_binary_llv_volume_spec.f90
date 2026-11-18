program test_binary_llv_volume_spec
   !! Binary LLV line equations specified by the liquid volume ratio
   !! (`ns = -1`, \(F_{10} = \ln(V_x/V_y) - S\)).
   !!
   !! A converged LLV point from `binary_llv_from_cep` must satisfy the
   !! equations with that specification, the analytical Jacobian must agree
   !! with finite differences, and solving with a perturbed volume ratio must
   !! give another LLV point with the requested ratio.
   use yaeos, only: pr, CubicEoS, SoaveRedlichKwong, CriticalLine, critical_line, &
      spec_CP, find_llcl
   use yaeos__equilibria_binaries, only: BinaryThreePhase, binary_llv_from_cep, &
      three_phase_line_F, three_phase_line_F_solve
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 2
   integer, parameter :: nv = 10
   type(CubicEoS) :: model
   type(CriticalLine) :: cl
   type(BinaryThreePhase) :: llv
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)
   real(pr) :: P
   real(pr) :: T
   real(pr) :: V
   real(pr) :: a
   real(pr) :: X(nv)
   real(pr) :: Xp(nv)
   real(pr) :: Xm(nv)
   real(pr) :: F(nv)
   real(pr) :: Fp(nv)
   real(pr) :: Fm(nv)
   real(pr) :: dF(nv, nv)
   real(pr) :: dFtmp(nv, nv)
   real(pr) :: num(nv)
   real(pr) :: S
   real(pr) :: S_new
   real(pr) :: h
   real(pr) :: Px
   real(pr) :: Py
   real(pr) :: Pw
   integer :: i
   integer :: j
   integer :: k
   integer :: iters

   write(*, *) test_title("Binary LLV equations with a volume-ratio specification")

   ! Water + ethanol with SRK (same system as test_binary_llv)
   model = SoaveRedlichKwong([647._pr, 514._pr], [218._pr, 63._pr], [0.344_pr, 0.644_pr])
   z0 = [0._pr, 1._pr]
   zi = [1._pr, 0._pr]
   P = 2000
   T = 500
   call find_llcl(model, z0, zi, P, a, V, T)
   cl = critical_line(model, a, z0, zi, ns0=spec_CP%P, S0=log(P), ds0=-0.01_pr, &
      v0=V, T0=T, p0=P, stability_analysis=.true.)
   llv = binary_llv_from_cep(model, cl%cep)
   call assert(size(llv%T) > 10, "LLV line computed")

   ! A point in the middle of the line
   k = size(llv%T)/2
   X = [ &
      log(llv%x1(k)), log(1 - llv%x1(k)), &
      log(llv%y1(k)), log(1 - llv%y1(k)), &
      log(llv%w1(k)), log(1 - llv%w1(k)), &
      log(llv%Vx(k)), log(llv%Vy(k)), log(llv%Vw(k)), log(llv%T(k))]
   S = log(llv%Vx(k)/llv%Vy(k))

   call three_phase_line_F(model, X, -1, S, F, dF)
   call assert(abs(F(10)) < 1e-12_pr, "ns=-1: F(10) = ln(Vx/Vy) - S vanishes on the point")
   call assert(maxval(abs(F)) < 1e-7_pr, "ns=-1: stored LLV point satisfies the equations")

   ! Jacobian (wrt the log variables) against central differences
   do j = 1, nv
      h = 1e-6_pr
      Xp = X
      Xm = X
      Xp(j) = X(j) + h
      Xm(j) = X(j) - h
      call three_phase_line_F(model, Xp, -1, S, Fp, dFtmp)
      call three_phase_line_F(model, Xm, -1, S, Fm, dFtmp)
      num = (Fp - Fm)/(2*h)
      do i = 1, nv
         call assert(abs(num(i) - dF(i, j)) < 1e-5_pr*max(1._pr, abs(dF(i, j))), &
            "ns=-1: analytical Jacobian matches finite differences")
      end do
   end do
   call assert(abs(dF(10, 7) - 1) < 1e-12_pr .and. abs(dF(10, 8) + 1) < 1e-12_pr, &
      "ns=-1: d F10/d lnVx = 1 and d F10/d lnVy = -1")

   ! Solve a new point with a slightly different liquid-volume ratio
   S_new = S + 0.02_pr*sign(1._pr, S)
   call three_phase_line_F_solve(model, X, -1, S_new, F, dF, iters)
   call assert(maxval(abs(F)) < 1e-9_pr, "ns=-1: new LLV point converged")
   call assert(abs(X(7) - X(8) - S_new) < 1e-9_pr, "ns=-1: new point has the requested volume ratio")

   ! The three phases are at the same pressure
   call model%pressure(exp(X(1:2)), exp(X(7)), exp(X(10)), Px)
   call model%pressure(exp(X(3:4)), exp(X(8)), exp(X(10)), Py)
   call model%pressure(exp(X(5:6)), exp(X(9)), exp(X(10)), Pw)
   call assert(abs(Px - Pw) < 1e-6_pr*abs(Pw) + 1e-8_pr .and. abs(Py - Pw) < 1e-6_pr*abs(Pw) + 1e-8_pr, &
      "ns=-1: mechanical equilibrium of the three phases")
   call assert(abs(sum(exp(X(1:2))) - 1) < 1e-9_pr .and. abs(sum(exp(X(5:6))) - 1) < 1e-9_pr, &
      "ns=-1: phase compositions are normalized")
end program test_binary_llv_volume_spec
