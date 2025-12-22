program main
   use yaeos
   use yaeos__equilibria_binaries, only: BinaryThreePhase, binary_llv_from_cep
   use testing_aux, only: test_title, assert
   implicit none
   integer, parameter :: nc=2

   real(pr) :: Tc(2), Pc(2), w(2)
   real(pr) :: z0(nc), zi(nc)
   real(pr) :: P, T, V, a
   type(CubicEoS) :: model
   type(CriticalLine) :: cl
   type(BinaryThreePhase) :: llv

   integer :: points

   write(*, *) test_title("Binary LLV from CEP")

   ! This is kind Water/Ethanol
   Tc = [647, 514]
   Pc = [218, 63]
   w = [0.344, 0.644]

   model = SoaveRedlichKwong(Tc, Pc, w)

   z0 = [0, 1]
   zi = [1, 0]

   P = 2000
   T = 500
   call find_llcl(model, z0, zi, P, a, V, T)
   cl = critical_line(&
      model, a, z0, zi, &
      ns0=spec_CP%P, S0=log(P), ds0=-0.1_pr, v0=V, T0=T, p0=P, stability_analysis=.true.)
   llv = binary_llv_from_cep(model, cl%cep)

   points = size(llv%P)

   ! Check that the first point of the LLV matches the CEP
   call assert(abs(llv%P(1)- cl%cep%P) &
         < 1.e-2_pr, "LLV first point pressure matches CEP")
   call assert(llv%P(points) < 1.e-2_pr, "LLV ends at low pressure")
end program