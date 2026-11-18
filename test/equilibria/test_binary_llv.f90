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

   call water_ethanol()
   call methane_co2()

contains

   subroutine water_ethanol
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
         ns0=spec_CP%P, S0=log(P), ds0=-0.01_pr, v0=V, T0=T, p0=P, stability_analysis=.true.)
      llv = binary_llv_from_cep(model, cl%cep)

      points = size(llv%P)

      ! Check that the first point of the LLV matches the CEP
      call assert(llv%P(points) < 1.e-2_pr, "LLV ends at low pressure")
      ! print *, "LLV first point P:", llv%P(1), "CEP P:", cl%cep%P
      ! call assert(abs(llv%P(1)- cl%cep%P) &
      !       < 1.e-2_pr, "LLV first point pressure matches CEP")

   end subroutine water_ethanol

   subroutine methane_co2
      real(pr) :: kij(nc, nc)
      integer :: i
      Tc = [190.6_pr, 304.1_pr]
      Pc = [46.0_pr, 73.8_pr]
      w = [0.011_pr, 0.225_pr]
      model = PengRobinson76(Tc, Pc, w)

      kij = 0
      kij(1, 2) = 0.01
      kij(2, 1) = 0.01
      z0 = [0, 1]
      zi = [1, 0]

      P = 800
      T = 500

      call find_llcl(model, z0, zi, P, a, V, T)
      cl = critical_line(&
         model, a, z0, zi, &
         ns0=spec_CP%P, S0=log(P), ds0=-0.1_pr, v0=V, T0=T, p0=P, stability_analysis=.true.)
      llv = binary_llv_from_cep(model, cl%cep)

      points = size(llv%P)
      call assert(llv%P(points) < 1.e-2_pr, "LLV ends at low pressure")
   end subroutine methane_co2
end program main
