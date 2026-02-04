program test_pcsaft
   use yaeos
   use yaeos__models_ar_saft_pcsaft
   use testing_aux, only: test_title, assert
   use auxiliar_functions, only: allclose
   implicit none

   write(*, *)  test_title("PC-SAFT EoS")

   call pt_methane_hexane

contains
   subroutine pt_methane_hexane
      integer, parameter :: nc = 2
      real(pr) :: Tc(2), Pc(2), w(2)
      real(pr) :: m(2), sigma(2), epsilon_k(2)
      type(PcSaft) :: model
      type(EquilibriumState) :: sat_point

      real(pr) :: P_test = 95.087
      real(pr) :: y_test(2) = [0.98185, 1-0.98185]

      real(pr) :: z(2)

      Tc = [190.564_pr, 507.6_pr]
      Pc = [45.99_pr, 30.25_pr]
      w = [0.01155, 0.30126]

      m = [1.0582, 3.3004]
      epsilon_k = [145.5257, 224.0780]
      sigma = [3.6316, 3.8639]
      model = init_pcsaft(m, sigma, epsilon_k)

      z = [0.5, 0.5]
      sat_point = saturation_pressure(model, z, T=300._pr, kind="bubble", P0=100._pr)
      call assert(abs(sat_point%P - P_test) < 1.0e-3_pr, &
         "PC-SAFT Methane/Hexane bubble point pressure test failed.")
      call assert(allclose(sat_point%y, y_test, rtol=1.0e-2_pr), &
         "PC-SAFT Methane/Hexane bubble point composition test failed.")
   end subroutine pt_methane_hexane
end program test_pcsaft
