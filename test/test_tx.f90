program main
   use yaeos
   use yaeos__math, only: interpol
   use yaeos__equilibria_boundaries_pure_saturation, only: PurePsat, pure_saturation_line
   use forsus, only: forsus_default_dir, forsus_dir, Substance
   use testing_aux, only: test_title, assert

   type(EquilibriumState) :: sat
   type(TXEnvel2) :: tx
   type(CubicEoS) :: model

   integer, parameter :: nc=2
   real(pr) :: z0(nc)= [1, 0]
   real(pr) :: zi(nc)= [0, 1]
   real(pr) :: z(nc), T, P
   real(pr) :: a
   real(pr) :: Tenv
   type(Substance) :: sus(nc)
   type(PurePsat) :: vp1, vp2

   integer :: i

   write(*, *) test_title("TX Envelope 2ph")

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("methane")
   sus(2) = Substance("propane")

   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )

   vp1 = pure_saturation_line(model, 1, 1._pr, 100._pr)
   vp2 = pure_saturation_line(model, 2, 1._pr, 100._pr)

   P = 10._pr
   T = vp1%get_T(P)

   a = 0.001
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="bubble", T0=T)
   tx = tx_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat)

   i = minloc(abs(270._pr - tx%points%T), dim=1)

   Tenv = interpol(&
      tx%points(i)%y(1), tx%points(i+1)%y(1), &
      tx%points(i)%T, tx%points(i+1)%T, 0.53_pr &
   )

   call assert(abs(Tenv - 270.3_pr) < 0.1_pr, "Predicted Dew Temperature")

   P = sum(model%components%Pc)/2
   T = vp1%get_T(P)
   sat = saturation_temperature(model, z, P=P, kind="bubble", T0=T)
   tx = tx_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat)

   call assert(abs(tx%points(size(tx%points))%T - 368) < 1, "Reach to critical point")
end program main
