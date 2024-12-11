program main
   use yaeos
   use yaeos__math, only: interpol
   use forsus, only: forsus_default_dir, forsus_dir, Substance
   use testing_aux, only: test_title, assert

   type(EquilibriumState) :: sat
   type(TXEnvel2) :: tx
   type(CubicEoS) :: model
   type(CriticalLine) :: cl, cl_ll
   type(EquilibriumState) :: cp_ll

   integer, parameter :: nc=2
   real(pr) :: z0(nc)= [1, 0]
   real(pr) :: zi(nc)= [0, 1]
   real(pr) :: z(nc), P
   real(pr) :: a
   real(pr) :: Tenv
   type(Substance) :: sus(nc)

   integer :: i, tf_brute, tf_criti, tf_envel

   open(newunit=tf_brute, file="test_tx_brute")
   open(newunit=tf_criti, file="test_tx_criti")
   open(newunit=tf_envel, file="test_tx_envel")


   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("methane")
   sus(2) = Substance("propane")

   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )

   p = 10

   a = 0.01
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="bubble", T0=150._pr)
   tx = tx_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat)

   i = minloc(abs(270._pr - tx%points%T), dim=1)

   Tenv = interpol(&
      tx%points(i)%y(1), tx%points(i+1)%y(1), &
      tx%points(i)%T, tx%points(i+1)%T, 0.53_pr &
      )

   call assert(abs(Tenv - 270.3_pr) < 0.1_pr, "Predicted Dew Temperature")
end program main
