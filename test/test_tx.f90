program main
   use yaeos
   use yaeos__equilibria_boundaries_phase_envelopes_tx, only: TXEnvel2, tx_envelope_2ph
   use forsus, only: forsus_default_dir, forsus_dir, Substance

   type(EquilibriumState) :: sat
   type(TXEnvel2) :: tx
   type(CubicEoS) :: model
   type(CriticalLine) :: cl, cl_ll
   type(EquilibriumState) :: cp_ll

   integer, parameter :: nc=2
   real(pr) :: zi(nc)= [0, 1]
   real(pr) :: z0(nc)= [1, 0]
   real(pr) :: z(nc)
   real(pr) :: a
   real(pr) :: P, maxPc
   type(Substance) :: sus(nc)
   real(pr) :: lnphix(nc), lnphiy(nc)

   integer :: i


   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("methane")
   sus(2) = Substance("propane")

   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )

   p = 10

   print *, "Pc1=", model%components%Pc(1)
   print *, "Pc2=", model%components%Pc(2)
   print *, "Tc1=", model%components%Tc(1)
   print *, "Tc2=", model%components%Tc(2)
   print *, "P=", P
   maxPc = maxval(model%components%Pc)


   do i=1,99
      a = real(i)/100
      z = z0*a + zi*(1-a)
      sat = saturation_temperature(model, z, P=P, kind="dew", T0=250._pr)
      !tx = tx_envelope_2ph(model, z, sat)
      write (1, *) sat%x(2), sat%y(2), sat%T
   end do

   a = 0.001
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="dew", T0=150._pr)
   tx = tx_envelope_2ph(model, z0=z0, z_injection=zi, alpha0=a, first_point=sat)

   print *, sat%iters, size(tx%points)
   do i=1,size(tx%points)
      write(2, *) tx%points(i)%x(2), tx%points(i)%y(2), tx%points(i)%T
   end do

   write(2, *)
   write(2, *)

   cp_ll = critical_point(model, z0=z0, zi=zi, spec=spec_CP%P, S=log(maxPc+10), max_iters=1500, a0=0.5_pr)
   cl_ll = critical_line(&
      model, z0=z0, zi=zi, a0=cp_ll%x(2), &
      ns=spec_CP%P, S=log(cp_ll%P), dS0=-0.05_pr, first_point=cp_ll)
   print *, cp_ll%iters, size(cl_ll%a), cl_ll%P(1), cl_ll%P(size(cl_ll%a))

   a = cp_ll%x(2)-0.01_pr
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="dew", T0=cp_ll%T)
   print *, "from X2", sat%iters, sat
   tx = tx_envelope_2ph(&
      model, z0=z0, z_injection=zi, alpha0=a, &
      first_point=sat, delta_0=0.01_pr, specified_variable_0=3&
      )

   print *, sat%iters, size(tx%points)
   do i=1,size(tx%points)
      write(2, *) tx%points(i)%x(2), tx%points(i)%y(2), tx%points(i)%T
   end do

   cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, ns=1, S=0.99_pr, dS0=-0.01_pr)
   print *, size(cl%a)

   do i=1,size(cl%a)
      write(3, *) cl%a(i), cl%T(i), cl%P(i)
   end do

end program main
