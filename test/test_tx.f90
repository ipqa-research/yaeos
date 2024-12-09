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

   integer :: i, tf_brute, tf_criti, tf_envel

   open(newunit=tf_brute, file="test_tx_brute")
   open(newunit=tf_criti, file="test_tx_criti")
   open(newunit=tf_envel, file="test_tx_envel")


   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("methane")
   sus(2) = Substance("carbon dioxide")

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

   print *, saturation_temperature(model, z0, P=P, kind="dew", T0=250._pr)
   print *, saturation_temperature(model, zi, P=P, kind="dew", T0=250._pr)

   ! ==============================================================
   ! Brute force calculation
   ! --------------------------------------------------------------
   do i=1,99
      a = real(i)/100
      z = z0*a + zi*(1-a)
      sat = saturation_temperature(model, z, P=P, kind="dew", T0=250._pr)
      write (tf_brute, *) sat%x(2), sat%y(2), sat%T
   end do
   call flush(tf_brute)

   ! ==============================================================
   ! TX from X2
   ! --------------------------------------------------------------
   print *, "TX from X2"
   a = 0.999
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="dew", T0=250._pr)
   print *, sat
   tx = tx_envelope_2ph(model, z0=z0, z_injection=zi, alpha0=a, first_point=sat, delta_0=-0.1_pr, specified_variable_0=nc+2)
   print *, sat%iters, size(tx%points)
   do i=1,size(tx%points)
      write(tf_envel, *) tx%points(i)%x(2), tx%points(i)%y(2), tx%points(i)%T
   end do
   call flush(tf_envel)

   ! ==============================================================
   ! Critical lines
   ! --------------------------------------------------------------
   
   ! Critical line
   cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, ns=1, S=0.99_pr, dS0=-0.01_pr)
   print *, "CL: ", size(cl%a)
   do i=1,size(cl%a)
      write(tf_criti, *) cl%a(i), cl%T(i), cl%P(i)
   end do
   write(tf_criti, *)
   write(tf_criti, *)

   ! HP LL
   print *, "C LL"
   maxPc = 500
   cp_ll = critical_point(model, z0=z0, zi=zi, spec=spec_CP%P, S=log(maxPc), max_iters=2500, a0=0.7_pr)
   cl_ll = critical_line(&
      model, z0=z0, zi=zi, a0=cp_ll%x(2), &
      ns=spec_CP%P, S=log(cp_ll%P), dS0=-0.05_pr, first_point=cp_ll)
      print *, cp_ll%iters
      print *, size(cl_ll%a), cl_ll%P(1), cl_ll%P(size(cl_ll%a))
   do i=1,size(cl_ll%a)
      write(tf_criti, *) cl_ll%a(i), cl_ll%T(i), cl_ll%P(i)
   end do
   
   ! ==============================================================
   ! TX From LL line
   ! --------------------------------------------------------------
   print *, "LL from X2", sat%iters, sat
   a = cp_ll%x(2)-0.01_pr
   z = zi*a + z0*(1-a)
   sat = saturation_temperature(model, z, P=P, kind="dew", T0=cp_ll%T)
   tx = tx_envelope_2ph(&
      model, z0=z0, z_injection=zi, alpha0=a, &
      first_point=sat, delta_0=0.01_pr, specified_variable_0=3&
      )
   write(tf_envel, *)
   write(tf_envel, *)
   print *, sat%iters, size(tx%points)
   do i=1,size(tx%points)
      write(tf_envel, *) tx%points(i)%x(2), tx%points(i)%y(2), tx%points(i)%T
   end do
   call flush(tf_envel)
end program main
