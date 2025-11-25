program main
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, saturation_temperature, HV_NRTL
   use fixtures_models, only: binary_PR76
   implicit none
   type(CubicEoS) :: model
   type(EquilibriumState) :: sat
   type(HV_NRTL) :: mr
   logical :: use_kij(3, 3)
   real(pr) :: z(3), gji(3,3), alpha(3,3), kij(3,3), Tc(3), Pc(3), w(3), P
   real(pr) :: y0(3)

   gji = 0
   gji(1, 2) = 20.53672267_pr
   gji(1, 3) = 471.015777_pr
   gji(2, 1) = 246.93953975_pr
   gji(2, 3) = -50.0_pr
   gji(3, 1) = 36.584873_pr
   gji(3, 2) = 37.5912979_pr

   alpha = 0
   alpha(1, 2) = 0.4_pr
   alpha(1, 3) = 0.193694337_pr
   alpha(2, 1) = 0.4_pr
   alpha(2, 3) = 0.05_pr
   alpha(3, 1) = 0.193694337_pr
   alpha(3, 2) = 0.05_pr

   use_kij = .false.
   kij = 0.0_pr

   Tc = [304.21, 512.5 , 720.  ]
   Pc = [73.83, 80.84, 82.  ]
   w = [0.223621, 0.565831, 0.506776]

   model = PengRobinson76(Tc, Pc, w)
   mr = HV_NRTL(b=model%b, del1=model%del1, gji=gji, alpha=alpha, use_kij=use_kij, kij=kij)
   call model%set_mixrule(mr)

   z = [0.95, 0.04, 0.01]
   y0 = 1-z
   P = 1e-1_pr
   do while(sat%kind /= "failed")
      P = 2e-6
      sat = saturation_temperature(model, n=z, P=P, kind="dew", y0=y0, t0=200._pr)
      P = P/10
      print *, sat
      print *, ""
      P = 1e-5
      sat = saturation_temperature(model, n=z, P=P, kind="dew", y0=y0, t0=200._pr)
      print *, sat
      print *, ""
      P = 1e-4
      sat = saturation_temperature(model, n=z, P=P, kind="dew", y0=y0, t0=200._pr)
      print *, sat
      print *, ""
      
      call exit
   end do


end program main
