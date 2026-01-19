program main
   use testing_aux, only: test_title, assert
   
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, saturation_temperature, HV_NRTL, Groups, PSRK
   use yaeos, only: PTEnvelMP, pt_envelope
   use fixtures_models, only: binary_PR76
   implicit none
   type(CubicEoS) :: model, model_psrk
   type(EquilibriumState) :: sat
   real(pr) :: z(3), Tc(3), Pc(3), w(3), P
   real(pr) :: y0(3)

   type(PTEnvelMP) :: env
   integer, parameter :: np = 1, nc=3
   real(pr) :: x_l0(np, nc), w0(nc), T0, P0, betas0(np)
   character(len=14) :: kinds_x(np), kind_w
   real(pr) :: a, zi(nc), z0(nc)

   type(Groups) :: molecules(3)

   write(*, *) test_title("Convergence of low P saturation point")

   molecules(1)%groups_ids = [117]
   molecules(1)%number_of_groups = [1]

   molecules(2)%groups_ids = [15]
   molecules(2)%number_of_groups = [1]

   molecules(3)%groups_ids = [62]
   molecules(3)%number_of_groups = [1]

   Tc = [304.21, 512.5 , 720.  ]
   Pc = [73.83, 80.84, 82.  ]
   w = [0.223621, 0.565831, 0.506776]

   model = PSRK(Tc, Pc, w, molecules=molecules)

   z = [0.22_pr, 0.624_pr, 0.156_pr]
   sat = saturation_temperature(model, n=z, P=0.1_pr, kind="dew", t0=300._pr)

   ! call assert(sat%kind == "dew", "Point must converge")
   ! call assert(abs(sat%T - 222._pr) < 1._pr, "Temperature must be close to initial guess")


   kinds_x = ["vapor"]
   kind_w = "liquid"
   x_l0(1, :) = [0.22_pr, 0.624_pr, 0.156_pr]
   w0 = sat%x
   T0 = sat%T
   P0 = sat%P
   betas0 = 1.0_pr
   env = pt_envelope(&
         model, z=z, np=np, x_l0=x_l0, w0=w0, T0=T0, P0=P0, betas0=betas0, &
        kinds_x=kinds_x, kind_w=kind_w, beta_w=0.0_pr, &
        ns0=nc*np+np+2, ds0=1e-2_pr, &
        points=300 &
        )
   open(1, file="psrkm1.dat", status="replace")
   call env%write(1)
   close(1)


   z0 = [0., 0.8, 0.2]
   zi = [1., 0., 0.]

   a = 0.4
   z = a * zi + (1.-a) * z0
   
   x_l0(1, :) = z
   w0 = sat%x
   T0 = sat%T
   P0 = sat%P
   betas0 = 1.0_pr
   env = pt_envelope(&
         model, z=z, np=np, x_l0=x_l0, w0=w0, T0=T0, P0=P0, betas0=betas0, &
        kinds_x=kinds_x, kind_w=kind_w, beta_w=0.0_pr, &
        ns0=nc*np+np+2, ds0=1e-2_pr, &
        points=300 &
        )
   open(1, file="psrkm2.dat", status="replace")
   call env%write(1)
   close(1)
   
   a = 0.689
   z = a * zi + (1.-a) * z0
   
   x_l0(1, :) = z
   w0 = sat%x
   T0 = sat%T
   P0 = sat%P
   betas0 = 1.0_pr
   env = pt_envelope(&
         model, z=z, np=np, x_l0=x_l0, w0=w0, T0=T0, P0=P0, betas0=betas0, &
        kinds_x=kinds_x, kind_w=kind_w, beta_w=0.0_pr, &
        ns0=nc*np+np+2, ds0=1e-2_pr, &
        points=400 &
        )
   open(1, file="psrkm3.dat", status="replace")
   call env%write(1)
   close(1)
end program main
