program test
   use yaeos
   use testing_aux, only: test_title, assert

   write(*, *) "ENVELOPES DATABANK"
   call test_b3_should_not_get_back
   call test_double_cp
   call test_b3_should_not_get_back_2
contains
   subroutine test_b3_should_not_get_back
      use yaeos__extra_fluids, only: oil_gao, FixtureFluid
      integer, parameter :: nc = 6, np=2
      type(FixtureFluid) :: fluid
      class(ArModel), allocatable :: model
      real(pr) :: z(nc), x_l(np, nc), w(nc), T, P, betas(np)
      character(len=14) :: k_x(np)
      character(len=14) :: k_w
      real(pr) :: mintpd

      type(PTEnvelMP) :: env
      type(EquilibriumState) :: bub, fr, dew
      integer :: its, points

      fluid = oil_gao()
      z = fluid%z0 ! / sum(fluid%z0)
      bub = saturation_pressure(fluid%ar_model, z, T=200._pr, kind="bubble", p0=3._pr)
      dew = saturation_temperature(fluid%ar_model, z, P=1._pr, kind="dew", t0=300._pr)

      call min_tpd(fluid%ar_model, z, bub%P, bub%T, mintpd, w)

      fr = flash(fluid%ar_model, z=z, P_spec=bub%P, T=bub%T, k0=w/z, iters=its)

      x_l(1, :) = fr%x
      x_l(2, :) = fr%y
      w = bub%y
      betas = [1-fr%beta, fr%beta]
      P = bub%P
      T = bub%T

      k_x = "liquid"
      k_w = "vapor"
      env = pt_envelope(&
         model=fluid%ar_model, np=np, z=z, x_l0=x_l, w0=w, betas0=betas, p0=P, t0=T, &
         kinds_x=k_x, kind_w=k_w, &
         ns0=nc*np+np+2, ds0=1e-3_pr, beta_w=0._pr &
         )
      points = size(env%points)

      call assert(env%points(points)%T > 400, "PT3 envelope should break instead of getting back over itself")

   end subroutine test_b3_should_not_get_back

   subroutine test_double_cp
      integer, parameter :: nc=3, np=1
      real(pr) :: tc(nc), pc(nc), w(nc), z(nc)
      real(pr) :: kij(nc, nc)
      type(CubicEoS) :: model
      type(PTEnvelMP) :: dew
      type(EquilibriumState) :: sat
      real(pr) :: x_l0(np, nc), w0(nc), betas0(np), p0, t0
      character(len=14) :: kinds_x(np), kind_w
      integer :: ns0
      real(pr) :: ds0, beta_w
      integer :: idx


      tc = [304.21, 373.53, 190.564]
      pc = [73.83000000000001, 89.62910000000001, 45.99]
      w = [0.223621, 0.0941677, 0.0115478]

      kij = 0
      kij(1, 2) = 0.0974
      kij(2, 1) = 0.0974
      kij(1, 3) = 0.110
      kij(3, 1) = 0.110
      kij(2, 3) = 0.069
      kij(3, 2) = 0.069


      model = PengRobinson76(tc, pc, w, kij=kij)
      z = [0.0987, 0.4023, 0.4990]
      z = z / sum(z)

      p0 = 0.1
      sat = saturation_temperature(model, z, p0, kind="dew", T0=150._pr)

      x_l0(1, :) = sat%y
      w0 = sat%x
      betas0(1) = 1
      p0 = sat%P
      t0 = sat%T
      ns0 = np*nc+np+1
      ds0=1e-1
      beta_w = 0
      kinds_x = "vapor"
      kind_w = "liquid"
      dew = pt_envelope(model, z, np, kinds_x, kind_w, x_l0, w0, betas0, p0, t0, ns0, ds0, beta_w, max_pressure=1500._pr)
      idx = size(dew%points)
      call assert(size(dew%Tc) == 2, "Double CP: Two critical points not found")
      call assert(dew%points(idx)%P > 1000._pr, "Double CP: Envelope should end at high pressure")
   end subroutine test_double_cp

   subroutine test_multicomponent
      use fixtures_models, only: multicomponent_PR

      type(CubicEoS) :: eos
      type(PTEnvel2) :: env
      type(EquilibriumState) :: sat, cp
      integer, parameter :: nc=12
      real(pr) :: z0(nc), zi(nc)
      real(pr) :: z(nc), P, T

      real(pr) :: Tc=699.059
      real(pr) :: Pc=180.226

      eos = multicomponent_PR(z0, zi)

      z = z0
      P = 0.0001
      sat = saturation_temperature(eos, z, P, kind="dew")

      env = pt_envelope_2ph(eos, z, sat, points=1000)
      call assert(abs(env%cps(1)%T - Tc)/Tc < 1e-1, "multicomponent_PR: Critical Temperature")
      call assert(abs(env%cps(1)%P - Pc)/pc < 1e-1, "multicomponent_PR: Critical Pressure")
   end subroutine test_multicomponent

   subroutine test_b3_should_not_get_back_2
      use yaeos__extra_fluids, only: oil_b71, FixtureFluid
      use yaeos, only: pt_envelope, PTEnvelMP, ArModel
      type(FixtureFluid) :: fluid
      integer, parameter :: nc = 16, np=2
      real(pr) :: z(nc), x_l(np, nc), w(nc), T, P, betas(np)
      character(len=14) :: k_x(np)
      character(len=14) :: k_w
      
      type(PTEnvelMP) :: env

      fluid = oil_b71()

      z = fluid%z0
      x_l(1, :) = [0.9416741701541934_pr, 0.0011213305701257231_pr, 0.03701631257908194_pr, 0.006399527236767095_pr, 0.003817604608926774_pr, 0.0002966023186278051_pr, 0.0030372769086758223_pr, 0.0010069560277233447_pr, 0.0013898065767137563_pr, 0.0018678423174985478_pr, 0.0022016092602370452_pr, 0.0001685391251062202_pr, 2.418477580604817e-06_pr, 3.8386836693575555e-09_pr, 5.812460558552299e-14_pr, 8.213791709817468e-25_pr]
      x_l(2, :) = [0.35045087301301103_pr, 0.0018130206354868726_pr, 0.06321939626556232_pr, 0.019911553263846114_pr, 0.016072137206215923_pr, 0.0022065614570470613_pr, 0.019657405203985616_pr, 0.010142906105283357_pr, 0.013771597003497274_pr, 0.021698762924659325_pr, 0.1330484520274068_pr, 0.12119564721915836_pr, 0.09430425135410335_pr, 0.07180261856002143_pr, 0.04358195796601319_pr, 0.017122859794701805_pr]
      w = [0.4735720485900767_pr, 0.052429758529525014_pr, 0.46949954475431144_pr, 0.004083537850152927_pr, 0.0003471229043833454_pr, 9.691898117208192e-06_pr, 4.853080578189653e-05_pr, 4.72601393083613e-06_pr, 3.8040860361621804e-06_pr, 1.121841630098261e-06_pr, 1.1258574146363489e-07_pr, 1.4620244416616193e-10_pr, 4.196413708634265e-14_pr, 2.64232988450528e-18_pr, 2.4455204481216936e-23_pr, 1.3740466369262361e-30_pr]
      betas = [1 - 0.39131605820152016_pr, 0.39131605820152016_pr]
      P = 7.344294594110742
      T = 199.99999999999991
      k_x = "liquid"
      k_w = "vapor"

      env = pt_envelope(&
         model=fluid%ar_model, np=np, z=z, x_l0=x_l, w0=w, betas0=betas, p0=P, t0=T, &
         kinds_x=k_x, kind_w=k_w, &
         ns0=nc*np+np+2, ds0=1e-3_pr, beta_w=0._pr &
         )
      open(1, file="TEST_ENVELOPE.dat")
      call env%write(1)
      close(1)
   end subroutine test_b3_should_not_get_back_2
end program test
