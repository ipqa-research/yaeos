program main
   use yaeos
   use fixtures_models, only: multicomponent_PR
   use testing_aux, only: assert, test_title
   implicit none

   type(CubicEoS) :: eos
   type(PTEnvel2) :: env
   type(EquilibriumState) :: sat, cp
   integer, parameter :: nc=12
   real(pr) :: z0(nc), zi(nc)
   real(pr) :: z(nc), P, T

   real(pr) :: Tc=699.059
   real(pr) :: Pc=180.226


   write(*, *) test_title("PT envelope test multicomponent")

   eos = multicomponent_PR(z0, zi)

   z = z0
   P = 0.0001
   sat = saturation_temperature(eos, z, P, kind="dew")

   env = pt_envelope_2ph(eos, z, sat)
   cp = critical_point(eos, z0, zi=0*z0, spec=spec_CP%a, S=0._pr, max_iters=100)
   call assert(abs(env%cps(1)%T - Tc)/Tc < 1e-1, "Critical Temperature")
   call assert(abs(env%cps(1)%P - Pc)/pc < 1e-1, "Critical Pressure")

   call test_1

contains

   subroutine test_1
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

      sat = saturation_temperature(model, z, 0.01_pr, kind="dew")

      x_l0(1, :) = sat%y
      w0 = sat%x
      betas0(1) = 1
      p0 = sat%P
      t0 = sat%T
      ns0 = np*nc+np+2
      ds0=1e-1_pr
      beta_w = 0
      kinds_x = "vapor"
      kind_w = "liquid"
      dew = pt_envelope(model, z, np, kinds_x, kind_w, x_l0, w0, betas0, p0, t0, ns0, ds0, beta_w, max_pressure=1500._pr)
      idx = size(dew%points)
      call assert(size(dew%Tc) == 2, "Two critical points found")
      call assert(dew%points(idx)%P > 1000._pr, "Envelope should end at high pressure")
   end subroutine test_1

end program main
