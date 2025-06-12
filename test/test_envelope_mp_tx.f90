program main
   !! Test for multi-phase envelope calculation
   !! In this test we calculate two phase envelope from a previously known
   !! double saturation point.
   use testing_aux, only: assert, test_title
   use yaeos__constants, only: pr
   use yaeos
   implicit none
   integer, parameter :: nc=15

   type(CubicEoS) :: model

   real(pr) :: z(nc)

   print *, test_title("Multi-phase TX envelope test")

   model = get_model()

   call calc_2ph
contains

   type(CubicEoS) function get_model()
      real(pr) :: tc(15), pc(15), w(15), kij(15, 15)
      z = [0.0048, 0.00919, 0.43391, 0.1101, 0.06544, 0.00789, 0.03787, &
         0.01279, 0.02248, 0.02698, 0.22738, 0.03747, 0.0023, 0.00054, 0.00086]
      Tc = [126.2, 304.2, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, &
         469.6, 507.4, 691.81, 956.77, 1118.6, 1325.03, 1445.73]
      Pc = [33.94, 73.76, 46., 48.84, 42.46, 36.48, 38., 33.84, 33.74, 29.69, 19.46, 13.08, 10.66, 10.28, 17.3]
      w = [0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.68, 1.208, 0.949, 0.182, 1.274]

      kij(1, :) = [0., -0.032, 0.028, 0.041, 0.076, 0.094, 0.07, 0.087, 0.088, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08]
      kij(2, :) = [-0.032, 0., 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.1, 0.1, 0.1, 0.1, 0.1]
      kij(3, :) = [0.028, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(4, :) = [0.041, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(5, :) = [0.076, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(6, :) = [0.094, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(7, :) = [0.07, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(8, :) = [0.087, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(9, :) = [0.088, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(10, :) = [0.08, 0.12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(11, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(12, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(13, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(14, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(15, :) = [0.08, 0.1, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0., 0., 0., 0., 0.]

      get_model = PengRobinson76(Tc, Pc, w, kij=kij)
   end function get_model

   subroutine calc_2ph
      integer, parameter :: np=1
      real(pr) :: x_l0(np, nc), w0(nc), z0(nc), zi(nc)
      type(EquilibriumState) :: sat
      type(txenvelmp) :: tx
      type(PTEnvelMP) :: PT

      character(len=14) :: kind_w
      character(len=14) :: kinds_x(np)

      real(pr) :: T, P, beta
      integer :: i

      sat = saturation_temperature(model, z, p=1e-3_pr, kind="dew", t0=600._pr)
      x_l0(1, :) = z
      P = 100

      sat = saturation_temperature(model, z, P=P, kind="dew", t0=800._pr)

      w0 = sat%t

      z0 = z
      zi = 0
      zi(2) = 1

      x_l0(1, :) = z
      
      sat = saturation_temperature(model, z, P=P, kind="bubble", t0=300._pr)
      kind_w = "vapor"
      kinds_x(1) = "liquid"
      
      w0 = sat%y
      tx = tx_envelope(&
         model, z0, zi, np, sat%P, kinds_x=kinds_x, kind_w=kind_w, x_l0=x_l0, w0=w0, betas0=[1._pr], &
         T0=sat%t, alpha0=0._pr, ns0=np*nc+np+2, ds0=0.005_pr, &
         beta_w=0.0_pr, points=1000&
         )

      sat = saturation_temperature(model, z, P=P, kind="dew", t0=800._pr)
      w0 = sat%x
      kind_w = "liquid"
      kinds_x(1) = "vapor"
      tx = tx_envelope(&
         model, z0, zi, np, sat%P, kinds_x=kinds_x, kind_w=kind_w, x_l0=x_l0, w0=w0, betas0=[1._pr], &
         T0=sat%t, alpha0=0._pr, ns0=np*nc+np+2, ds0=0.005_pr, &
         beta_w=0.0_pr, points=1000&
         )

      call assert(tx%alpha(size(tx%alpha)) > 0.99, "final alpha")

   end subroutine calc_2ph
end program main
