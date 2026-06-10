program test
   use yaeos
   use testing_aux, only: test_title, assert

   write(*, *) "ENVELOPES DATABANK"
   call test_b3_should_not_get_back
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
end program test
