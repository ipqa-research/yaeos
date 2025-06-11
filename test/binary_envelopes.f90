program main
   use yaeos
   use testing_aux, only: assert, test_title
   use fixtures_models, only: binary_PR76

   implicit none

   integer, parameter :: nc=2
   type(CubicEoS) :: eos

   write(*, *) test_title("Phase envelopes of binary systems")

   call pt
   call pxy
   call txy

contains

   subroutine pt
      type(PTEnvelMP) :: env
      type(EquilibriumState) :: sat
      real(pr) :: z(nc), P, T

      real(pr) :: betas0(1), x_l0(1, nc), w0(nc)
      character(len=14) :: kinds_x(1), kind_w

      eos = binary_PR76()


      z = [0.5_pr, 0.5_pr]
      P = 0.0001_pr
      sat = saturation_temperature(eos, z, P, kind="dew")

      w0 = sat%x
      x_l0(1, :) = z
      betas0(1) = 1._pr
      kinds_x = "liquid"
      kind_w = "vapor"
      env = pt_envelope(model=eos, z=z, np=1, kinds_x=kinds_x, kind_w=kind_w, &
         x_l0=x_l0, w0=w0, betas0=betas0, P0=sat%P, T0=sat%T, ns0=nc+1+1, dS0=1e-3_pr, beta_w=0.0_pr, points=300, &
         max_pressure=2500._pr)
   end subroutine pt

   subroutine pxy
      type(PXEnvelMP) :: env
      type(EquilibriumState) :: sat
      real(pr) :: zi(nc), z0(nc), z(nc), P, T

      real(pr) :: alpha0, betas0(1), x_l0(1, nc), w0(nc)
      character(len=14) :: kinds_x(1), kind_w

      integer :: i

      eos = binary_PR76()

      zi = [1, 0]
      z0 = [0, 1]
      alpha0 = 1e-5
      z = alpha0 * zi + (1.0_pr - alpha0) * z0

      T = 200
      sat = saturation_pressure(eos, z, T=T, kind="bubble")

      w0 = sat%y
      x_l0(1, :) = z
      betas0(1) = 1._pr
      kinds_x = "liquid"
      kind_w = "vapor"

      env = px_envelope(&
         eos, z0=z0, zi=zi, np=1, T=T, &
         x_l0=x_l0, w0=w0, betas0=betas0, P0=sat%P, alpha0=alpha0, &
         ns0=nc+1+2, ds0=1e-5_pr, beta_w=0.0_pr,  points=800, kinds_x=kinds_x, kind_w=kind_w&
         )

      i = size(env%alpha)
      call assert(all(abs(env%points(i)%x_l(1, :) - env%points(i)%w) < 1e-3), "Pxy: End a critical point")
   end subroutine pxy

   subroutine txy
      type(TXEnvelMP) :: env
      type(EquilibriumState) :: sat
      type(PurePsat) :: psat
      real(pr) :: zi(nc), z0(nc), z(nc), P, T

      real(pr) :: alpha0, betas0(1), x_l0(1, nc), w0(nc)
      character(len=14) :: kinds_x(1), kind_w

      integer :: i

      eos = binary_PR76()

      zi = [1, 0]
      z0 = [0, 1]
      alpha0 = 1e-3
      z = alpha0 * zi + (1.0_pr - alpha0) * z0

      P = 20

      psat = pure_saturation_line(eos, 2, 1e-5_pr, 100._pr)
      T = psat%get_T(P)

      sat = saturation_temperature(eos, z, P=P, kind="bubble", T0=T)

      w0 = sat%y
      x_l0(1, :) = z
      betas0(1) = 1._pr
      kinds_x = "liquid"
      kind_w = "vapor"

      env = tx_envelope(&
         eos, z0=z0, zi=zi, np=1, P=P, &
         x_l0=x_l0, w0=w0, betas0=betas0, T0=sat%T, alpha0=alpha0, &
         ns0=nc+1+2, ds0=1e-5_pr, beta_w=0.0_pr,  points=800, &
         kinds_x=kinds_x, kind_w=kind_w&
         )
      i = size(env%alpha)
      call assert(all(abs(env%points(i)%x_l(1, :) - env%points(i)%w) < 1e-2), "Txy: End a critical point")
   end subroutine txy
end program main
