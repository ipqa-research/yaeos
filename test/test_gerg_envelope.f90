program main 
   use yaeos
   use testing_aux, only: test_title, assert
   implicit none
   integer, parameter :: nc = 5
   integer, parameter :: np = 1

   type(PTEnvelMP) :: env
   type(GERG2008) :: model
   integer :: ids(nc)
   type(EquilibriumState) :: sat
   real(pr) :: z(nc), x_l(1, nc), w(nc)
   character(len=14) :: kinds_x(np), kind_w
   integer :: i

   real(pr) :: last_temperature

   write(*, *) test_title("PT envelope with GERG-2008")

   ids = [ &
      g2008components%carbon_dioxide, &
      g2008components%methane, &
      g2008components%ethane, &
      g2008components%propane, &
      g2008components%nbutane & 
   ]
   model = gerg_2008(ids)

   z = [0.1356_pr, 0.6636_pr, 0.09975_pr, 0.06677_pr, 0.0343_pr]
   z = z/sum(z)

   sat = saturation_temperature(model, z, P=1._pr, kind="dew", t0=150._pr)

   x_l(1, :) = z
   w = sat%x
   kinds_x = "vapor"
   kind_w = "liquid"

   env = pt_envelope(model, z, 1, kinds_x, kind_w, x_l, w, [1._pr], sat%P, sat%T, ns0=nc+3, ds0=1e-3_pr, beta_w=0._pr)

   last_temperature = env%points(size(env%points))%T

   call assert(last_temperature < 210._pr, "Envelope temperature should be less than 210 K")

end program