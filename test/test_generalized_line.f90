program generalized_line
   use testing_aux, only: test_title, assert
   use yaeos
   use yaeos__equilibria_boundaries_generalized_isopleths, only: create_generalized_isoz_line, GeneralizedIsoZLine
   implicit none
   integer, parameter :: nc = 3
   real(pr) :: Tc(nc), Pc(nc), w(nc)
   real(pr) :: kij(nc,nc)
   type(CubicEoS) :: model
   type(PTEnvelMP) :: env
   type(EquilibriumState) :: sat
   type(EquilibriumState) :: bub
   real(pr) :: T, P, z(nc)

   integer :: i
   
   write(*, *) test_title("GENERALIZED ISOPLETHS")

   Tc = [190.564, 304.1282, 768.0]
   Pc = [45.992, 73.773, 10.7]
   w = [0.01142, 0.22394, 0.8805]
   kij(1, :) = [0.,  0.1, 0. ]
   kij(2, :) = [0.1, 0.,  0.2]
   kij(3, :) = [0.,  0.2, 0. ]


   model = PengRobinson78(Tc, Pc, w, kij=kij)

   z = [0.2, 0.4, 0.4]

   P = 1
   sat = saturation_temperature(model, z, P, kind="dew")
   bub = saturation_temperature(model, z, P=1._pr, kind="bubble", t0=150._pr)
   call isoP

contains
   subroutine isoP
      integer, parameter :: np=1
      real(pr) :: x_l0(np, nc), w0(nc), betas0(np+1), P0, T0
      type(GeneralizedIsoZLine) :: line
      character(len=14) :: kinds_x(np), kind_w
      integer :: spec_variable
      real(pr) :: spec_variable_value
      integer :: ns0
      real(pr) :: S0, dS0
      integer :: nstab=2
      real(pr) :: ws_stab(2, nc)

      x_l0(1, :) = sat%y
      w0 = sat%x
      betas0 = [1-1e-10, 1e-10]
      P0 = sat%P
      T0 = sat%T
      kinds_x = "vapor"
      kind_w = "liquid"

      spec_variable = (nc*np) + (np + 1) + 1
      spec_variable_value = log(P0)

      ns0 = (nc*np) + (np+1) + 2
      S0 = log(T0) - 0.001
      dS0 = -0.1

      ! ws_stab(1, :) = [4.79003592e-01, 5.20996408e-01, 1.30001898e-15]
      ! ws_stab(2, :) = [5.63039080e-01, 4.36960920e-01, 4.17032026e-16]
      ws_stab(1, :) = bub%y
      ws_stab(2, :) = 1 - sat%x

      print *, sat%y
      print *, sat%x
      print *, sat%P
      print *, sat%T
      print *, "======="

      line = create_generalized_isoz_line(&
         model, nc, np, nstab, kinds_x, kind_w, z, &
         x_l0, w0, betas0, P0, T0, &
         spec_variable, spec_variable_value, ns0, S0, dS0, &
         ws_stab &
         )

      i = size(line%points)
      print *, i, line%points(1)%T, line%points(i)%T
      call assert(line%points(1)%T > 570._pr, "Line sart at high T")
      call assert(line%points(i)%T < 150._pr, "Line stop at low T")

   end subroutine isoP
end program generalized_line
