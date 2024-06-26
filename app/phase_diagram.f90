program phase_diagram
   use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos
   use yaeos__phase_equilibria_auxiliar, only: k_wilson
   implicit none

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: init
   type(PTEnvel2) :: envelope

   integer, parameter :: nc=2
   integer :: i
   real(pr) :: n(nc), tc(nc), pc(nc), w(nc), kij(nc, nc), lij(nc, nc), T
   type(Substance) :: sus(nc)

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
   sus(1) = Substance("methane")
   sus(2) = Substance("n-hexane")

   n = [0.9_pr, 0.1_pr]       ! Composition
   tc = sus%critical%critical_temperature%value
   pc = sus%critical%critical_pressure%value/1e5
   w = sus%critical%acentric_factor%value

   kij = 0
   kij(1, 2) = 0.0
   kij(2, 1) = kij(1, 2)

   lij =  0

   model = PengRobinson76(tc, pc, w, kij, lij)

   ! Calculate a dew point
   T = 150
   do i=0, 6
      kij(1, 2) = 0._pr + real(i,pr)/10
      kij(2, 1) = kij(1, 2)
      lij = kij/10

      model = PengRobinson78(tc, pc, w, kij, lij)

      init = saturation_temperature(model, n, P=1._pr, kind="dew", t0=150._pr)
      init%x = 1/k_wilson(model, init%T, init%P) * init%y
      envelope = pt_envelope_2ph(model, n, init, delta_0=0.01_pr, points=1000)
      write(1, *) envelope

      init = saturation_pressure(model, n, T=150._pr, kind="bubble", p0=40._pr)
      envelope = pt_envelope_2ph(model, n, init, points=1000)
      write(1, *) envelope

   end do
  
end program phase_diagram
