program phase_diagram
   use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos
   use yaeos__phase_equilibria_auxiliar, only: k_wilson
   implicit none

   integer, parameter :: nc=2
   class(ArModel), allocatable :: model
   type(EquilibriaState) :: init
   type(PTEnvel2) :: envelope

   real(pr) :: n(nc), tc(nc), pc(nc), w(nc), T
   type(Substance) :: sus(nc)

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   sus(1) = Substance("methane")
   sus(2) = Substance("n-hexane")
   call get_critical_constants(sus, tc, pc, w)

   model = PengRobinson76(tc, pc, w)

   n = [0.9_pr, 0.1_pr]       ! Composition
   T = 150
   init = saturation_temperature(model, n, P=1._pr, kind="dew", t0=150._pr)
   envelope = pt_envelope_2ph(model, n, init, delta_0=0.01_pr, points=1000)
   write(*, *) envelope

contains

   subroutine get_critical_constants(subs, tc, pc, w)
      type(Substance) :: subs(:)
      real(pr) :: tc(:), pc(:), w(:)
      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5
      w = sus%critical%acentric_factor%value
   end subroutine
end program phase_diagram
