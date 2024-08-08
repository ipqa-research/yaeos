program phase_diagram
   !! Program for calculation of phase diagrams. 
   use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos, only: pr, &
      SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, &
      EquilibriumState, ArModel, PTEnvel2, &
      pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=2            
   class(ArModel), allocatable :: model ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point   ! Init
   type(PTEnvel2) :: envelope           ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)    ! Component's critical constants
   real(pr) :: n(nc)                    ! Termodynamic variables
   type(Substance) :: sus(nc)           ! Substances to use
   ! ===========================================================================

   ! forsus database directory
   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   ! Find the selected substances on the database and extract their
   ! critical constants
   sus(1) = Substance("methane")
   sus(2) = Substance("n-hexane")
   call get_critical_constants(sus, tc, pc, w)

   ! Model definition
   model = PengRobinson76(tc, pc, w)

   ! Composition vector
   n = [0.9_pr, 0.1_pr]
   
   ! Calculate a dew point at low pressure to later 
   ! initialize the phase envelope
   sat_point = saturation_temperature(model, n, P=1._pr, kind="dew", t0=150._pr)

   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, n, sat_point)

   ! Write the phase envelope to screen
   write(*, *) envelope

contains

   subroutine get_critical_constants(subs, tc_in, pc_in, w_in)
      type(Substance) :: subs(:)
      real(pr), intent(out) :: tc_in(:), pc_in(:), w_in(:)

      tc_in = subs%critical%critical_temperature%value
      pc_in = subs%critical%critical_pressure%value/1e5
      w_in = subs%critical%acentric_factor%value
   end subroutine
end program phase_diagram
