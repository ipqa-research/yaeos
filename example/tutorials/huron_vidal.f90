program main
   !! Example of using CubicEoS with Huron-Vidal mixing rules with an
   !! NRTL model as the Ge model
   use forsus, only: Substance, forsus_dir
   use yaeos, only: &
      pr, SoaveRedlichKwong, CubicEoS, NRTL, saturation_pressure, &
         pt_envelope_2ph, EquilibriaState, PTEnvel2, UNIFAC, setup_unifac, Groups
   use yaeos__models_cubic_mixing_rules_huron_vidal, only: MHV

   implicit none
   integer, parameter :: nc = 2

   real(pr) :: n(nc)

   real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc) ! NRTL parameters
   real(pr) :: tc(nc), pc(nc), w(nc) ! Cubic EoS parameters

   type(NRTL) :: ge_model ! Excess Gibbs model that will be used
   type(UNIFAC) :: ge_model_unifac
   type(CubicEoS) :: model ! Main model
   type(MHV) :: mixrule
   type(Groups) :: molecules(nc)

   type(EquilibriaState) :: sat

   type(Substance)  :: sus(nc)

   integer :: i, j

   molecules(1)%groups_ids = [16]
   molecules(1)%number_of_groups = [1]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   forsus_dir = "./build/dependencies/forsus/data/json"


   sus(1) = Substance("water")
   sus(2) = Substance("ethanol")

   tc = sus%critical%critical_temperature%value
   w = sus%critical%acentric_factor%value
   pc = sus%critical%critical_pressure%value/1e5

   a = 0; b = 0; c = 0

   ! NRTL model parameters
   a(1, 2) = 3.458
   a(2, 1) = -0.801

   b(1, 2) = -586.1
   b(2, 1) = 246.2

   c(1, 2) = 0.3
   c(2, 1) = 0.3

   ge_model = NRTL(a, b, c)

   n = [0.9, 0.1]
   ! n = [0.8, 0.2]
   ! Define the model to be SRK
   model = SoaveRedlichKwong(tc, pc, w)
   call phase_envel(1)
   
   mixrule = MHV(ge=ge_model, q=-0.593_pr, b=model%b)
   deallocate (model%mixrule)
   model%mixrule = mixrule
   call phase_envel(2)

   ge_model_unifac = setup_unifac(molecules)
  
   mixrule = MHV(ge=ge_model_unifac, q=-0.593_pr, b=model%b)
   deallocate (model%mixrule)
   model%mixrule = mixrule
   call phase_envel(3)

   do i=1,99
      n(2) = real(i,pr)/100
      n(1) = 1 - n(2)
      sat = saturation_pressure(model, n, T=473._pr, kind="bubble")
      ! write (*, *) sat
   end do

contains
   
   subroutine phase_envel(fu)
      use yaeos, only: EquilibriaState, PTEnvel2, pt_envelope_2ph, saturation_temperature
      integer :: fu
      type(EquilibriaState) :: sat
      type(PTEnvel2) :: env

      sat = saturation_pressure(model, n, T=300._pr, kind="bubble")
      write (*, *) sat, sat%iters

      env = pt_envelope_2ph(model, n, sat, specified_variable_0=nc + 1, delta_0=0.001_pr)
      write (fu, *) env
   end subroutine
end program
