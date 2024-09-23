module yaeos__models_ge_group_contribution_psrk
   use yaeos__constants, only: pr
   use yaeos__models_ge_gc_td, only: QuadraticPsi
   use yaeos__models_ge_group_contribution_groups, only: Groups
   use yaeos__models_ge_group_contribution_unifac, only: UNIFAC, setup_unifac
   use yaeos__models_ge_group_contribution_psrk_parameters, only: PSRKParameters, GeGCModelParameters

   implicit none

contains

   type(UNIFAC) function setup_psrk(molecules, parameters)
      use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      type(Groups), intent(in) :: molecules(:)
      type(GeGCModelParameters), optional, intent(in) :: parameters
      
      type(GeGCModelParameters) :: params
      type(QuadraticPsi) :: psi_function

      if (present(parameters)) then
         params = parameters
      else
         params = PSRKParameters()
      end if

      setup_psrk = setup_unifac(molecules, params)
      
      psi_function%Aij = params%maingroups_aij
      psi_function%Bij = params%maingroups_bij
      psi_function%Cij = params%maingroups_cij
      deallocate(setup_psrk%psi_function)
      setup_psrk%psi_function = psi_function
   end function
end module