module yaeos__models_ge_group_contribution_dortmund
   use yaeos__constants, only: pr
   use yaeos__models_ge_gc_td, only: QuadraticPsi
   use yaeos__models_ge_group_contribution_groups, only: Groups
   use yaeos__models_ge_group_contribution_unifac, only: UNIFAC, setup_unifac
   use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
   use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

   implicit none

contains

   type(UNIFAC) function setup_dortmund(molecules, parameters)
      use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      type(Groups), intent(in) :: molecules(:)
      type(GeGCModelParameters), optional, intent(in) :: parameters
      
      type(GeGCModelParameters) :: params
      type(QuadraticPsi) :: psi_function

      real(pr), allocatable :: Aij(:, :), Bij(:, :), Cij(:, :)
      type(Groups) :: soup
      integer :: i, j, ng

      if (present(parameters)) then
         params = parameters
      else
         params = DortmundParameters()
      end if

      setup_dortmund = setup_unifac(molecules, params)
      
      ! ========================================================================
      ! Build Aij, Bij and Cij matrix (interaction of the soup's subgroups)
      ! ------------------------------------------------------------------------
      soup = setup_dortmund%groups_stew
      ng = size(soup%groups_ids)

      allocate(Aij(ng, ng), Bij(ng, ng), Cij(ng, ng))

      Aij = 0
      Bij = 0
      Cij = 0

      do i=1,size(soup%groups_ids)
         do j=1,size(soup%groups_ids)
            Aij(i, j) = params%get_subgroups_aij(&
               soup%groups_ids(i), soup%groups_ids(j) &
               )
            Bij(i, j) = params%get_subgroups_bij(&
               soup%groups_ids(i), soup%groups_ids(j) &
               )
            Cij(i, j) = params%get_subgroups_cij(&
               soup%groups_ids(i), soup%groups_ids(j) &
               )
         end do
      end do
      
      psi_function%Aij = Aij
      psi_function%Bij = Bij
      psi_function%Cij = Cij
      
      deallocate(setup_dortmund%psi_function)
      setup_dortmund%psi_function = psi_function

      ! Important, the parameter d (exponent of the r params in Flory-Huggins)
      setup_dortmund%d = 3.0_pr / 4.0_pr
   end function
end module