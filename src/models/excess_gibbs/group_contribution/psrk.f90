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
      
      ! ========================================================================
      ! Build Aij, Bij and Cij matrix (interaction of the soup's subgroups)
      ! ------------------------------------------------------------------------
      matrices: block
      real(pr) :: Aij(size(setup_psrk%groups_stew%groups_ids), size(setup_psrk%groups_stew%groups_ids))
      real(pr) :: Bij(size(setup_psrk%groups_stew%groups_ids), size(setup_psrk%groups_stew%groups_ids))
      real(pr) :: Cij(size(setup_psrk%groups_stew%groups_ids), size(setup_psrk%groups_stew%groups_ids))
      type(Groups) :: soup
      integer :: i, j

      soup = setup_psrk%groups_stew

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

      end block matrices

      
      deallocate(setup_psrk%psi_function)
      setup_psrk%psi_function = psi_function
   end function
end module