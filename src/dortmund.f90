module yaeos__models_ge_group_contribution_dortmund
   !! # Dortmund UNIFAC model
   !!
   !! # Description
   !! Modified UNIFAC (Dortmund) group-contribution model. Use
   !! [[setup_dortmund]] to build the model from the groups of each molecule.
   !! The constructor loads the Dortmund parameter tables and installs the
   !! Dortmund temperature dependence on the shared UNIFAC implementation.
   !!
   !! # Examples
   !! ```fortran
   !! model = setup_dortmund(molecules)
   !! call model%ln_activity_coefficient(n, T, ln_gamma)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__models_ge_gc_td, only: QuadraticPsi
   use yaeos__models_ge_group_contribution_groups, only: Groups
   use yaeos__models_ge_group_contribution_unifac, only: UNIFAC, setup_unifac
   use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
   use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

   implicit none

contains

   type(UNIFAC) function setup_dortmund(molecules, parameters)
      !! Create a Dortmund UNIFAC model.
      !!
      !! Uses the Dortmund parameter set by default. An alternative parameter
      !! object can be supplied with `parameters`.
      use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      type(Groups), intent(in) :: molecules(:) !! Groups of each molecule
      type(GeGCModelParameters), optional, intent(in) :: parameters
      !! Optional parameter set
      
      type(GeGCModelParameters) :: params
      type(QuadraticPsi) :: psi_function

      real(pr), allocatable :: Aij(:, :)
      real(pr), allocatable :: Bij(:, :)
      real(pr), allocatable :: Cij(:, :)
      type(Groups) :: soup
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: ng

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
