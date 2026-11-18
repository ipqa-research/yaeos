module yaeos__models_ge_group_contribution_psrk
   !! # PSRK UNIFAC model
   !!
   !! # Description
   !! Group-contribution model used by the PSRK equation of state. Use
   !! [[setup_psrk]] to build the model from the groups of each molecule.
   !! Uses the PSRK group parameter set with the UNIFAC combinatorial and
   !! residual contributions.
   !!
   !! # Examples
   !! ```fortran
   !! model = setup_psrk(molecules)
   !! call model%ln_activity_coefficient(n, T, ln_gamma)
   !! ```
   use yaeos__constants, only: pr
   use yaeos__models_ge_gc_td, only: QuadraticPsi
   use yaeos__models_ge_group_contribution_groups, only: Groups
   use yaeos__models_ge_group_contribution_unifac, only: UNIFAC, setup_unifac
   use yaeos__models_ge_group_contribution_psrk_parameters, only: PSRKParameters, GeGCModelParameters

   implicit none

contains

   type(UNIFAC) function setup_psrk(molecules, parameters)
      !! Create a PSRK UNIFAC model.
      !!
      !! Uses the PSRK parameter set by default. An alternative parameter
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
         params = PSRKParameters()
      end if

      setup_psrk = setup_unifac(molecules, params)
      
      ! ========================================================================
      ! Build Aij, Bij and Cij matrix (interaction of the soup's subgroups)
      ! ------------------------------------------------------------------------
      soup = setup_psrk%groups_stew
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
      
      deallocate(setup_psrk%psi_function)
      setup_psrk%psi_function = psi_function
   end function
end module
