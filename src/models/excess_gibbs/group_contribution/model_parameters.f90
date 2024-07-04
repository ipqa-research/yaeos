module yaeos__models_ge_group_contribution_model_parameters
   !! # \(G^E\) group contribution model parameters
   !! \(G^E\) group contribution model parameters module.
   !!
   !! # Description
   !! This module contrains the GeGCModelParameters type that allows to store
   !! the subgroups ids, maingroups ids, subgroups Rs, subgroups Qs,
   !! subgroups maingroups, and maingroups interaction parameters for UNIFAC
   !! like models (UNIFAC, LL-UNIFAC, Dortmund UNIFAC, PSRK, etc)
   !!
   use yaeos__constants, only: pr
   implicit none

   type :: GeGCModelParameters
      !! # GeGCModelParameters
      !! \(G^E\) group contribution model parameters container
      !!
      !! # Description
      !! Type to represent a UNIFAC like models parameters. The type must be
      !! provided with the subgroups ids, maingroups ids, subgroups Rs,
      !! subgroups Qs, subgroups maingroups, and maingroups interaction
      !! parameters. Specifically, the type requires \(a_{ij}\), \(b_{ij}\), and
      !! \(c_{ij}\) for the maingroups interaction parameters. In the case of
      !! the classic UNIFAC model that only requires \(a_{ij}\) parameters, the
      !! \(b_{ij}\) and \(c_{ij}\) must be set as null matrixes.
      !! The documentation and source code of `yaeos` [[UNIFACParameters]]
      !! function could be consulted to understand how to instantiate a
      !! [[GeGCModelParameters]] object with the classic liquid-vapor UNIFAC
      !! parameters defined in DDBST.
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      integer, allocatable :: subgroups_ids(:)
      !! ID of each model's subgroup
      integer, allocatable :: maingroups_ids(:)
      !! ID of each model's maingroup
      integer, allocatable :: subgroups_maingroups(:)
      !! Maingroup of each subgroup
      real(pr), allocatable :: subgroups_Rs(:)
      !! \(R \) value of each subgroup
      real(pr), allocatable :: subgroups_Qs(:)
      !! \(Q \) value of each subgroup
      real(pr), allocatable :: maingroups_aij(:,:)
      !! Maingroup \(a_{ij} \) interaction parameters matrix
      real(pr), allocatable :: maingroups_bij(:,:)
      !! Maingroup \(b_{ij} \) interaction parameters matrix
      real(pr), allocatable :: maingroups_cij(:,:)
      !! Maingroup \(c_{ij} \) interaction parameters matrix
   contains
      procedure :: get_subgroup_index => get_subgroup_index
      procedure :: get_maingroup_index => get_maingroup_index
      procedure :: get_subgroup_maingroup => get_subgroup_maingroup
      procedure :: get_subgroup_R => get_subgroup_R
      procedure :: get_subgroup_Q => get_subgroup_Q
      procedure :: get_maingroups_aij => get_maingroups_aij
      procedure :: get_maingroups_bij => get_maingroups_bij
      procedure :: get_maingroups_cij => get_maingroups_cij
      procedure :: get_subgroups_aij => get_subgroups_aij
      procedure :: get_subgroups_bij => get_subgroups_bij
      procedure :: get_subgroups_cij => get_subgroups_cij
   end type GeGCModelParameters

contains
   function get_subgroup_index(self, subgroup_id) result(subgroup_idx)
      !! # get_subgroup_index
      !! Get index of the subgroup with id: `subgroup_id`
      !!
      !! # Description
      !! Get index of the subgroup with id: `subgroup_id`. Gets the index of the
      !! subgroup in the `self%subgroups_ids` vector.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  ! Default parameters of UNIFAC (ddbst)
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get index of the subgroup with id 178 (IMIDAZOL)
      !!  print *, parameters%get_subgroup_index(178) ! Will print: 112
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      !! ID of the subgroup
      integer :: subgroup_idx
      !! Index of the subgroup on the `self%subgroups_ids` vector

      subgroup_idx = findloc(self%subgroups_ids, subgroup_id, dim=1)
   end function get_subgroup_index

   function get_maingroup_index(self, maingroup_id) result(maingroup_idx)
      !! # get_maingroup_index
      !! Get index of the maingoup with id: `maingoup_id`
      !!
      !! # Description
      !! Get index of the maingoup with id: `maingoup_id`. Gets the index of the
      !! maingoup in the `self%maingoups_ids` vector.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get index of the maingroup with id 55 (Sulfones: [118](CH2)2SU [119]CH2CHSU)
      !!  print *, parameters%get_maingroup_index(55) ! Will print: 52
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_id
      !! ID of the subgroup
      integer :: maingroup_idx
      !! Index of the maingroup on the `self%maingroups_ids` vector

      maingroup_idx = findloc(self%maingroups_ids, maingroup_id, dim=1)
   end function get_maingroup_index

   function get_subgroup_maingroup(self, subgroup_id) result(subgroup_maingroup)
      !! # get_subgroup_maingroup
      !! Get the subgroup's maingroup
      !!
      !! # Description
      !! Uses the `self%subgroups_maingroups` attribute to locate the maingroup
      !! where the subgroup with id `subgroup_id` belongs
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the maingroup of the subgroup with id 16 (H2O)
      !!  print *, parameters%get_subgroup_maingroup(16) ! Will print: 7
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      !! ID of the subgroup
      integer :: subgroup_maingroup
      !! Maingroup of the subgroup

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_maingroup = self%subgroups_maingroups(subgroup_idx)
   end function get_subgroup_maingroup

   function get_subgroup_R(self, subgroup_id) result(subgroup_R)
      !! # get_subgroup_R
      !! Get the subgroup's \(R \) value
      !!
      !! # Description
      !! Uses the `self%subgroups_Rs` attribute to locate the subgroup \(R \)
      !! value.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the subgroup, with id 1 (CH3), R value
      !!  print *, parameters%get_subgroup_R(1) ! Will print: 0.9011
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      !! ID of the subgroup
      real(pr) :: subgroup_R
      !! \(R \) value of the subgroup

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_R = self%subgroups_Rs(subgroup_idx)
   end function get_subgroup_R

   function get_subgroup_Q(self, subgroup_id) result(subgroup_Q)
      !! # get_subgroup_Q
      !! Get the subgroup's \(Q \) value
      !!
      !! # Description
      !! Uses the `self%subgroups_Qs` attribute to locate the subgroup \(Q \)
      !! value.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the subgroup, with id 1 (CH3), Q value
      !!  print *, parameters%get_subgroup_Q(1) ! Will print: 0.8480
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      !! ID of the subgroup
      real(pr) :: subgroup_Q
      !! \(Q \) value of the subgroup

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_Q = self%subgroups_Qs(subgroup_idx)
   end function get_subgroup_Q

   function get_maingroups_aij(self, maingroup_i_id, maingroup_j_id) result(aij)
      !! # get_maingroups_aij
      !! Get the interaction parameter \(a_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(a_{ij}\) of the maingroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the maingroups i:1, j:7 interaction parameter aij (CH2-H2O)
      !!  print *, parameters%get_maingroups_aij(1, 7) ! prints: 1318.0000
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      !! ID of the maingroup `i`
      integer, intent(in) :: maingroup_j_id
      !! ID of the maingroup `j`
      real(pr) :: aij
      !! Interaction parameter \(a_{ij}\)

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      aij = self%maingroups_aij(i, j)
   end function get_maingroups_aij

   function get_maingroups_bij(self, maingroup_i_id, maingroup_j_id) result(bij)
      !! # get_maingroups_bij
      !! Get the interaction parameter \(b_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(b_{ij}\) of the maingroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the maingroups i:1, j:7 interaction parameter bij (CH2-H2O)
      !!  print *, parameters%get_maingroups_bij(1, 7) ! prints: 0.0
      !! ```
      !!
      !! In the example we obtain 0.0 because UNIFAC only have \(a_{ij}\)
      !! parameters
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      !! ID of the maingroup `i`
      integer, intent(in) :: maingroup_j_id
      !! ID of the maingroup `j`
      real(pr) :: bij
      !! Interaction parameter \(b_{ij}\)

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      bij = self%maingroups_bij(i, j)
   end function get_maingroups_bij

   function get_maingroups_cij(self, maingroup_i_id, maingroup_j_id) result(cij)
      !! # get_maingroups_cij
      !! Get the interaction parameter \(c_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(c_{ij}\) of the maingroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the maingroups i:1, j:7 interaction parameter cij (CH2-H2O)
      !!  print *, parameters%get_maingroups_cij(1, 7) ! prints: 0.0
      !! ```
      !!
      !! In the example we obtain 0.0 because UNIFAC only have \(a_{ij} \)
      !! parameters
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      !! ID of the maingroup `i`
      integer, intent(in) :: maingroup_j_id
      !! ID of the maingroup `j`
      real(pr) :: cij
      !! Interaction parameter \(c_{ij}\)

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      cij = self%maingroups_cij(i, j)
   end function get_maingroups_cij

   function get_subgroups_aij(self, subgroup_i_id, subgroup_j_id) result(aij)
      !! # get_subgroups_aij
      !! Get the interaction parameter \(a_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(a_{ij}\) of the subgroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the subgroups i:1, j:16 interaction parameter aij (CH3-H2O)
      !!  ! with maingroups 1 and 7 respectively.
      !!  print *, parameters%get_subgroups_aij(1, 16) ! prints: 1318.0000
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      !! ID of the subgroup `i`
      integer, intent(in) :: subgroup_j_id
      !! ID of the subgroup `j`
      real(pr) :: aij
      !! Interaction parameter \(a_{ij}\)

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      aij = self%maingroups_aij(i, j)
   end function get_subgroups_aij

   function get_subgroups_bij(self, subgroup_i_id, subgroup_j_id) result(bij)
      !! # get_subgroups_bij
      !! Get the interaction parameter \(b_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(b_{ij}\) of the subgroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the subgroups i:1, j:16 interaction parameter bij (CH3-H2O)
      !!  ! with maingroups 1 and 7 respectively.
      !!  print *, parameters%get_subgroups_bij(1, 16) ! prints: 0.0000
      !! ```
      !!
      !! In the example we obtain 0.0 because UNIFAC only have \(a_{ij} \)
      !! parameters
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      !! ID of the subgroup `i`
      integer, intent(in) :: subgroup_j_id
      !! ID of the subgroup `j`
      real(pr) :: bij
      !! Interaction parameter \(b_{ij}\)

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      bij = self%maingroups_bij(i, j)
   end function get_subgroups_bij

   function get_subgroups_cij(self, subgroup_i_id, subgroup_j_id) result(cij)
      !! # get_subgroups_cij
      !! Get the interaction parameter \(c_{ij}\)
      !!
      !! # Description
      !! Get the interaction parameter \(c_{ij}\) of the subgroups `i` and `j`
      !! ids.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  use yaeos__models_ge_group_contribution_unifac_parameters, only: UNIFACParameters
      !!  use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      !!
      !!  type(GeGCModelParameters) :: parameters
      !!
      !!  parameters = UNIFACParameters()
      !!
      !!  ! Get the subgroups i:1, j:16 interaction parameter cij (CH3-H2O)
      !!  ! with maingroups 1 and 7 respectively.
      !!  print *, parameters%get_subgroups_cij(1, 16) ! prints: 0.0000
      !! ```
      !!
      !! In the example we obtain 0.0 because UNIFAC only have \(a_{ij} \)
      !! parameters
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.dd
      !! bst.com/published-parameters-unifac.html)
      !!
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      !! ID of the subgroup `i`
      integer, intent(in) :: subgroup_j_id
      !! ID of the subgroup `j`
      real(pr) :: cij
      !! Interaction parameter \(c_{ij}\)

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      cij = self%maingroups_cij(i, j)
   end function get_subgroups_cij
end module yaeos__models_ge_group_contribution_model_parameters
