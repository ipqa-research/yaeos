module yaeos__models_ge_group_contribution_model_parameters
   use yaeos__constants, only: pr
   implicit none

   type :: GeGCModelParameters
      integer, allocatable :: subgroups_ids(:)
      integer, allocatable :: maingroups_ids(:)
      integer, allocatable :: subgroups_maingroups(:)
      real(pr), allocatable :: subgroups_Rs(:)
      real(pr), allocatable :: subgroups_Qs(:)
      real(pr), allocatable :: maingroups_aij(:,:)
      real(pr), allocatable :: maingroups_bij(:,:)
      real(pr), allocatable :: maingroups_cij(:,:)
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
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      integer :: subgroup_idx

      subgroup_idx = findloc(self%subgroups_ids, subgroup_id, dim=1)
   end function get_subgroup_index

   function get_maingroup_index(self, maingroup_id) result(maingroup_idx)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_id
      integer :: maingroup_idx

      maingroup_idx = findloc(self%maingroups_ids, maingroup_id, dim=1)
   end function get_maingroup_index

   function get_subgroup_maingroup(self, subgroup_id) result(subgroup_maingroup)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      integer :: subgroup_maingroup

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_maingroup = self%subgroups_maingroups(subgroup_idx)
   end function get_subgroup_maingroup

   function get_subgroup_R(self, subgroup_id) result(subgroup_R)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      real(pr) :: subgroup_R

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_R = self%subgroups_Rs(subgroup_idx)
   end function get_subgroup_R

   function get_subgroup_Q(self, subgroup_id) result(subgroup_Q)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_id
      real(pr) :: subgroup_Q

      integer :: subgroup_idx

      subgroup_idx = self%get_subgroup_index(subgroup_id)

      subgroup_Q = self%subgroups_Qs(subgroup_idx)
   end function get_subgroup_Q

   function get_maingroups_aij(self, maingroup_i_id, maingroup_j_id) result(aij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      integer, intent(in) :: maingroup_j_id
      real(pr) :: aij

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      aij = self%maingroups_aij(i, j)
   end function get_maingroups_aij

   function get_maingroups_bij(self, maingroup_i_id, maingroup_j_id) result(bij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      integer, intent(in) :: maingroup_j_id
      real(pr) :: bij

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      bij = self%maingroups_bij(i, j)
   end function get_maingroups_bij

   function get_maingroups_cij(self, maingroup_i_id, maingroup_j_id) result(cij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: maingroup_i_id
      integer, intent(in) :: maingroup_j_id
      real(pr) :: cij

      integer :: i, j

      i = self%get_maingroup_index(maingroup_i_id)
      j = self%get_maingroup_index(maingroup_j_id)

      cij = self%maingroups_cij(i, j)
   end function get_maingroups_cij

   function get_subgroups_aij(self, subgroup_i_id, subgroup_j_id) result(aij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      integer, intent(in) :: subgroup_j_id
      real(pr) :: aij

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      aij = self%maingroups_aij(i, j)
   end function get_subgroups_aij

   function get_subgroups_bij(self, subgroup_i_id, subgroup_j_id) result(bij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      integer, intent(in) :: subgroup_j_id
      real(pr) :: bij

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      bij = self%maingroups_bij(i, j)
   end function get_subgroups_bij

   function get_subgroups_cij(self, subgroup_i_id, subgroup_j_id) result(cij)
      class(GeGCModelParameters) :: self

      integer, intent(in) :: subgroup_i_id
      integer, intent(in) :: subgroup_j_id
      real(pr) :: cij

      integer :: mi_id, mj_id, i, j

      mi_id = self%get_subgroup_maingroup(subgroup_i_id)
      mj_id = self%get_subgroup_maingroup(subgroup_j_id)

      i = self%get_maingroup_index(mi_id)
      j = self%get_maingroup_index(mj_id)

      cij = self%maingroups_cij(i, j)
   end function get_subgroups_cij
end module yaeos__models_ge_group_contribution_model_parameters
