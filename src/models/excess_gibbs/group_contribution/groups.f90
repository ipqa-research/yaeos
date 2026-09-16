module yaeos__models_ge_group_contribution_groups
   use yaeos__constants, only: pr
   implicit none

   type :: Groups
      !! # Groups
      !! Derived type used to represent a molecule and its UNIFAC groups.
      !!
      !! # Description
      !! Derived type used to represent a molecule and its UNIFAC groups. Is
      !! necessary to specify the subgroups ids and the subgroups on each
      !! molecule as shown in the example.
      !!
      !! # Examples
      !!
      !! ```fortran
      !!  ! Define toluene molecule groups
      !!  use yaeos, only: Groups
      !!
      !!  type(Groups) :: toluene
      !!
      !!  ! Toluene [ACH, ACCH3]
      !!  toluene%groups_ids = [9, 11] ! Subgroups ids
      !!  toluene%number_of_groups = [5, 1] ! Subgroups occurrences
      !! ```
      !!
      !! # References
      !! 1. [Dortmund Data Bank Software & Separation Technology](https://www.ddbst
      !! .com/published-parameters-unifac.html)
      integer, allocatable :: groups_ids(:)
      !! Indexes (ids) of each subgroup in the main group matrix
      integer, allocatable :: number_of_groups(:)
      !! Occurrences of each subgroup in the molecule
      real(pr) :: surface_area
      !! Molecule surface area \(q\)
      real(pr) :: volume
      !! Molecule volume \(r\)
   end type Groups
end module yaeos__models_ge_group_contribution_groups
