program main
   use yaeos, only: pr, setup_dortmund, UNIFAC
   use yaeos__models_ge_group_contribution_groups, only: Groups

   type(Groups) :: molecules(3)
   type(UNIFAC) :: model

   real(pr) :: lngamma(3), z(3), T

   molecules(1)%groups_ids = [1, 2]
   molecules(1)%number_of_groups = [2, 4]

   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   molecules(3)%groups_ids = [1, 7, 8, 78, 79]
   molecules(3)%number_of_groups = [2, 1, 1, 3, 1]

   model = setup_dortmund(molecules)

   T = 303.15_pr
   z = [2.0_pr, 5.0_pr, 3.0_pr]

   call model%ln_activity_coefficient(z, T, lngamma)

   print *, lngamma

end program main