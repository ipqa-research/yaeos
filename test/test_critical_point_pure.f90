program test_critical_point_pure
   use yaeos
   use fixtures_models, only: binary_PR76
   use yaeos__critical_pure_point_solver, only: find_critical_points_all_components
   use testing_aux, only: test_title, assert
   use auxiliar_functions, only: allclose
   implicit none
   class(ArModel), allocatable :: model
   real(pr) :: Tc(2), Pc(2), w(2), Vc(2)
   logical :: converged(2)

   write (*, *) test_title("Pure CP solver")

   model = binary_PR76()

   call find_critical_points_all_components(model, 2, Vc, Tc, Pc, converged)

   call assert(allclose(model%components%Pc, Pc, 1e-5_pr), "Pc")
   call assert(allclose(model%components%Tc, Tc, 1e-5_pr), "Tc")

end program