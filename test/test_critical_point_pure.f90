program test_critical_point_pure
   use yaeos
   use yaeos__extra_fluids, only: FixtureFluid, oil_b71
   use yaeos__equilibria_critical, only: get_critical_constants
   use yaeos__critical_pure_point_solver, only: find_critical_points_all_components
   use testing_aux, only: test_title, assert
   use auxiliar_functions, only: allclose
   implicit none
   class(ArModel), allocatable :: model
   type(FixtureFluid) :: fluid

   real(pr), allocatable :: Tc(:), Pc(:), w(:), Vc(:)
   logical :: converged(2)

   write (*, *) test_title("Pure CP solver")
   fluid = oil_b71()
   model = fluid%ar_model 

   allocate(Tc(fluid%nc), Pc(fluid%nc), Vc(fluid%nc), w(fluid%nc))
   call find_critical_points_all_components(model, fluid%nc, Vc, Tc, Pc, converged)

   call assert(allclose(model%components%Pc, Pc, 1e-5_pr), "Pc")
   call assert(allclose(model%components%Tc, Tc, 1e-5_pr), "Tc")

end program