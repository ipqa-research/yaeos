program test_feos_api
   use yaeos, only: pr
   use yaeos__models_external_apis_feos, only: PCSAFT, pc_saft_to_str
   use json_module, only: json_file, json_value, json_core

   type(PCSAFT) :: model
   character(len=:), allocatable :: str

   type(json_core) :: json
   type(json_value), pointer :: p, inp

   model%sigma = [3.7039_pr, 1._pr]
   model%m = [1.0_pr, 5._pr]
   model%epsilon_k = [150.03_pr, 7._pr]

   str = pc_saft_to_str(model)

   print *, str
end program