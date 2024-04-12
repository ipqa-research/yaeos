module test_ge_models
   use yaeos_constants, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose
   implicit none

   real(pr) :: absolute_tolerance = 1e-4_pr

contains

   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("NRTL", test_nrtl) &
         ]
   end subroutine collect_suite

   subroutine test_nrtl(error)
      use yaeos_constants, only: pr
      use fixtures_models, only: binary_NRTL_tape
      use yaeos, only: GeModel
      type(error_type), allocatable, intent(out) :: error

      class(GeModel), allocatable :: model
      integer, parameter :: n = 2
      real(pr) :: z(n), T
      real(pr) :: Ge, GeT, GeT2
      real(pr) :: Gen(n), GeTn(n), Gen2(n, n)

      real(pr) :: Ge_val, GeT_val, GeT2_val
      real(pr) :: Gen_val(n), GeTn_val(n), Gen2_val(n**2)
      Ge_val = 0.73578738104034147
      GeT_val = 6.2478507144254576E-002
      GeT2_val = -8.6660745337171748E-004
      Gen_val = [1.4672850471801369, 0.42228836346995569]
      GeTn_val = [0.13831602133231552, 2.9976713504363872E-002]
      Gen2_val = [-3.9849140290093517, 1.7078203950935533,&
                   1.7078203950935529, -0.73192306801724871 ]

      z = [0.3, 0.7]
      T = 150

      model = binary_NRTL_tape()
      call model%excess_gibbs( &
         z, T, Ge=Ge, GeT=GeT, &
         GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
         )

      call check(error, allclose([Ge], [Ge_val], absolute_tolerance))
      call check(error, allclose([GeT], [GeT_val], absolute_tolerance))
      call check(error, allclose([GeT2], [GeT2_val], absolute_tolerance))

      call check(error, allclose([GeTn], [GeTn_val], absolute_tolerance))
      call check(error, allclose([Gen2], [Gen2_val], absolute_tolerance))
   end subroutine test_nrtl
end module test_ge_models

