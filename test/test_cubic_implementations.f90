module test_cubic_implementations
   use yaeos__constants, only: pr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: allclose
   implicit none

   real(pr) :: absolute_tolerance = 1e-4_pr

contains

   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("SRK", test_srk), &
         new_unittest("PR76", test_pr76), &
         new_unittest("PR78", test_pr78), &
         new_unittest("RKPR", test_RKPR) &
         ]
   end subroutine collect_suite

   subroutine test_srk(error)
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_SRK
      use yaeos, only: ArModel
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val = -9.4849231269705072
      ArV_val = 9.0478810077323164
      ArT_val = 3.0631941020155939E-002
      ArT2_val = -1.0589478951539604E-004
      ArV2_val = -17.255504598247207
      ArTV_val = -3.0039878324119831E-002
      Arn_val = [-14.710404803520872, -20.170975369630906]
      ArVn_val = [13.488065586019152, 18.870121409429380]
      ArTn_val = [5.7833039664119255E-002, 6.1888439276263030E-002]
      Arn2_val(1, :) = [-11.980899399513874, -14.133993988331257]
      Arn2_val(2, :) = [-14.133993988331257, -20.899890419408248]

      eos = binary_SRK()
      z = [0.3, 0.7]
      v = 1
      T = 150
      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
      call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
      call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
      call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
      call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
      call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

      call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
      call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
      call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
   end subroutine test_srk

   subroutine test_pr76(error)
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_PR76
      use yaeos, only: ArModel
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val = -9.5079213387597061
      ArV_val = 8.8348105702414230
      ArT_val = 2.5288760006412853E-002
      ArT2_val = -8.1263714911056052E-005
      ArV2_val = -16.452712169871607
      ArTV_val = -2.4354181554918298E-002
      Arn_val = [-14.760083989416412, -19.878152533126190]
      ArVn_val = [12.970846906902654, 17.944940224423746]
      ArTn_val = [4.7299709855544367E-002, 5.0647183777961201E-002]
      Arn2_val(1, :) = [-11.697767407192650, -13.516452437750393]
      Arn2_val(2, :) = [-13.516452437750393, -19.842863669307611]

      eos = binary_PR76()
      z = [0.3, 0.7]
      v = 1
      T = 150
      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
      call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
      call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
      call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
      call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
      call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

      call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
      call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
      call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
   end subroutine test_pr76

   subroutine test_pr78(error)
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_PR78
      use yaeos, only: ArModel
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val = -9.5079213387597061
      ArV_val = 8.8348105702414230
      ArT_val = 2.5288760006412853E-002
      ArT2_val = -8.1263714911056052E-005
      ArV2_val = -16.452712169871607
      ArTV_val = -2.4354181554918298E-002
      Arn_val = [-14.760083989416412, -19.878152533126190]
      ArVn_val = [12.970846906902654, 17.944940224423746]
      ArTn_val = [4.7299709855544367E-002, 5.0647183777961201E-002]
      Arn2_val(1, :) = [-11.697767407192650, -13.516452437750393]
      Arn2_val(2, :) = [-13.516452437750393, -19.842863669307611]

      eos = binary_PR78()
      z = [0.3, 0.7]
      v = 1
      T = 150
      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
      call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
      call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
      call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
      call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
      call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

      call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
      call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
      call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
   end subroutine test_pr78

   subroutine test_RKPR(error)
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_RKPR
      use yaeos, only: ArModel
      type(error_type), allocatable, intent(out) :: error

      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val         = -9.6482477518604082
      ArV_val        = 8.5009523009964276
      ArT_val        = 1.3361597729919862E-002
      ArT2_val       = -1.9770475708937217E-005
      ArV2_val       = -15.089728149827863
      ArTV_val       = -1.2540180323248612E-002
      Arn_val        = [-15.232342932903979, -19.399281746656321]
      ArVn_val       = [12.154119716077535, 16.347846121425317]
      ArTn_val       = [2.0758309799901429E-002,   2.8106121715415409E-002]
      Arn2_val(1, :) = [-11.397459730966577, -12.478402585247725]
      Arn2_val(2, :) = [-12.478402585247725, -18.006179159638897]

      eos = binary_RKPR()
      z = [0.3, 0.7]
      v = 1
      T = 150
      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call check(error, allclose([Ar], [Ar_val], absolute_tolerance))
      call check(error, allclose([ArV], [ArV_val], absolute_tolerance))
      call check(error, allclose([ArT], [ArT_val], absolute_tolerance))
      call check(error, allclose([ArTV], [ArTV_val], absolute_tolerance))
      call check(error, allclose([ArV2], [ArV2_val], absolute_tolerance))
      call check(error, allclose([ArT2], [ArT2_val], absolute_tolerance))

      call check(error, allclose([ArVn], [ArVn_val], absolute_tolerance))
      call check(error, allclose([ArTn], [ArTn_val], absolute_tolerance))
      call check(error, allclose([Arn2], [Arn2_val], absolute_tolerance))
   end subroutine test_RKPR
end module test_cubic_implementations
