module test_consistency_checkers_models
   !! Analytic mock models used to check that the consistency utilities
   !! detect inconsistencies between "all at once" and individual calls.
   !!
   !! Both models use a quadratic form Q(n) = n^T K n with symmetric K:
   !!   Ar(n, V, T) = Q(n) / (V T)
   !!   Ge(n, T)    = Q(n) / T
   !! When `corrupt` = k > 0 the k-th derivative returned by a call that
   !! requests every output is shifted by `shift`, while individual calls stay
   !! exact.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use yaeos__models_ge, only: GeModel
   implicit none

   real(pr), parameter :: shift = 1.0e-3_pr
   integer :: corrupt = 0

   type, extends(ArModel) :: MockAr
      real(pr) :: K(2, 2)
   contains
      procedure :: residual_helmholtz => mock_ar
      procedure :: get_v0 => mock_v0
   end type MockAr

   type, extends(GeModel) :: MockGe
      real(pr) :: K(2, 2)
   contains
      procedure :: excess_gibbs => mock_ge
   end type MockGe

contains

   subroutine mock_ar(self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
      class(MockAr), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: v
      real(pr), intent(in) :: t
      real(pr), optional, intent(out) :: Ar
      real(pr), optional, intent(out) :: ArV
      real(pr), optional, intent(out) :: ArT
      real(pr), optional, intent(out) :: ArT2
      real(pr), optional, intent(out) :: ArTV
      real(pr), optional, intent(out) :: ArV2
      real(pr), optional, intent(out) :: Arn(size(n))
      real(pr), optional, intent(out) :: ArVn(size(n))
      real(pr), optional, intent(out) :: ArTn(size(n))
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))

      real(pr) :: Q
      real(pr) :: Kn(size(n))
      real(pr) :: d(10)
      logical :: all_present

      all_present = present(Ar) .and. present(ArV) .and. present(ArT) &
         .and. present(ArT2) .and. present(ArTV) .and. present(ArV2) &
         .and. present(Arn) .and. present(ArVn) .and. present(ArTn) .and. present(Arn2)
      d = 0.0_pr
      if (all_present .and. corrupt > 0) d(corrupt) = shift

      Kn = matmul(self%K, n)
      Q = dot_product(n, Kn)

      if (present(Ar)) Ar = Q / (v*t) + d(1)
      if (present(ArV)) ArV = -Q / (v**2 * t) + d(2)
      if (present(ArT)) ArT = -Q / (v * t**2) + d(3)
      if (present(Arn)) Arn = 2 * Kn / (v*t) + d(4)
      if (present(ArV2)) ArV2 = 2 * Q / (v**3 * t) + d(5)
      if (present(ArT2)) ArT2 = 2 * Q / (v * t**3) + d(6)
      if (present(ArTV)) ArTV = Q / (v**2 * t**2) + d(7)
      if (present(ArVn)) ArVn = -2 * Kn / (v**2 * t) + d(8)
      if (present(ArTn)) ArTn = -2 * Kn / (v * t**2) + d(9)
      if (present(Arn2)) Arn2 = 2 * self%K / (v*t) + d(10)
   end subroutine mock_ar

   function mock_v0(self, n, p, t)
      class(MockAr), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: mock_v0
      mock_v0 = 0.1_pr * sum(n)
   end function mock_v0

   subroutine mock_ge(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      class(MockGe), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: T
      real(pr), optional, intent(out) :: Ge
      real(pr), optional, intent(out) :: GeT
      real(pr), optional, intent(out) :: GeT2
      real(pr), optional, intent(out) :: Gen(size(n))
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      real(pr) :: Q
      real(pr) :: Kn(size(n))
      real(pr) :: d(6)
      logical :: all_present

      all_present = present(Ge) .and. present(GeT) .and. present(GeT2) &
         .and. present(Gen) .and. present(GeTn) .and. present(Gen2)
      d = 0.0_pr
      if (all_present .and. corrupt > 0) d(corrupt) = shift

      Kn = matmul(self%K, n)
      Q = dot_product(n, Kn)

      if (present(Ge)) Ge = Q / T + d(1)
      if (present(GeT)) GeT = -Q / T**2 + d(2)
      if (present(GeT2)) GeT2 = 2 * Q / T**3 + d(3)
      if (present(Gen)) Gen = 2 * Kn / T + d(4)
      if (present(GeTn)) GeTn = -2 * Kn / T**2 + d(5)
      if (present(Gen2)) Gen2 = 2 * self%K / T + d(6)
   end subroutine mock_ge
end module test_consistency_checkers_models

program test_consistency_checkers
   use yaeos__constants, only: pr, R
   use yaeos__consistency_armodel, only: individual_ar_calls, numeric_ar_derivatives
   use yaeos__consistency_gemodel, only: individual_ge_calls, numeric_ge_derivatives, ge_consistency
   use test_consistency_checkers_models
   use testing_aux, only: assert, test_title
   implicit none

   character(len=*), parameter :: ar_names(10) = [character(len=4) :: &
      "Ar", "ArV", "ArT", "Arn", "ArV2", "ArT2", "ArTV", "ArVn", "ArTn", "Arn2"]
   character(len=*), parameter :: ge_names(6) = [character(len=4) :: &
      "Ge", "GeT", "GeT2", "Gen", "GeTn", "Gen2"]

   type(MockAr) :: ar
   type(MockGe) :: ge
   real(pr) :: n(2)
   real(pr) :: V
   real(pr) :: T
   logical :: passed
   integer :: k

   write(*, *) test_title("Consistency checkers on analytic mock models")

   ar%K = reshape([1.0_pr, 0.3_pr, 0.3_pr, 2.0_pr], [2, 2])
   ge%K = ar%K
   n = [0.4_pr, 0.6_pr]
   V = 2.0_pr
   T = 3.0_pr

   ! ---- ArModel -------------------------------------------------------------
   corrupt = 0
   call individual_ar_calls(ar, n, V, T, passed)
   call assert(passed, "individual_ar_calls: passes for a consistent model")

   do k = 1, 10
      corrupt = k
      call individual_ar_calls(ar, n, V, T, passed)
      call assert(.not. passed, "individual_ar_calls: detects a mismatch in " // trim(ar_names(k)))
   end do
   corrupt = 0

   call check_numeric_ar()

   ! ---- GeModel -------------------------------------------------------------
   corrupt = 0
   call individual_ge_calls(ge, n, T, passed)
   call assert(passed, "individual_ge_calls: passes for a consistent model")

   do k = 1, 6
      corrupt = k
      call individual_ge_calls(ge, n, T, passed)
      call assert(.not. passed, "individual_ge_calls: detects a mismatch in " // trim(ge_names(k)))
   end do
   corrupt = 0

   call check_numeric_ge()

   write(*, *) ""

contains

   subroutine check_numeric_ar()
      !! Central differences of the mock Ar against its analytic derivatives.
      real(pr) :: Ar_num
      real(pr) :: ArV_num
      real(pr) :: ArT_num
      real(pr) :: Arn_num(2)
      real(pr) :: ArV2_num
      real(pr) :: ArT2_num
      real(pr) :: ArTV_num
      real(pr) :: ArVn_num(2)
      real(pr) :: ArTn_num(2)
      real(pr) :: Arn2_num(2, 2)
      real(pr) :: Ar_an
      real(pr) :: ArV_an
      real(pr) :: ArT_an
      real(pr) :: Arn_an(2)
      real(pr) :: ArV2_an
      real(pr) :: ArT2_an
      real(pr) :: ArTV_an
      real(pr) :: ArVn_an(2)
      real(pr) :: ArTn_an(2)
      real(pr) :: Arn2_an(2, 2)
      real(pr), parameter :: tol = 1e-5_pr

      call ar%residual_helmholtz(n, V, T, Ar=Ar_an, ArV=ArV_an, ArT=ArT_an, &
         ArTV=ArTV_an, ArV2=ArV2_an, ArT2=ArT2_an, Arn=Arn_an, ArVn=ArVn_an, &
         ArTn=ArTn_an, Arn2=Arn2_an)
      call numeric_ar_derivatives(ar, n, V, T, d_n=1e-4_pr, d_v=1e-4_pr, d_t=1e-4_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, Arn=Arn_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, ArTV=ArTV_num, ArVn=ArVn_num, ArTn=ArTn_num, Arn2=Arn2_num)

      call assert(abs(Ar_num - Ar_an) < 1e-14_pr, "numeric_ar_derivatives: Ar at the point")
      call assert(abs(ArV_num - ArV_an) < tol, "numeric_ar_derivatives: ArV")
      call assert(abs(ArT_num - ArT_an) < tol, "numeric_ar_derivatives: ArT")
      call assert(maxval(abs(Arn_num - Arn_an)) < tol, "numeric_ar_derivatives: Arn")
      call assert(abs(ArV2_num - ArV2_an) < tol, "numeric_ar_derivatives: ArV2")
      call assert(abs(ArT2_num - ArT2_an) < tol, "numeric_ar_derivatives: ArT2")
      call assert(abs(ArTV_num - ArTV_an) < tol, "numeric_ar_derivatives: ArTV")
      call assert(maxval(abs(ArVn_num - ArVn_an)) < tol, "numeric_ar_derivatives: ArVn")
      call assert(maxval(abs(ArTn_num - ArTn_an)) < tol, "numeric_ar_derivatives: ArTn")
      call assert(maxval(abs(Arn2_num - Arn2_an)) < tol, "numeric_ar_derivatives: Arn2")

      ! Only second derivatives requested: the first-derivative branch is skipped
      ArV2_num = 0.0_pr
      ArT2_num = 0.0_pr
      call numeric_ar_derivatives(ar, n, V, T, d_n=1e-4_pr, d_v=1e-4_pr, d_t=1e-4_pr, &
         Ar=Ar_num, ArV2=ArV2_num, ArT2=ArT2_num)
      call assert(abs(ArV2_num - ArV2_an) < tol, "numeric_ar_derivatives: ArV2 without ArV")
      call assert(abs(ArT2_num - ArT2_an) < tol, "numeric_ar_derivatives: ArT2 without ArT")
   end subroutine check_numeric_ar

   subroutine check_numeric_ge()
      !! Central differences of the mock Ge against its analytic derivatives,
      !! and the Michelsen-Mollerup identities (Eqs. 60-61) on the mock.
      real(pr) :: Ge_num
      real(pr) :: GeT_num
      real(pr) :: Gen_num(2)
      real(pr) :: GeT2_num
      real(pr) :: GeTn_num(2)
      real(pr) :: Gen2_num(2, 2)
      real(pr) :: Ge_an
      real(pr) :: GeT_an
      real(pr) :: Gen_an(2)
      real(pr) :: GeT2_an
      real(pr) :: GeTn_an(2)
      real(pr) :: Gen2_an(2, 2)
      real(pr) :: eq60(2, 2)
      real(pr) :: eq61(2)
      real(pr), parameter :: tol = 1e-5_pr

      call ge%excess_gibbs(n, T, Ge=Ge_an, GeT=GeT_an, GeT2=GeT2_an, Gen=Gen_an, &
         GeTn=GeTn_an, Gen2=Gen2_an)
      call numeric_ge_derivatives(ge, n, T, d_n=1e-4_pr, d_t=1e-4_pr, Ge=Ge_num, &
         GeT=GeT_num, Gen=Gen_num, GeT2=GeT2_num, GeTn=GeTn_num, Gen2=Gen2_num)

      call assert(abs(Ge_num - Ge_an) < 1e-14_pr, "numeric_ge_derivatives: Ge at the point")
      call assert(abs(GeT_num - GeT_an) < tol, "numeric_ge_derivatives: GeT")
      call assert(abs(GeT2_num - GeT2_an) < tol, "numeric_ge_derivatives: GeT2")
      call assert(maxval(abs(Gen_num - Gen_an)) < tol, "numeric_ge_derivatives: Gen")
      call assert(maxval(abs(GeTn_num - GeTn_an)) < tol, "numeric_ge_derivatives: GeTn")
      call assert(maxval(abs(Gen2_num - Gen2_an)) < tol, "numeric_ge_derivatives: Gen2")

      ! Eq. 60 (Hessian symmetry) holds for a symmetric K. Eq. 61 (Gibbs-Duhem
      ! on Gen2) is n^T Gen2 = 2 K n / T, nonzero for this non-homogeneous-of-
      ! degree-one mock: the checker must report exactly that value.
      call ge_consistency(ge, n, T, eq60=eq60, eq61=eq61)
      call assert(maxval(abs(eq60)) < 1e-14_pr, "ge_consistency: Eq. 60 symmetry is zero")
      call assert(maxval(abs(eq61 - matmul(n, 2*ge%K/T) / (R*T))) < 1e-12_pr, &
         "ge_consistency: Eq. 61 equals n^T Gen2 / (R T)")
   end subroutine check_numeric_ge
end program test_consistency_checkers
