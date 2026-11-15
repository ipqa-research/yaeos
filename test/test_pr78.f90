program test_pr78
   use yaeos, only: pr
   use testing_aux, only: assert, test_title
   use auxiliar_functions, only: allclose, rel_error
   implicit none

   write(*, *) test_title("PengRobinson78")

   call test_pr78_individual_calls()
   call test_pr78_cons_mixture()
   call test_pr78_cons_pure()

contains

   subroutine test_pr78_individual_calls()
      use yaeos, only: pr, PengRobinson78, ArModel
      use yaeos__consistency, only: individual_ar_calls

      class(ArModel), allocatable :: model, model_kij
      real(pr) :: tc(4), pc(4), w(4), kij(4, 4), lij(4, 4)

      real(pr) :: n(4), t, v

      logical :: passed

      n = [1.5, 0.2, 0.7, 2.3]
      tc = [190.564, 425.12, 300.11, 320.25]
      pc = [45.99, 37.96, 39.23, 40.21]
      w = [0.0115478, 0.200164, 0.3624, 0.298]

      t = 600_pr
      v = 0.5_pr

      ! =======================================================================
      ! Without kij
      ! -----------------------------------------------------------------------
      model = PengRobinson78(tc, pc, w)

      call individual_ar_calls(model, n, v, t, passed=passed)
      call assert(passed, "PR78 Ar individual calls (no kij): individual Ar calls match full call")

      ! =======================================================================
      ! With kij
      ! -----------------------------------------------------------------------
      kij = reshape([&
         0.0_pr, 0.1_pr, 0.2_pr, 0.1_pr, &
         0.1_pr, 0.0_pr, 0.3_pr, 0.25_pr, &
         0.2_pr, 0.3_pr, 0.0_pr, 0.18_pr, &
         0.1_pr, 0.25_pr, 0.18_pr, 0.0_pr], [size(n), size(n)])

      lij = reshape([&
         0.0_pr, 0.001_pr, 0.002_pr, 0.001_pr, &
         0.001_pr, 0.0_pr, 0.003_pr, 0.0025_pr, &
         0.002_pr, 0.003_pr, 0.0_pr, 0.0018_pr, &
         0.001_pr, 0.0025_pr, 0.0018_pr, 0.0_pr], [size(n), size(n)])


      model_kij = PengRobinson78(tc, pc, w, kij, lij)

      call individual_ar_calls(model_kij, n, v, t, passed=passed)
      call assert(passed, "PR78 Ar individual calls (kij/lij): individual Ar calls match full call")
   end subroutine test_pr78_individual_calls

   subroutine test_pr78_cons_mixture()
      use yaeos, only: pr, PengRobinson78, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency

      class(ArModel), allocatable :: model, model_kij
      real(pr) :: tc(4), pc(4), w(4)

      real(pr) :: n(4), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      real(pr) :: kij(size(n), size(n)), lij(size(n), size(n))

      n = [1.5, 0.2, 0.7, 2.3]
      tc = [190.564, 425.12, 300.11, 320.25]
      pc = [45.99, 37.96, 39.23, 40.21]
      w = [0.0115478, 0.200164, 0.3624, 0.298]

      t = 600_pr
      v = 0.5_pr

      ! ========================================================================
      ! Model without kij and lij
      ! ------------------------------------------------------------------------
      model = PengRobinson78(tc, pc, w)

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call assert(rel_error(Ar, Ar_num) < 1e-6, "PR78 Ar consistency mix (no kij): Ar vs numeric derivative")
      call assert(rel_error(ArV, ArV_num) < 1e-6, "PR78 Ar consistency mix (no kij): ArV vs numeric derivative")
      call assert(rel_error(ArT, ArT_num) < 1e-6, "PR78 Ar consistency mix (no kij): ArT vs numeric derivative")
      call assert(allclose(Arn, Arn_num, 1e-6_pr), "PR78 Ar consistency mix (no kij): Arn vs numeric derivative")
      call assert(rel_error(ArV2, ArV2_num) < 1e-6, "PR78 Ar consistency mix (no kij): ArV2 vs numeric derivative")
      call assert(rel_error(ArT2, ArT2_num) < 1e-5, "PR78 Ar consistency mix (no kij): ArT2 vs numeric derivative")
      call assert(rel_error(ArTV, ArTV_num) < 1e-6, "PR78 Ar consistency mix (no kij): ArTV vs numeric derivative")
      call assert(allclose(ArVn, ArVn_num, 1e-6_pr), "PR78 Ar consistency mix (no kij): ArVn vs numeric derivative")
      call assert(allclose(ArTn, ArTn_num, 1e-6_pr), "PR78 Ar consistency mix (no kij): ArTn vs numeric derivative")
      call assert(maxval(rel_error(Arn2, Arn2_num)) < 1e-6, &
         "PR78 Ar consistency mix (no kij): Arn2 vs numeric derivative")

      ! Consistency tests
      call assert(abs(eq31) <= 1e-14, "PR78 Ar consistency mix (no kij): eq31 (sum n_i lnphi_i = Gr/RT)")
      call assert(maxval(abs(eq33)) < 1e-15, "PR78 Ar consistency mix (no kij): eq33 (symmetric dlnphi/dn)")
      call assert(maxval(abs(eq34)) < 1e-15, "PR78 Ar consistency mix (no kij): eq34 (sum_i n_i dlnphi_i/dn_j = 0)")
      call assert(abs(eq36) <= 1e-15, "PR78 Ar consistency mix (no kij): eq36 (sum n_i dlnphi_i/dP = (Z-1)n/P)")
      call assert(abs(eq37) <= 1e-15, "PR78 Ar consistency mix (no kij): eq37 (sum n_i dlnphi_i/dT = -Hr/RT^2)")

      ! ========================================================================
      ! Model with kij and lij
      ! ------------------------------------------------------------------------
      kij = reshape([&
         0.0_pr, 0.1_pr, 0.2_pr, 0.1_pr, &
         0.1_pr, 0.0_pr, 0.3_pr, 0.25_pr, &
         0.2_pr, 0.3_pr, 0.0_pr, 0.18_pr, &
         0.1_pr, 0.25_pr, 0.18_pr, 0.0_pr], [size(n), size(n)])

      lij = reshape([&
         0.0_pr, 0.001_pr, 0.002_pr, 0.001_pr, &
         0.001_pr, 0.0_pr, 0.003_pr, 0.0025_pr, &
         0.002_pr, 0.003_pr, 0.0_pr, 0.0018_pr, &
         0.001_pr, 0.0025_pr, 0.0018_pr, 0.0_pr], [size(n), size(n)])


      model_kij = PengRobinson78(tc, pc, w, kij, lij)

      call model_kij%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model_kij, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model_kij, n, v, t, eq31=eq31, eq33=eq33, &
         eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call assert(rel_error(Ar, Ar_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): Ar vs numeric derivative")
      call assert(rel_error(ArV, ArV_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): ArV vs numeric derivative")
      call assert(rel_error(ArT, ArT_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): ArT vs numeric derivative")
      call assert(allclose(Arn, Arn_num, 1e-6_pr), "PR78 Ar consistency mix (kij/lij): Arn vs numeric derivative")
      call assert(rel_error(ArV2, ArV2_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): ArV2 vs numeric derivative")
      call assert(rel_error(ArT2, ArT2_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): ArT2 vs numeric derivative")
      call assert(rel_error(ArTV, ArTV_num) < 1e-6, "PR78 Ar consistency mix (kij/lij): ArTV vs numeric derivative")
      call assert(allclose(ArVn, ArVn_num, 1e-6_pr), "PR78 Ar consistency mix (kij/lij): ArVn vs numeric derivative")
      call assert(allclose(ArTn, ArTn_num, 1e-6_pr), "PR78 Ar consistency mix (kij/lij): ArTn vs numeric derivative")
      call assert(maxval(rel_error(Arn2, Arn2_num)) < 1e-5, &
         "PR78 Ar consistency mix (kij/lij): Arn2 vs numeric derivative")

      ! Consistency tests
      call assert(abs(eq31) <= 1e-14, "PR78 Ar consistency mix (kij/lij): eq31 (sum n_i lnphi_i = Gr/RT)")
      call assert(maxval(abs(eq33)) < 1e-15, "PR78 Ar consistency mix (kij/lij): eq33 (symmetric dlnphi/dn)")
      call assert(maxval(abs(eq34)) < 1e-14, "PR78 Ar consistency mix (kij/lij): eq34 (sum_i n_i dlnphi_i/dn_j = 0)")
      call assert(abs(eq36) <= 1e-15, "PR78 Ar consistency mix (kij/lij): eq36 (sum n_i dlnphi_i/dP = (Z-1)n/P)")
      call assert(abs(eq37) <= 1e-15, "PR78 Ar consistency mix (kij/lij): eq37 (sum n_i dlnphi_i/dT = -Hr/RT^2)")
   end subroutine test_pr78_cons_mixture

   subroutine test_pr78_cons_pure()
      use yaeos, only: pr, PengRobinson78, ArModel
      use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency

      class(ArModel), allocatable :: model
      real(pr) :: tc(1), pc(1), w(1)

      real(pr) :: n(1), t, v

      real(pr) :: Ar, ArV, ArT, Arn(size(n)), ArV2, ArT2, ArTV
      real(pr) :: ArVn(size(n)), ArTn(size(n)), Arn2(size(n), size(n))

      real(pr) :: Ar_num, ArV_num, ArT_num, Arn_num(size(n)), ArV2_num, ArT2_num
      real(pr) :: ArTV_num, ArVn_num(size(n)), ArTn_num(size(n))
      real(pr) :: Arn2_num(size(n), size(n))
      real(pr) :: eq31, eq33(size(n), size(n)), eq34(size(n)), eq36, eq37

      n = [5.0]
      tc = [190.564]
      pc = [45.99]
      w = [0.0115478]

      model = PengRobinson78(tc, pc, w)

      t = 600_pr
      v = 0.5_pr

      call model%residual_helmholtz(&
         n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, &
         ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

      call numeric_ar_derivatives(&
         model, n, v, t, d_n = 0.0001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
         Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
         ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, &
         Arn2=Arn2_num &
         )

      call ar_consistency(&
         model, n, v, t, eq31=eq31, eq33=eq33, &
         eq34=eq34, eq36=eq36, eq37=eq37 &
         )

      ! Numeric derivatives
      call assert(rel_error(Ar, Ar_num) < 1e-6, "PR78 Ar consistency pure: Ar vs numeric derivative")
      call assert(rel_error(ArV, ArV_num) < 1e-6, "PR78 Ar consistency pure: ArV vs numeric derivative")
      call assert(rel_error(ArT, ArT_num) < 1e-6, "PR78 Ar consistency pure: ArT vs numeric derivative")
      call assert(allclose(Arn, Arn_num, 1e-6_pr), "PR78 Ar consistency pure: Arn vs numeric derivative")
      call assert(rel_error(ArV2, ArV2_num) < 1e-6, "PR78 Ar consistency pure: ArV2 vs numeric derivative")
      call assert(rel_error(ArT2, ArT2_num) < 1e-6, "PR78 Ar consistency pure: ArT2 vs numeric derivative")
      call assert(rel_error(ArTV, ArTV_num) < 1e-6, "PR78 Ar consistency pure: ArTV vs numeric derivative")
      call assert(allclose(ArVn, ArVn_num, 1e-6_pr), "PR78 Ar consistency pure: ArVn vs numeric derivative")
      call assert(allclose(ArTn, ArTn_num, 1e-6_pr), "PR78 Ar consistency pure: ArTn vs numeric derivative")
      call assert(maxval(rel_error(Arn2, Arn2_num)) < 1e-5, "PR78 Ar consistency pure: Arn2 vs numeric derivative")

      ! Consistency tests
      call assert(abs(eq31) <= 1e-14, "PR78 Ar consistency pure: eq31 (sum n_i lnphi_i = Gr/RT)")
      call assert(maxval(abs(eq33)) < 1e-15, "PR78 Ar consistency pure: eq33 (symmetric dlnphi/dn)")
      call assert(maxval(abs(eq34)) < 1e-15, "PR78 Ar consistency pure: eq34 (sum_i n_i dlnphi_i/dn_j = 0)")
      call assert(abs(eq36) <= 1e-15, "PR78 Ar consistency pure: eq36 (sum n_i dlnphi_i/dP = (Z-1)n/P)")
      call assert(abs(eq37) <= 1e-15, "PR78 Ar consistency pure: eq37 (sum n_i dlnphi_i/dT = -Hr/RT^2)")
   end subroutine test_pr78_cons_pure

end program test_pr78

