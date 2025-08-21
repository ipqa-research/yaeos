program test_autodiff_apis
   use yaeos__constants, only: pr
   use fixtures_models, only: binary_PR76_hd, binary_PR76_tape
   use yaeos, only: ArModel
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   implicit none

   real(pr) :: absolute_tolerance = 1e-6_pr

   write(*, *) test_title("AUTODIFF APIs TESTS")

   call test_pr76_hd()
   call test_pr76_tape()

   write(*, *) " "

contains

   subroutine test_pr76_hd()
      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val = -9.5079006412803206
      ArV_val = 8.8347920054119555
      ArT_val = 2.5288703990562308E-002
      ArT2_val = -8.1263537125878752E-005
      ArV2_val = -16.452678772166539
      ArTV_val = -2.4354128400887243E-002
      Arn_val = [-14.760052698527073, -19.878109854494507]
      ArVn_val = [12.970821145430932, 17.944903554047062]
      ArTn_val = [4.7299605834274866E-002, 5.0647072401531434E-002]
      Arn2_val(1, :) = [-11.697744457813224, -13.516425471096190]
      Arn2_val(2, :) = [-13.516425471096190, -19.842822840192490]

      eos = binary_PR76_hd()
      z = [0.3, 0.7]
      v = 1
      T = 150
      call eos%residual_helmholtz(z, V, T, ArV=ArV)     
      call eos%residual_helmholtz(z, V, T, ArTV=ArTV, ArV2=ArV2, ArVn=ArVn)     

      call eos%residual_helmholtz(z, V, T, ArT=ArT)     
      call eos%residual_helmholtz(z, V, T, ArTV=ArTV, ArV2=ArT2, ArVn=ArTn)     

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn &
         )

      call assert(allclose([Ar], [Ar_val], absolute_tolerance), "hdPR76: Ar")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), "hdPR76: ArV")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), "hdPR76: ArT")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), "hdPR76: ArTV")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), "hdPR76: ArV2")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), "hdPR76: ArT2")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), "hdPR76: ArVn")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), "hdPR76: ArTn")

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call assert(allclose([Ar], [Ar_val], absolute_tolerance), "hdPR76: Ar (with Arn2)")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), "hdPR76: ArV (with Arn2)")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), "hdPR76: ArT (with Arn2)")
      call assert(allclose([Arn], [Arn_val], absolute_tolerance), "hdPR76: Arn")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), "hdPR76: ArTV (with Arn2)")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), "hdPR76: ArV2 (with Arn2)")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), "hdPR76: ArT2 (with Arn2)")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), "hdPR76: ArVn (with Arn2)")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), "hdPR76: ArTn (with Arn2)")
      call assert(allclose([Arn2], [Arn2_val], absolute_tolerance), "hdPR76: Arn2")

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, ArVn=ArVn, ArTn=ArTn &
         )
      call assert(allclose([Ar], [Ar_val], absolute_tolerance), "hdPR76: Ar (final)")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), "hdPR76: ArV (final)")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), "hdPR76: ArT (final)")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), "hdPR76: ArTV (final)")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), "hdPR76: ArV2 (final)")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), "hdPR76: ArT2 (final)")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), "hdPR76: ArVn (final)")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), "hdPR76: ArTn (final)")
   end subroutine test_pr76_hd

   subroutine test_pr76_tape()
      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2
      real(pr) :: z(n), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(n), ArVn(n), ArTn(n), Arn2(n, n)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(n), ArVn_val(n), ArTn_val(n), Arn2_val(n, n)

      Ar_val = -9.5079006412803206
      ArV_val = 8.8347920054119555
      ArT_val = 2.5288703990562308E-002
      ArT2_val = -8.1263537125878752E-005
      ArV2_val = -16.452678772166539
      ArTV_val = -2.4354128400887243E-002
      Arn_val = [-14.760052698527073, -19.878109854494507]
      ArVn_val = [12.970821145430932, 17.944903554047062]
      ArTn_val = [4.7299605834274866E-002, 5.0647072401531434E-002]
      Arn2_val(1, :) = [-11.697744457813224, -13.516425471096190]
      Arn2_val(2, :) = [-13.516425471096190, -19.842822840192490]

      eos = binary_PR76_tape()
      z = [0.3, 0.7]
      v = 1
      T = 150

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn &
         )

      call assert(allclose([Ar], [Ar_val], absolute_tolerance), "tapePR76: Ar")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), "tapePR76: ArV")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), "tapePR76: ArT")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), "tapePR76: ArTV")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), "tapePR76: ArV2")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), "tapePR76: ArT2")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), "tapePR76: ArVn")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), "tapePR76: ArTn")

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call assert(allclose([Ar], [Ar_val], absolute_tolerance), "tapePR76: Ar (with Arn2)")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), "tapePR76: ArV (with Arn2)")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), "tapePR76: ArT (with Arn2)")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), "tapePR76: ArTV (with Arn2)")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), "tapePR76: ArV2 (with Arn2)")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), "tapePR76: ArT2 (with Arn2)")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), "tapePR76: ArVn (with Arn2)")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), "tapePR76: ArTn (with Arn2)")
      call assert(allclose([Arn2], [Arn2_val], absolute_tolerance), "tapePR76: Arn2")
   end subroutine test_pr76_tape
end program test_autodiff_apis
