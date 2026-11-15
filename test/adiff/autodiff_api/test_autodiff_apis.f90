program test_autodiff_apis
   use yaeos__constants, only: pr
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   implicit none

   real(pr) :: absolute_tolerance = 1e-6_pr

   write(*, *) test_title("Autodiff APIs")
   call test_pr76_hd()
   call test_pr76_tape()
   write(*, *) ""

contains

   subroutine test_pr76_hd()
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_PR76_hd
      use yaeos, only: ArModel

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


      call assert(allclose([Ar], [Ar_val], absolute_tolerance), &
         "hdPR76: Ar vs reference value")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), &
         "hdPR76: ArV vs reference value")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), &
         "hdPR76: ArT vs reference value")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), &
         "hdPR76: ArTV vs reference value")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), &
         "hdPR76: ArV2 vs reference value")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), &
         "hdPR76: ArT2 vs reference value")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), &
         "hdPR76: ArVn vs reference value")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), &
         "hdPR76: ArTn vs reference value")

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call assert(allclose([Ar], [Ar_val], absolute_tolerance), &
         "hdPR76: Ar vs reference value")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), &
         "hdPR76: ArV vs reference value")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), &
         "hdPR76: ArT vs reference value")
      call assert(allclose([Arn], [Arn_val], absolute_tolerance), &
         "hdPR76: Arn vs reference value")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), &
         "hdPR76: ArTV vs reference value")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), &
         "hdPR76: ArV2 vs reference value")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), &
         "hdPR76: ArT2 vs reference value")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), &
         "hdPR76: ArVn vs reference value")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), &
         "hdPR76: ArTn vs reference value")
      call assert(allclose([Arn2], [Arn2_val], absolute_tolerance), &
         "hdPR76: Arn2 vs reference value")

      call eos%residual_helmholtz( &
         z, V, T, Ar=Ar, ArV=ArV, ArV2=ArV2, ArT=ArT, ArTV=ArTV, &
         ArT2=ArT2, ArVn=ArVn, ArTn=ArTn &
         )
      call assert(allclose([Ar], [Ar_val], absolute_tolerance), &
         "hdPR76: Ar vs reference value")
      call assert(allclose([ArV], [ArV_val], absolute_tolerance), &
         "hdPR76: ArV vs reference value")
      call assert(allclose([ArT], [ArT_val], absolute_tolerance), &
         "hdPR76: ArT vs reference value")
      call assert(allclose([ArTV], [ArTV_val], absolute_tolerance), &
         "hdPR76: ArTV vs reference value")
      call assert(allclose([ArV2], [ArV2_val], absolute_tolerance), &
         "hdPR76: ArV2 vs reference value")
      call assert(allclose([ArT2], [ArT2_val], absolute_tolerance), &
         "hdPR76: ArT2 vs reference value")
      call assert(allclose([ArVn], [ArVn_val], absolute_tolerance), &
         "hdPR76: ArVn vs reference value")
      call assert(allclose([ArTn], [ArTn_val], absolute_tolerance), &
         "hdPR76: ArTn vs reference value")

   end subroutine test_pr76_hd

   subroutine test_pr76_tape()
      use yaeos__constants, only: pr
      use fixtures_models, only: binary_PR76_tape
      use yaeos, only: ArModel

      class(ArModel), allocatable :: eos
      integer, parameter :: nc = 2
      real(pr) :: z(nc), V, T
      real(pr) :: Ar, ArV, ArV2, ArT, ArTV, ArT2
      real(pr) :: Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)

      real(pr) :: Ar_val, ArV_val, ArV2_val, ArT_val, ArTV_val, ArT2_val
      real(pr) :: Arn_val(nc), ArVn_val(nc), ArTn_val(nc), Arn2_val(nc, nc)

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

      call pr76_tape_numdiff(eos)
   end subroutine test_pr76_tape

   subroutine pr76_tape_numdiff(eos)
      !! Compare the tapenade PR76 derivatives against numerical ones
      !! (formerly an internal procedure of `test_pr76_tape`).
      use yaeos, only: ArModel, individual_ar_calls, numeric_ar_derivatives
      class(ArModel), intent(in) :: eos
      integer, parameter :: nc = 2
      real(pr) :: V, T, z(nc)
      real(pr) :: Ar_num, ArV_num, ArT_num, ArV2_num, ArT2_num, ArTV_num
      real(pr), dimension(nc) :: Arn_num, ArVn_num, ArTn_num
      real(pr), dimension(nc, nc) :: Arn2_num
      real(pr) :: Ar, ArV, ArT, ArV2, ArT2, ArTV
      real(pr), dimension(nc) :: Arn, ArVn, ArTn
      real(pr), dimension(nc, nc) :: Arn2

      real(pr) :: dn = 1e-5, dv=1e-6, dt=1e-2
      real(pr), parameter :: tol=1e-2
      integer :: i, j

      logical :: individual_calls


      z = [0.3, 0.7]
      v = 1
      T = 150
      
      call numeric_ar_derivatives(eos, z, V, T, dn, dv, dt, &
         Ar_num, ArV_num, ArT_num, Arn_num, ArV2_num, ArT2_num, ArTV_num, ArVn_num, ArTn_num, Arn2_num)

      call eos%residual_helmholtz(z, V, T, &
         Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, ArV2=ArV2, ArT2=ArT2, ArTV=ArTV, &
         ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )

      call assert(allclose([Ar_num], [Ar], tol), &
         "tapePR76: Ar vs numerical derivative")
      
      call assert(allclose([ArT_num] , [ArT], tol), &
         "tapePR76: ArT vs numerical derivative")
      call assert(allclose([ArT2_num], [ArT2], tol), &
         "tapePR76: ArT2 vs numerical derivative")

      call assert(allclose([ArV_num] , [ArV], tol), &
         "tapePR76: ArV vs numerical derivative")
      call assert(allclose([ArV2_num], [ArV2], tol), &
         "tapePR76: ArV2 vs numerical derivative")

      call assert(allclose([ArTV_num], [ArTV], tol), &
         "tapePR76: ArTV vs numerical derivative")

      call assert(allclose([ArVn_num], [ArVn], tol), &
         "tapePR76: ArVn vs numerical derivative")
      call assert(allclose([ArTn_num], [ArTn], tol), &
         "tapePR76: ArTn vs numerical derivative")
   end subroutine pr76_tape_numdiff
end program test_autodiff_apis
