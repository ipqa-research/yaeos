module test_flash
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("FlashPT", test_flash_pt), &
         new_unittest("FlashTV", test_flash_tv), &
         new_unittest("FlashPT Failed", test_flash_pt_failed), &
         new_unittest("FlashPT Bad Specification", test_flash_pt_bad_spec), &
         new_unittest("Stability tm minimization", test_tm) &
         ]
   end subroutine collect_suite

   subroutine test_tm(error)
      use forsus, only: Substance, forsus_dir
      use yaeos
      use yaeos__equilibria_stability, only: tm, min_tpd
      use yaeos, only: flash
      implicit none
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc=2

      class(ArModel), allocatable :: model
      type(Substance) :: sus(nc)
      real(pr) :: tc(nc), pc(nc), ac(nc)
      real(pr) :: z(nc), T, P
      real(pr) :: w(nc), mintpd

      integer :: i

      forsus_dir = "build/dependencies/forsus/data/json"
      sus(1) = Substance("methane")
      sus(2) = Substance("hydrogen sulfide")

      z = [0.13, 1-0.13]
      z = z/sum(z)

      P = 20.0_pr
      T = 190._pr

      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5_pr
      ac = sus%critical%acentric_factor%value

      model = SoaveRedlichKwong(tc, pc, ac)

      call min_tpd(model, z, P, T, mintpd, w)
      call check(error, abs(mintpd) < 1e-8_pr)
      
      P = 15
      call min_tpd(model, z, P, T, mintpd, w)
      call check(error, abs(mintpd - (-0.1726_pr)) < 1e-3)
      call check(error, abs(tm(model, z, w, p, t) - mintpd) < 1e-10_pr)
   end subroutine test_tm

   subroutine test_flash_pt(error)
      use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      real(pr) :: x(nc) = [0.32424471950363210, 0.67575531029866709]
      real(pr) :: y(nc) = [0.91683466155334536, 8.3165368249135715E-002]
      real(pr) :: vx = 8.4918883298198036E-002
      real(pr) :: vy = 0.32922132295944545

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(pr) :: tc(nc), pc(nc), w(nc)
      real(pr) :: n(nc), p, t, k0(nc)
      integer :: iters

      n = [0.4, 0.6]
      tc = [190.564, 425.12]
      pc = [45.99, 37.96]
      w = [0.0115478, 0.200164]
      model = PengRobinson76(tc, pc, w)

      P = 60
      t = 294
      k0 = (PC/P)*exp(5.373*(1 + w)*(1 - TC/T))

      flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iters)

      call check(error, maxval(abs(flash_result%x - x)) < 1e-5)
      call check(error, maxval(abs(flash_result%y - y)) < 1e-5)
      call check(error, (abs(flash_result%Vx - Vx)) < 1e-5)
      call check(error, (abs(flash_result%Vy - Vy)) < 1e-5)
   end subroutine test_flash_pt

   subroutine test_flash_tv(error)
      use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel
      use fixtures_models, only: binary_PR76
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      real(pr) :: x(nc) = [6.8598497458814592E-002,  0.93140153234346601]
      real(pr) :: y(nc) = [0.61055654015073813, 0.38944348965161085]
      real(pr) :: vx = 9.3682339483042124E-002
      real(pr) :: vy = 2.3935064128338039
      real(pr) :: beta = 0.61148923421371815
      real(pr) :: P = 6.097517429661468

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(pr) :: tc(nc), pc(nc), w(nc)
      real(pr) :: n(nc), t, k0(nc), v
      integer :: iters, i

      n = [0.4, 0.6]
      model = binary_PR76()

      V = 1.5_pr
      t = 210
      k0 = (model%components%Pc/10._pr) &
         * exp(5.373_pr*(1 + model%components%w)&
         * (1 - model%components%Tc/T))

      flash_result = flash(model, n, t=t, v_spec=v, k0=k0, iters=iters)

      call check(error, maxval(abs(flash_result%x - x)) < 1e-5)
      call check(error, maxval(abs(flash_result%y - y)) < 1e-5)
      call check(error, abs(flash_result%vx - vx)  < 1e-5)
      call check(error, abs(flash_result%vy - vy)  < 1e-5)
      call check(error, abs(flash_result%p - p)  < 1e-5)
      call check(error, abs(flash_result%beta - beta)  < 1e-5)
   end subroutine test_flash_tv

   subroutine test_flash_pt_failed(error)
      use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      real(pr) :: x(nc) = [0.32424471950363210, 0.67575531029866709]
      real(pr) :: y(nc) = [0.91683466155334536, 8.3165368249135715E-002]
      real(pr) :: vx = 8.4918883298198036E-002
      real(pr) :: vy = 0.32922132295944545

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(pr) :: tc(nc), pc(nc), w(nc)
      real(pr) :: n(nc), p, t, k0(nc)
      integer :: iters

      n = [0.4, 0.6]
      tc = [190.564, 425.12]
      pc = [45.99, 37.96]
      w = [0.0115478, 0.200164]
      model = PengRobinson76(tc, pc, w)

      P = 500
      t = 294
      k0 = (PC/P)*exp(5.373*(1 + w)*(1 - TC/T))

      flash_result = flash(model, n, t=t, p_spec=p, k0=k0, iters=iters)
      call check(error, flash_result%p < 0)
   end subroutine test_flash_pt_failed

   subroutine test_flash_pt_bad_spec(error)
      use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nc = 2

      real(pr) :: x(nc) = [0.32424471950363210, 0.67575531029866709]
      real(pr) :: y(nc) = [0.91683466155334536, 8.3165368249135715E-002]
      real(pr) :: vx = 8.4918883298198036E-002
      real(pr) :: vy = 0.32922132295944545

      class(ArModel), allocatable :: model
      type(EquilibriumState) :: flash_result

      real(pr) :: tc(nc), pc(nc), w(nc)
      real(pr) :: n(nc), p, t, k0(nc)
      integer :: iters

      n = [0.4, 0.6]
      tc = [190.564, 425.12]
      pc = [45.99, 37.96]
      w = [0.0115478, 0.200164]
      model = PengRobinson76(tc, pc, w)

      P = 500
      t = 294
      k0 = (PC/P)*exp(5.373*(1 + w)*(1 - TC/T))

      flash_result = flash(&
         model, n, t=t, p_spec=p, v_spec=2._pr, k0=k0, iters=iters&
         )
   end subroutine test_flash_pt_bad_spec
end module test_flash
