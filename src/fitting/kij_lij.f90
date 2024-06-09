module yaeos__fitting_fit_kij_lij
   !! Binary interaction parameters fitting problem
   use forsus, only: Substance, forsus_dir
   use yaeos__fitting, only: pr, FittingProblem, ArModel
   implicit none

   integer, parameter :: nc=2
   type, extends(FittingProblem) :: FitKijLij
      logical :: fit_lij = .false.
      logical :: fit_kij = .false.
   contains
      procedure :: get_model_from_X => model_from_X
   end type FitKijLij

contains

   subroutine init_model(problem, sus)
      use yaeos, only: R, ArModel, CubicEoS, PengRobinson78, RKPR, SoaveRedlichKwong
      class(FitKijLij), intent(in out) :: problem
      type(Substance), intent(in) :: sus(2)
      real(pr) :: tc(nc), pc(nc), w(nc), vc(nc), zc(nc)

      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5
      w = sus%critical%acentric_factor%value
      vc = sus%critical%critical_volume%value
      zc = pc * vc / (R * tc)

      problem%model = PengRobinson78(tc, pc, w)
   end subroutine init_model

   function model_from_X(problem, X) result(model)
      use yaeos, only: R, RKPR, PengRobinson78, ArModel, QMR, CubicEoS
      real(pr), intent(in) :: X(:)
      class(FitKijLij), intent(in) :: problem
      class(ArModel), allocatable :: model

      real(pr) :: kij(nc, nc), lij(nc, nc)

      kij = 0
      kij(1, 2) = X(1)
      kij(2, 1) = kij(1, 2)

      lij = 0
      lij(1, 2) = X(2)
      lij(2, 1) = X(2)

      model = problem%model

      select type (model)
       class is (CubicEoS)
         associate (mr => model%mixrule)
            select type(mr)
             class is (QMR)
               if (problem%fit_kij) mr%k = kij
               if (problem%fit_lij) mr%l = lij
            end select
         end associate
      end select
   end function model_from_X
end module yaeos__fitting_fit_kij_lij

