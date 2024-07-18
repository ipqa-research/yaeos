module yaeos__fitting_fit_nrtl_mhv
   use yaeos__constants, only: pr
   use yaeos__fitting, only: FittingProblem
   use yaeos__models, only: ArModel, NRTL, CubicEoS, MHV
   use forsus, only: Substance
   implicit none

   integer, parameter :: nc = 2

   type, extends(FittingProblem) :: FitMHVNRTL
      logical :: fit_nrtl = .false.
      logical :: fit_lij = .false.
   contains
      procedure :: get_model_from_X => model_from_X
   end type FitMHVNRTL

contains

   subroutine model_from_X(problem, X)
      use yaeos, only: R, RKPR, PengRobinson78, ArModel, QMR, CubicEoS
      use yaeos__models_ar_cubic_quadratic_mixing, only: RKPR_D1mix
      class(FitMHVNRTL), intent(in out) :: problem
      real(pr), intent(in) :: X(:)
      type(NRTL) :: ge

      real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc)

      a=0; b=0; c=0

      a(1, 2) = x(1)
      a(2, 1) = x(2)

      b(1, 2) = x(3)
      b(2, 1) = x(4)

      c(1, 2) = x(5)
      c(2, 1) = x(6)

      ge = NRTL(a, b, c)

      associate (model => problem%model)
         select type(model)
          class is (CubicEoS)
            associate(mr => model%mixrule)
               select type (mr)
                class is (MHV)
                  if (problem%fit_lij) mr%l(1, 2) = x(7)
                  if (problem%fit_lij) mr%l(2, 1) = x(7)
                  if (problem%fit_nrtl) mr%ge = ge
               end select
            end associate
         end select
      end associate
   end subroutine model_from_X
end module yaeos__fitting_fit_nrtl_mhv
