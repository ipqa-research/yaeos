module yaeos__fitting_fit_nrtl_mhv
   use yaeos__constants, only: pr
   use yaeos__fitting, only: FittingProblem
   use yaeos__models, only: ArModel, NRTL, CubicEoS, MHV
   use forsus, only: Substance
   implicit none

   integer, parameter :: nc = 2

   type, extends(FittingProblem) :: FitMHVNRTL
   contains
      procedure :: get_model_from_X => model_from_X
   end type FitMHVNRTL

contains

   subroutine init_model(problem, sus)
      use yaeos, only: R, ArModel, CubicEoS, PengRobinson78, RKPR, SoaveRedlichKwong
      class(FitMHVNRTL), intent(in out) :: problem
      type(Substance), intent(in) :: sus(2)
      type(MHV) :: mixrule
      type(NRTL) :: ge
      real(pr) :: tc(nc), pc(nc), w(nc), vc(nc), zc(nc)
      real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc), bi(nc)

      a=0; b=0; c=0

      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5
      w = sus%critical%acentric_factor%value
      vc = sus%critical%critical_volume%value
      zc = pc*vc/(R*tc)

      ge = NRTL(a, b, c)

      allocate(CubicEoS :: problem%model)
      problem%model = SoaveRedlichKwong(tc, pc, w)

      associate(m => problem%model)
         select type(m)
          type is (CubicEoS)
            bi = m%b
            mixrule = MHV(ge=ge, q=-0.593_pr, b=bi)
            deallocate(m%mixrule)
            m%mixrule = mixrule
         end select
      end associate
   end subroutine init_model

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
                  mr%l(1, 2) = x(7)
                  mr%l(2, 1) = x(7)
                  mr%ge = ge
                  model%del1 = x(8:)
                  model%del2 = (1._pr - model%del1)/(1._pr + model%del1)
               end select
            end associate
         end select
      end associate
   end subroutine model_from_X
end module yaeos__fitting_fit_nrtl_mhv
