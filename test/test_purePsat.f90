program main
   use yaeos
   use yaeos__equilibria_boundaries_pure_saturation, only: pure_saturation_line, PurePsat, Psat => solve_point
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   use testing_aux, only: assert, test_title
   implicit none
   type(CubicEoS) :: model
   type(EquilibriumState) :: sat
   type(PurePsat) :: pt
   type(Substance) :: sus(2)

   real(pr) :: a, z(2), P, T

   real(Pr) :: Tc(2), Pc(2)
   integer :: i

   print *, test_title("Pure Psat curve")

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   sus(1) = Substance("methane")
   sus(2) = Substance("propane")

   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )


   ! ==========================================================================
   ! Test the pure saturation line of propane
   ! -------------------------------------------------------------------------
   pt = pure_saturation_line(model, 2, 0.001_pr, 100._pr)
   call assert(abs(pt%get_P(200._pr) - 0.2068) < 0.1_Pr, "Propane Psat at 140K")
   call assert(abs(pt%get_T(10._pr) - 300.08) < 0.1_Pr, "Propane Psat at 10 bar")
   
   
   ! ==========================================================================
   ! Test the pure saturation line of methane
   ! -------------------------------------------------------------------------
   pt = pure_saturation_line(model, 1, 1._pr, 100._pr)
   call assert(abs(pt%get_P(140._pr) - 6.45) < 0.1_Pr, "Methane Psat at 140K")
   call assert(abs(pt%get_T(10._pr) - 148.970) < 0.1_Pr, "Methane Psat at 10 bar")
end program main
