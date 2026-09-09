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

   real(Pr) :: Tc(2), Pc(2), w(2)
   integer :: i

   write(*, *) test_title("Pure Psat curve")

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


   ! ==========================================================================
   ! Propane-Heptane mixture
   ! -------------------------------------------------------------------------
   Tc = [369.83, 540.2]
   Pc = [42.48, 27.4]
   w = [0.152291, 0.349469]

   model = PengRobinson76(Tc, Pc, w)
   pt = pure_saturation_line(model, 1, 0.001_pr, 100._pr)
   call assert(pt%T(1) < 141._Pr, "Reach to the lowest temperature of propane saturation line")
   
   pt = pure_saturation_line(model, 2, 0.001_pr, 100._pr)
   call assert(pt%P(1) < 1.e-3_pr, "Reach to the lowest pressure of Heptane saturation line")

   call test_gerg2008_nitrogen

contains
   subroutine test_gerg2008_nitrogen
      type(GERG2008) :: eos
      eos = gerg_2008([G2008Components%nitrogen])
      pt = pure_saturation_line(eos, 1, 0.001_pr, 100._pr)
      call assert(pt%T(1) < 101._Pr, "Reach to the lowest temperature of Nitrogen saturation line")
   end subroutine
end program main
