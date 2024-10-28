program gpec
   !! Implementation of the Global Phase Equilibrium (GPEC) algorithm
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   use yaeos
   implicit none

   integer, parameter :: nc = 2 !! Number of components

   type(Substance) :: sus(nc) !! Substances
   class(ArModel), allocatable :: model !! Thermodynamic model to use

   type(EquilibriumState) :: sat_point !! Saturation point
   type(PTEnvel2) :: psats(2) !! Saturation curves
   type(CriticalLine) :: cl

   real(pr) :: z(nc) !! Molar fractions
   real(pr) :: a !! Fraction between component 1 and 2
   real(pr), parameter :: z0(nc) = [1-epsilon(1.0), epsilon(1.0)] !! Component 1 molar fractions
   real(pr), parameter :: zi(nc) = [epsilon(1.0), 1-epsilon(1.0)] !! Component 2 molar fractions
   real(pr) :: P !! Pressure [bar]
   real(pr) :: T !! Temperature [K]
   real(pr) :: V !! Volume [L/mol]

   integer :: diagram_type !! Diagram type

   integer :: i, j

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   ! ===========================================================================
   ! Set up the model
   ! ---------------------------------------------------------------------------
   model = get_model_nrtl_mhv()
   ! sus(1) = Substance("methane")
   ! sus(2) = Substance("n-butane")
   ! model = PengRobinson76(&
   !    Tc=sus%critical%critical_temperature%value, &
   !    Pc=sus%critical%critical_pressure%value/1e5, &
   !    w=sus%critical%acentric_factor%value &
   !    )

   ! ===========================================================================
   ! Calculate both saturation curves
   ! ---------------------------------------------------------------------------
   do i=0,1
      z = i*zi + (1-i)*z0
      sat_point = saturation_temperature(model, z, P=0.1_pr, kind="bubble")
      psats(i+1) = pt_envelope_2ph(model, z, sat_point)
   end do
   print *, psats(1)
   print *, ""
   print *, ""
   print *, psats(2)

   ! ===========================================================================
   ! Calculate the first critical line (2 -> 1)
   ! ---------------------------------------------------------------------------
   cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, dS0=-0.01_pr)

   do i=1,size(cl%a)
      print *, cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do
   call exit

   if (cl%a(size(cl%a)) < 1e-3) then
      type_1_or_2 : block
         ! Search for LLV
         ! IF LLV
         diagram_type = 1
         ! ELSE
         diagram_type = 2
      end block type_1_or_2
   else

   end if

contains

type(CubicEoS) function get_model_nrtl_mhv() result(model)
      type(MHV) :: mr
      type(NRTL) :: ge
      real(pr) :: a(nc,nc), b(nc,nc), c(nc,nc)
      real(pr) :: tc(nc), pc(nc), w(nc)

      a=0; b=0; c=0

      tc = [304.21_pr, 727.0_pr]
      pc = [73.83_pr, 25.6_pr]
      w = [0.223621_pr, 0.427556_pr]

      a(1, 2) = -2.8089495558489754
      a(2, 1) = -1.8212821725264361
      b(1, 2) = 1230.0987703604858
      b(2, 1) = 313.79150482742654
      c(1, 2) = 0.49412348866462708
      c(2, 1) = 0.49412348866462708

      model = PengRobinson76(tc, pc, w)
      ge = NRTL(a, b, c)
      mr = MHV(ge=ge, b=model%b, q=-0.53_pr)
      deallocate(model%mixrule)
      model%mixrule = mr
   end function get_model_nrtl_mhv
end program gpec
