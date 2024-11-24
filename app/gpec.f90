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
   real(pr), parameter :: z0(nc) = [1, 0] !! Component 1 molar fractions
   real(pr), parameter :: zi(nc) = [0, 1] !! Component 2 molar fractions
   real(pr) :: P !! Pressure [bar]
   real(pr) :: T !! Temperature [K]
   real(pr) :: V !! Volume [L/mol]
   real(pr) :: S

   integer :: diagram_type !! Diagram type

   integer :: i

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   ! ===========================================================================
   ! Set up the model
   ! ---------------------------------------------------------------------------
   ! model = get_model_nrtl_mhv()
   sus(1) = Substance("methane")
   sus(2) = Substance("n-butane")
   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )

   ! ===========================================================================
   ! Calculate both saturation curves
   ! ---------------------------------------------------------------------------
   print *, model%components%Tc
   print *, model%components%Pc
   do i=0,1
      sat_point = critical_point(model, z0, zi, S=S, spec=CPSpec%a, max_iters=1000, a0=S)

      print *, sat_point%iters, sat_point
      psats(i+1) = pt_envelope_2ph(model, z, sat_point, delta_0=-0.1_pr)
   end do
   call exit

   print *, "Psat 1"
   do i=1,size(psats(1)%points)
      print *, psats(1)%points(i)%T, psats(1)%points(i)%P
   end do

   print *, ""
   print *, ""
   
   print *, "Psat 2"
   do i=1,size(psats(2)%points)
      print *, psats(2)%points(i)%T, psats(2)%points(i)%P
   end do

   ! ===========================================================================
   ! Calculate the first critical line (2 -> 1)
   ! ---------------------------------------------------------------------------
   cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, dS0=-0.01_pr)
   do i=1,size(cl%a)
      write(1, *) cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do

   call exit
   
   ! cl = critical_line(model, a0=0.001_pr, z0=z0, zi=zi, dS0=0.001_pr)
   ! do i=1,size(cl%a)
   !    write(2, *) cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   ! end do
   
   call plot_pts([(real(i,pr)/100, i=1,99,10)])

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

   subroutine plot_pts(zs)
      real(pr), intent(in) :: zs(:)
      type(EquilibriumState) :: sat
      type(PTEnvel2) :: env
      integer :: i
      real(pr) :: z(nc)

      do i=1,size(zs)
         z = z0*zs(i) + zi*(1-zs(i))
         sat = saturation_pressure(model, z, T=200._pr, kind="bubble")
         env = pt_envelope_2ph(model, z, sat)
         write(i+10, *) env
         
         sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
         env = pt_envelope_2ph(model, z, sat)
         write(i+10, *) env
      end do
   end subroutine
end program gpec
