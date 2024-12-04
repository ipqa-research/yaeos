program gpec
   !! Implementation of the Global Phase Equilibrium (GPEC) algorithm
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   use yaeos
   implicit none

   integer, parameter :: nc = 2 !! Number of components

   type(Substance) :: sus(nc) !! Substances
   class(ArModel), allocatable :: model !! Thermodynamic model to use

   type(EquilibriumState) :: sat_point !! Saturation point
   type(EquilibriumState) :: cp !! Individual critical point
   type(PTEnvel2) :: psats(2) !! Saturation curves
   type(CriticalLine) :: cl

   real(pr) :: z(nc) !! Molar fractions
   real(pr) :: a !! Fraction between component 1 and 2
   real(pr), parameter :: eps = 1e-300
   real(pr), parameter :: z0(nc) = [1-eps, eps] !! Component 1 molar fractions
   real(pr), parameter :: zi(nc) = [eps, 1-eps] !! Component 2 molar fractions
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
   sus(2) = Substance("carbon dioxide")
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
      S = i
      z = zi*S + z0*(1-S)
      sat_point = critical_point(model, z0, zi, S=S, spec=spec_CP%a, max_iters=1000, a0=S)

      select case(i)
       case(0)
         sat_point = saturation_pressure(model, z, T=100._pr, kind="bubble")
       case(1)
         sat_point = saturation_pressure(model, z, T=100._pr, kind="bubble")
      end select
      
      psats(i+1) = pt_envelope_2ph(model, z, sat_point, delta_0=0.001_pr)
   end do

   print *, "Psat 1"
   open(unit=1, file="gpec_psat1.dat")
   do i=1,size(psats(1)%points)
      write(1, *) psats(1)%points(i)%T, psats(1)%points(i)%P
   end do
   close(1)

   print *, ""
   print *, ""

   print *, "Psat 2"
   open(unit=1, file="gpec_psat2.dat")
   do i=1,size(psats(2)%points)
      print *, psats(2)%points(i)%T, psats(2)%points(i)%P
   end do
   close(2)

   ! ===========================================================================
   ! Calculate the first critical line (2 -> 1)
   ! ---------------------------------------------------------------------------
   cl = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, ns=spec_CP%a, S=0.99_pr, dS0=-0.01_pr)
   open(unit=1, file="gpec_cl2.dat")
   do i=1,size(cl%a)
      write(1, *) cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do
   close(1)


   cl = critical_line(model, a0=0.001_pr, z0=z0, zi=zi, ns=spec_CP%a, S=0.001_pr, dS0=0.001_pr)
   open(unit=1, file="gpec_cl1.dat")
   do i=1,size(cl%a)
      write(1, *) cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do
   close(1)
   
   cp = critical_point(model, z0, zi, S=log(200._pr), spec=spec_CP%P, max_iters=1000, a0=0.5_pr)
  
   !TODO: Si inicializo con S=200 converge a otro lado, ver por qué pasa!
   cl = critical_line(model, a0=cp%x(2), z0=z0, zi=zi, ns=spec_CP%P, S=log(cp%P), dS0=-0.01_pr)
   open(unit=1, file="gpec_cl3.dat")
   do i=1,size(cl%a)
      write(1, *) cl%a(i), cl%T(i), cl%P(i), cl%V(i)
   end do
   close(1)
   print *, cp%iters, cp



   
   call exit

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
   end subroutine plot_pts
end program gpec
