program gpec
   !! Implementation of the Global Phase Equilibrium (GPEC) algorithm
   use forsus, only: Substance, forsus_default_dir, forsus_dir
   use yaeos
   implicit none

   integer, parameter :: nc = 2 !! Number of components

   type(Substance) :: sus(nc) !! Substances
   class(ArModel), allocatable :: model !! Thermodynamic model to use

   type(EquilibriumState) :: cp !! Individual critical point
   type(CriticalLine) :: cl21 !! Critical line 2 -> 1
   type(CriticalLine) :: cl12 !! Critical line 1 -> 2
   type(CriticalLine) :: clll !! Critical line LL
   type(PurePsat) :: psats(2) !! Pure component saturation lines

   real(pr) :: z(nc) !! Molar fractions
   real(pr) :: a !! Fraction between component 1 and 2
   real(pr), parameter :: eps = 1e-3
   real(pr), parameter :: z0(nc) = [1, 0] !! Component 1 molar fractions
   real(pr), parameter :: zi(nc) = [0, 1] !! Component 2 molar fractions
   real(pr) :: P !! Pressure [bar]
   real(pr) :: T !! Temperature [K]
   real(pr) :: V !! Volume [L/mol]
   real(pr) :: S !! Specification
   real(pr) :: HPLL_P = 1000 !! High pressure limit for LLV

   integer :: diagram_type !! Diagram type

   integer :: i

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   call system("rm gpec_*")

   ! ===========================================================================
   ! Set up the model
   ! ---------------------------------------------------------------------------
   sus(1) = Substance("methane")
   sus(2) = Substance("propane")
   model = PengRobinson76(&
      Tc=sus%critical%critical_temperature%value, &
      Pc=sus%critical%critical_pressure%value/1e5, &
      w=sus%critical%acentric_factor%value &
      )
   !model = get_model_nrtl_mhv()
   model = get_modelgerg()

   ! ===========================================================================
   ! Calculate both saturation curves
   ! ---------------------------------------------------------------------------
   psats(1) = pure_saturation_line(model, 1, 1._pr, 100._pr)
   psats(2) = pure_saturation_line(model, 2, 1._pr, 100._pr)

   open(unit=1, file="gpec_Psat1.dat")
   do i=1,size(psats(1)%T)
      write(1, *) psats(1)%T(i), psats(1)%P(i), psats(1)%Vx(i), psats(1)%Vy(i)
   end do
   close(1)

   open(unit=1, file="gpec_Psat2.dat")
   do i=1,size(psats(2)%T)
      write(1, *) psats(2)%T(i), psats(2)%P(i), psats(2)%Vx(i), psats(2)%Vy(i)
   end do
   close(1)

   ! ===========================================================================
   ! Calculate the first critical line (2 -> 1)
   ! ---------------------------------------------------------------------------
   cl21 = critical_line(model, a0=0.99_pr, z0=z0, zi=zi, ns=spec_CP%a, S=0.99_pr, dS0=-0.01_pr, maxP=HPLL_P)
   open(unit=1, file="gpec_cl21.dat")
   call write_cl(cl21)
   close(1)

   if (abs(cl21%P(size(cl21%a)) - model%components%Pc(1)) > 0.1_pr) then
      ! ========================================================================
      ! Calculate the second critical line (1 -> 2)
      ! ------------------------------------------------------------------------
      cl12 = critical_line(model, a0=0.001_pr, z0=z0, zi=zi, ns=spec_CP%a, S=0.001_pr, dS0=0.001_pr)
      open(unit=1, file="gpec_cl12.dat")
      call write_cl(cl12)
      close(1)
   else
      diagram_type = 1
   end if


   if (diagram_type == 1) then
      type_1_or_2 : block
         ! Search for LLV
         !TODO: Si inicializo con S=200 converge a otro lado, ver por qué pasa!
         ! Converge a critical point at high pressure
         cp = critical_point(model, z0, zi, S=log(HPLL_P), spec=spec_CP%P, max_iters=1000, a0=0.5_pr)
         clll = critical_line(model, a0=cp%x(2), z0=z0, zi=zi, ns=spec_CP%P, S=log(cp%P), dS0=-0.01_pr)
         open(unit=1, file="gpec_clll.dat")
         call write_cl(clll)
         close(1)
         if (size(clll%a) > 1) then
            diagram_type = 2
         else
         end if
      end block type_1_or_2
   else
      diagram_type = 3
   end if

   ! ===========================================================================
   ! Now with the global diagram information, calculate specific points
   ! ---------------------------------------------------------------------------

   call plot_pts([(real(i,pr)/100, i=1,99,10)])
   call plot_pxs([(real(i, pr), i=150,int(model%components%Tc(2)),20)])
   call plot_txs([(real(i, pr), i=0,int(model%components%Pc(1)),50)])
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

   type(GERG2008) function get_modelgerg() result(model)
      integer :: ids(2)
      ids = [G2008Components%methane, G2008Components%carbon_dioxide]
      model = gerg_2008(ids)
   end function get_modelgerg

   subroutine plot_pts(zs)
      real(pr), intent(in) :: zs(:)
      type(EquilibriumState) :: sat
      type(PTEnvel2) :: env
      integer :: i
      real(pr) :: z(nc)

      open(1, file="gpec_pts.dat")
      write(1, *) "kind T P beta x1 x2 y1 y2 z1 z2"
      do i=1,size(zs)

         z = [1-zs(i), zs(i)]

         sat = saturation_pressure(model, z, T=150._pr, kind="bubble")
         env = pt_envelope_2ph(model, z, sat, points=1000)
         call write_pt(env, z)

         sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
         env = pt_envelope_2ph(model, z, sat, points=1000)
         call write_pt(env, z)

         write(1, *)
         write(1, *)
      end do
      close(1)
   end subroutine plot_pts

   subroutine plot_pxs(Ts)
      real(pr), intent(in) :: Ts(:)
      real(pr) :: a, z(nc)
      integer :: i, j
      integer :: loc
      type(EquilibriumState) :: sat
      type(PXEnvel2) :: px


      open(unit=1, file="gpec_px.dat")
      write(1, *) "kind T P beta x1 x2 y1 y2"
      do i=1, size(Ts)
         T = Ts(i)

         if (Ts(i) < model%components%Tc(1)) then
            a = 0.0001_pr
            z = a*zi + (1-a)*z0
            P = psats(1)%get_P(T)
            sat = saturation_pressure(model, z, T, kind="bubble", P0=P)
            px = px_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat, delta_0=0.001_pr)
         else if (Ts(i) < model%components%Tc(2)) then
            a = 0.9999_pr
            z = a*zi + (1-a)*z0
            P = psats(2)%get_P(T)
            sat = saturation_pressure(model, z, T, kind="bubble", P0=P)
            px = px_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat, delta_0=-0.001_pr)
         else
            loc = minloc(abs(cl12%T - T), dim=1)
            a = cl12%a(loc) + eps
            z = a*zi + (1-a)*z0

         end if

         do j=1,size(px%points)
            write(1, *) px%points(j)
         end do
         write(1, *)
         write(1, *)
      end do
      close(1)

   end subroutine plot_pxs

   subroutine plot_txs(Ps)
      real(pr), intent(in) :: Ps(:)
      real(pr) :: a, z(nc)
      integer :: i, j
      type(EquilibriumState) :: sat
      type(TXEnvel2) :: tx

      open(unit=1, file="gpec_tx.dat")
      write(1, *) "kind T P beta x1 x2 y1 y2"

      do i=1, size(Ps)
         P = Ps(i)

         if (Ps(i) < model%components%Pc(2)) then
            a = 0.9999_pr
            z = a*zi + (1-a)*z0
            T = psats(2)%get_T(P)
            sat = saturation_temperature(model, z, P, kind="dew", T0=T)
            tx = tx_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat, delta_0=-0.0001_pr)
            call write_tx(tx)


            if (tx%alpha(size(tx%alpha)) > 0.001_pr) then
               ! In the case the line did not reach the light component
               a = 0.0001_pr
               z = a*zi + (1-a)*z0
               T = psats(1)%get_T(P)
               sat = saturation_temperature(model, z, T, kind="dew", T0=T)
               tx = tx_envelope_2ph(model, z0=z0, alpha0=a, z_injection=zi, first_point=sat, delta_0=0.0001_pr)
               call write_tx(tx)
            end if
         end if

      end do

      close(1)
   end subroutine plot_txs

   subroutine write_pt(pt, z)
      type(PTEnvel2), intent(in) :: pt
      real(pr), intent(in) :: z(nc)
      integer :: j
      do j=1,size(pt%points)
         write(1, *) pt%points(j), z
      end do
      write(1, *)
      write(1, *)
   end subroutine write_pt

   subroutine write_tx (tx)
      type(TXEnvel2), intent(in) :: tx
      integer :: j
      do j=1,size(tx%points)
         write(1, *) tx%points(j)
      end do
      write(1, *)
      write(1, *)
   end subroutine write_tx

   subroutine write_cl(cl)
      type(CriticalLine), intent(in) :: cl
      integer :: j
      do j=1,size(cl%a)
         write(1, *) cl%a(j), cl%T(j), cl%P(j), cl%V(j)
      end do
      write(1, *)
      write(1, *)
   end subroutine write_cl
end program gpec
