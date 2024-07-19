module yaeos__HardSpheresCubicEoS
   use yaeos, only: pr, R, Substances
   use yaeos__ar_models_hyperdual, only: ArModelAdiff
   use yaeos__autodiff
   implicit none

   type, extends(ArModelAdiff) :: HardSpheresCubicEoS
      real(pr), allocatable :: ac(:), b(:), k(:), kij(:, :), lij(:, :)
   contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
   end type HardSpheresCubicEoS

   real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
   real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)

contains

   function arfun(self, n, v, t) result(Ar)
      class(HardSpheresCubicEoS) :: self !! Model
      type(hyperdual), intent(in) :: n(:) !! Number of moles vector
      type(hyperdual), intent(in) :: v !! Volume [L/mol]
      type(hyperdual), intent(in) :: t !! Temperature [K]
      type(hyperdual) :: Ar !! Residual Helmholtz energy

      type(hyperdual), dimension(size(n)) :: a
      type(hyperdual) :: nij

      ! Cubic Parameters
      type(hyperdual) :: Ar_att
      type(hyperdual) :: amix, bmix
      type(hyperdual) :: b_v

      ! HardMixing parameters
      type(hyperdual) :: Ar_rep
      type(hyperdual) :: lambda(0:3), eta
      type(hyperdual) :: l1l2_l3l0
      type(hyperdual) :: l23_l0l32
      type(hyperdual) :: logContribution
      real(pr), parameter    :: xi = 4.0_pr

      integer :: i, j
      integer :: nc

      nc = size(n)

      ! Alpha function
      associate(&
         ac => self%ac, b => self%b, k => self%k, &
         tc => self%components%tc, &
         kij => self%kij, lij => self%lij)

         ! Attractive parameter
         a = self%ac*(1.0_pr + self%k*(1.0_pr - sqrt(t/self%components%tc)))**2

         ! Mixing rule
         amix = 0.0_pr
         bmix = 0.0_pr
         lambda = 0.0_pr
         do i = 1, nc
            do j = 1, nc
               nij = n(i)*n(j)
               amix = amix + nij * sqrt(a(i)*a(j)) * (1 - kij(i, j))
               bmix = bmix + nij * 0.5_pr * (b(i) + b(j)) *(1 - lij(i, j))
            end do
            lambda(1) = lambda(1) + n(i)*b(i)**(1.0_pr/3.0_pr)
            lambda(2) = lambda(2) + n(i)*b(i)**(2.0_pr/3.0_pr)
         end do
         bmix = bmix/sum(n)

         lambda(0) = sum(n)
         lambda(3) = bmix
         eta = bmix/v/xi
         l1l2_l3l0 = lambda(1)*lambda(2)/lambda(3)/lambda(0)
         l23_l0l32 = lambda(2)**3/lambda(0)/lambda(3)**2
         logContribution = log(1._pr - xi*eta)/xi

         Ar_rep = -sum(n)*((1._pr + 3._pr*l1l2_l3l0)*logContribution &
                          + 3._pr/xi*(l23_l0l32 - 1._pr/2._pr - l1l2_l3l0/2._pr) &
                          *(eta + logContribution))
         b_v = bmix/v
         Ar_att = (- sum(n) * log(1.0_pr - b_v) &
             - amix / (R*t*bmix)*1.0_pr / (del1 - del2) &
             * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
             ) * (R * t)

         ar = Ar_rep + Ar_att
      end associate
   end function arfun

   function v0(self, n, P, T)
      class(HardSpheresCubicEoS), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr) :: v0

      v0 = sum(n * self%b)
   end function v0


   subroutine main
      use yaeos
      use forsus, only: Substance, forsus_dir, forsus_default_dir
      use hyperdual_pr76, only: hPr76 => PR76, set_hpr => setup

      integer, parameter :: nc=2

      real(pr) :: n(nc), V, P, Phs, T
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc), lij(nc,nc)
      type(Substance) :: sus(nc)

      type(CubicEoS) :: pr76
      type(hPR76) :: hdpr76
      type(HardSpheresCubicEoS) :: hspr76

      type(EquilibriumState) :: eq
      type(PTEnvel2) :: env

      real(pr) :: Ar, ArT, ArV, ArTV, ArT2, ArV2, Arn(nc), Artn(nc), ArVn(nc), arn2(nc,nc)

      integer :: i

      forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

      sus(1) = Substance("methane")
      sus(2) = Substance("n-decane")

      tc = sus%critical%critical_temperature%value
      pc = sus%critical%critical_pressure%value/1e5
      w = sus%critical%acentric_factor%value

      tc(2) = 874.0
      pc(2) = 6.8
      w(2) = 1.52596

      kij = 0
      lij = 0

      pr76 = PengRobinson76(tc, pc, w, kij, lij)
      hdpr76 = set_hpr(tc, pc, w, kij, lij)

      ! Copy PR76 into HSPR76
      hspr76%components%tc = tc
      hspr76%components%pc = pc
      hspr76%components%w = w
      hspr76%ac = pr76%ac
      hspr76%b = pr76%b

      associate(alpha => pr76%alpha)
         select type(alpha)
          type is (AlphaSoave)
            hspr76%k = alpha%k
         end select
      end associate

      hspr76%kij = kij
      hspr76%lij = lij

      n = [0.3, 0.7]
      V = 1
      T = 400

      block
         real(pr) :: dPdV, k, khs
         do i=1,100
            P = real(i, pr)
            call volume(pr76, n, P, T, V, root_type="stable")
            call pressure(pr76, n, V, T, P, dPdV)
            k = -1/V * 1/dPdV
            
            call volume(hspr76, n, P, T, V, root_type="stable")
            call pressure(hspr76, n, V, T, P, dPdV)
            khs = -1/V * 1/dPdV

            print *, P, k, khs
         end do
      end block

      P = 50
      do i=1,99
         n(1) = real(i)/100
         n(2) = 1 - n(1)

         eq = saturation_pressure(pr76, n, T=T, kind="bubble", P0=P)
         P = eq%p
         if (eq%iters < 1000) write(1, *) eq%x(1), eq%y(1), eq%p

         eq = saturation_pressure(hspr76, n, T=T, kind="bubble", p0=eq%P)
         if (eq%iters < 1000) write(2, *) eq%x(1), eq%y(1), eq%p
      end do
      call exit

      T = 150
      eq = saturation_pressure(pr76, n, T=T, kind="bubble")
      print *, eq%iters, eq
      env = pt_envelope_2ph(pr76, n, eq)
      
      print *, size(env%points)
      write(1, *) env
      
      eq = saturation_pressure(hspr76, n, T=T, kind="bubble", p0=eq%P)
      print *, eq%iters, eq
      
      env = pt_envelope_2ph(hspr76, n, eq)
      print *, size(env%points)
      write(2, *) env
   end subroutine main
end module yaeos__HardSpheresCubicEoS
