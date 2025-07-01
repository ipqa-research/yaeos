program main
   use yaeos, only: pr, R
   use yaeos__consistency_gemodel, only: numeric_ge_derivatives
   use hyperdual_mod
   implicit none

   integer, parameter :: nc = 3
   real(pr) :: n(nc), T, Ge, Gen(nc), GenT(nc), dgendn_num(nc), Gen2(nc, nc), dgendn2(nc, nc), GeT, GeT2

   real(pr) :: alpha(nc, nc), b(nc), gij(nc, nc)

   type(hyperdual) :: t_hd, n_hd(nc)


   integer :: i

   gij = reshape([ &
      0.0, 0.1, 0.2, &
      0.3, 0.0, 0.4, &
      0.5, 0.6, 0.0], [nc, nc]) * 100
   alpha = reshape([ &
      0.0, 1.1, 1.2, &
      1.3, 0.0, 1.5, &
      1.6, 1.7, 0.0], [nc, nc])

   b = [0.1, 0.2, 0.3]
   n = [0.1, 0.3, 0.5]
   T = 150.0

   dgendn_num = dgedn(n, T)
   dgendn2 = dgedn2(n, T)
   call excess_gibbs(n, T, Ge, Gen, GeT, GeT2, GenT, Gen2)
   print *, "Excess Gibbs Energy Function:", ge_f(n, T)
   print *, "Excess Gibbs Energy:", Ge
   print *, "Excess Gibbs Energy Derivative:", Gen
   print *, "Numerical Derivative of Excess Gibbs Energy:", dgendn_num

   ! print *, "Anal:"
   ! do i=1,nc
   !    print *, Gen2(i, :)
   ! end do

   ! print *, "Numd:"
   ! do i=1,nc
   !    print *, dgendn2(i, :)
   ! end do

   print *, "GenT:    ", GenT
   print *, "numGenT: ", dgednt(n, T)
   n_hd = n
   t_hd = t
   t_hd%f1 = 1._pr
   t_hd%f2 = 1._pr
   print *, "Hyperdual Ge GeT and GeT2:", ge_hd(n_hd, t_hd)
   print *, "Analyitic Ge GeT and GeT2:", Ge, GeT, GeT2


contains
   subroutine excess_gibbs(n, T, Ge, Gen, GeT, GeT2, GenT, Gen2)
      real(pr), intent(in) :: n(:), T
      real(pr), intent(out) :: Ge, Gen(:), GeT, GeT2, GenT(:), Gen2(:, :)

      real(pr) :: E(nc, nc), U, D
      real(pr) :: dEdT(nc, nc), dEdT2(nc, nc)

      real(pr) :: T1(nc), T2(nc), T3(nc), Dn(nc)

      real(pr) :: xi(nc, nc), theta(nc, nc), omega(nc, nc), eta(nc, nc)
      real(pr) :: xiT(nc), etaT(nc, nc), thetaT(nc), omegaT(nc, nc)
      real(pr) :: xiTT(nc), thetaTT(nc)

      real(pr) :: tau(nc, nc), dtaudt(nc, nc), dtaudt2(nc, nc)
      real(pr) :: denom

      integer :: i, j, k, l

      real(pr) :: dT1dn(nc, nc), dDndn(nc, nc)
      real(pr) :: tau1(nc, nc)
      integer :: m

      call tdep(T, gij, tau, dtaudt, dtaudt2)

      E = exp(-alpha * tau)
      dEdT = - alpha * dtaudt * E
      dEdT2 = -alpha * (dtaudt2 * E + dtaudt * dEdT)

      do i=1,nc
         xi(:, i)     = E(:, i) * b * tau(:, i) * n
         eta(:, i)    = E(:, i) * b * tau(:, i)
         theta(:, i)  = E(:, i) * b * n
         omega(:, i)  = E(:, i) * b

         xiT(i)    = sum(b * n * (dEdT(:, i) * tau(:, i) + E(:, i)*dtaudt(:, i)))
         etaT(:, i)   = b * (dEdT(:, i) * tau(:, i) + E(:, i) * dtaudt(:, i))
         thetaT(i) = sum(dEdT(:, i) * b * n)
         omegaT(:, i) = dEdT(:, i) * b

         xiTT(i) = sum ( b * n *  (&
            (dEdT(:, i) * dtaudt(:, i) + dEdT(:, i) * dtaudt2(:, i)) &
            + (E(:, i) * dtaudt2(:, i) + dEdT(:, i) * dtaudt(:, i)) &
            ) )
         thetaTT(i) = sum(b * n * (dEdT2(:, i)))
      end do

      Ge= 0
      Gen = 0

      do i=1,nc
         U = 0
         D = 0
         do j=1,nc
            U = U + n(j) * b(j) * E(j, i) * tau(j, i)
            D = D + n(j) * b(j) * E(j, i)
         end do
         Ge = Ge + n(i) * U/D
      end do


      do i=1,nc
         Dn(i) = sum(theta(:, i))
         T1(i) = sum(xi(:, i))
      end do

      Gen = 0.0
      do i=1,nc
         Gen(i) = sum(xi(:, i)) / sum(theta(:, i))
         do k=1,nc
            Gen(i) = Gen(i) + n(k) * (&
               eta(i, k)/sum(theta(:, k)) &
               - omega(i, k) * sum(xi(:, k))/sum(theta(:,k))**2 &
               )
         end do
      end do

      Gen2 = 0.0
      do i=1,nc
         do j=1,nc
            Gen2(i, j) = &
               - omega(j, i) * sum(xi(:, i)) / sum(theta(:, i))**2 &
               - omega(i, j) * sum(xi(:, j)) / sum(theta(:, j))**2 &
               + eta(j, i) / sum(theta(:, i)) &
               + eta(i, j) / sum(theta(:, j))

            do k=1,nc
               denom = sum(theta(:, k))
               Gen2(i, j) = Gen2(i, j) + &
                  2 * n(k) * omega(i, k) * omega(j, k) * sum(xi(:, k)) / denom**3 &
                  - n(k) * omega(i, k) * eta(j, k) / denom**2 &
                  - n(k) * omega(j, k) * eta(i, k) / denom**2
            end do
         end do
      end do

      GenT = 0.0
      do i=1,nc
         GenT(i) = xiT(i)/Dn(i) - sum(xi(:, i)) * thetaT(i)/Dn(i)**2
         do k=1,nc
            GenT(i) = GenT(i) + n(k) * (&
               etaT(i, k)/Dn(k) - eta(i, k) * thetaT(k) / Dn(k)**2 &
               - (omega(i, k) * xiT(k) + omegaT(i, k) * sum(xi(:, k))) / Dn(k)**2 &
               + 2 * thetaT(k) * omega(i, k) * sum(xi(:, k)) / Dn(k)**3 &
               )
         end do
      end do

      Ge = Ge * R * T
      Gen = Gen * R * T
      GenT = Gen/T + R * T * GenT
      Gen2 = Gen2 * R * T

      GeT = 0
      do i=1,nc
         GeT = GeT + n(i) * (xiT(i)/sum(theta(:, i)) - sum(xi(:, i))*(thetaT(i))/Dn(i)**2)
      end do

      GeT2 = GeT
      do i=1,nc
         GeT2 = GeT2 + T * n(i) * (&
            xiTT(i)/Dn(i) - xiT(i) * sum(thetaT(:, i))/Dn(i)**2 &
            - xiT(i)*sum(thetaT(:, i))) + sum(xi(:, i)) * thetaTT(i))/Dn(i)**2 &
            + 2 * sum(thetaT(:, i)**2  * sum(xi(:, i))) / Dn(i)**3 &
         )
      end do


      GeT = GeT * R * T + Ge/T
      GeT2 = R * (-Ge/T**2 + GeT/T) + R*GeT2

   end subroutine excess_gibbs

   real(pr) function ge_f(n, T)
      real(pr), intent(in) :: n(:), T

      real(pr) :: E(nc, nc), U, D
      real(pr) :: tau(nc, nc), dtaudt(nc, nc), dtaudt2(nc, nc)

      integer :: i, j

      call tdep(T, gij, tau, dtaudt, dtaudt2)
      E = exp(-alpha * tau)

      ge_f = 0
      do i=1,nc
         U = 0
         D = 0
         do j=1,nc
            U = U + n(j) * b(j) * E(j, i) * tau(j, i)
            D = D + n(j) * b(j) * E(j, i)
         end do
         ge_f = ge_f + n(i) * U/D
      end do

      ge_f = ge_f * R * T
   end function ge_f

   type(hyperdual) function ge_hd(n, T)
      type(hyperdual), intent(in) :: n(:), T

      type(hyperdual) :: E(nc, nc), U, D
      type(hyperdual) :: tau(nc, nc)

      real(pr) :: tin

      integer :: i, j
      tin = t%f0

      ! call tdep(Tin, gij, tau, dtaudt)
      tau = gij/(R*T)
      E = exp(-alpha * tau)

      ge_hd = 0._pr
      do i=1,nc
         U = 0._pr
         D = 0._pr
         do j=1,nc
            U = U + n(j) * b(j) * E(j, i) * tau(j, i)
            D = D + n(j) * b(j) * E(j, i)
         end do
         ge_hd = ge_hd + n(i) * U/D
      end do

      ge_hd = ge_hd * R * T
   end function ge_hd

   function dgedn(n, T)
      real(pr), intent(in) :: n(:), T
      real(pr) :: dgedn(size(n))

      real(pr) :: dx(size(n))

      integer :: i

      do i=1,size(n)
         dx = 0.0
         dx(i) = 1e-4
         dgedn(i) = (ge_f(n + dx, T) - ge_f(n - dx, T)) / (2 * dx(i))
      end do
   end function dgedn

   function dgedn2(n, T)
      real(pr), intent(in) :: n(:), T
      real(pr) :: dgedn2(size(n), size(n))
      real(pr) :: dx(size(n)), eps
      integer :: i, j

      real(pr) :: f, f1, f2, f3, f4, f5, f6

      f = ge_f(n, T)
      eps = 1e-6
      do i=1,nc
         do j=1,nc
            dx = 0
            dx(i) = eps
            dx(j) = eps
            f1 = ge_f(n + dx, T)

            dx(i) = eps
            dx(j) = 0
            f2 = ge_f(n + dx, T)

            dx(i) = 0
            dx(j) = eps
            f3 = ge_f(n + dx, T)

            dx(i) = -eps
            dx(j) = 0
            f4 = ge_f(n + dx, T)

            dx(i) = 0
            dx(j) = -eps
            f5 = ge_f(n + dx, T)

            dx(i) = -eps
            dx(j) = -eps
            f6 = ge_f(n + dx, T)
            dgedn2(i, j) = (f1 - f2 - f3 + 2*f - f4 - f5 + f6) / (2 * eps**2)
         end do

         dx = 0
         dx(i) = 1e-2
         f1 = ge_f(n + dx, T)
         f2 = ge_f(n - dx, T)

         dgendn2(i, i) = (f1 - 2*f + f2) / (dx(i)**2)
      end do

   end function dgedn2

   function dgednt(n, T)
      real(pr), intent(in) :: n(:), T
      real(pr) :: dgednt(size(n))
      real(pr) :: dx(size(n)), eps
      integer :: i, j

      real(pr) :: f, f1, f2, f3, f4, f5, f6, dt

      f = ge_f(n, T)
      eps = 1e-6
      do i=1,nc
         dx = 0
         dx(i) = eps
         dt = eps
         f1 = ge_f(n + dx, T + dt)

         dx(i) = eps
         dt = 0
         f2 = ge_f(n + dx, T + dt)

         dx(i) = 0
         dt = eps
         f3 = ge_f(n + dx, T + dt)

         dx(i) = -eps
         dt = 0
         f4 = ge_f(n + dx, T + dt)

         dx(i) = 0
         dt = -eps
         f5 = ge_f(n + dx, T + dt)

         dx(i) = -eps
         dt = -eps
         f6 = ge_f(n + dx, T + dt)
         dgednt(i) = (f1 - f2 - f3 + 2*f - f4 - f5 + f6) / (2 * eps**2)
      end do

   end function dgednt

   elemental subroutine tdep(T, U, f, df, df2)
      real(pr), intent(in) :: T, U
      real(pr), intent(out) :: f, df, df2
      f = U/(R*T)
      df = -U/(R*T**2)
      df2 = 2*U/(R*T**3)
   end subroutine tdep
end program main


! program main
!    use yaeos, only: pr, R
!    implicit none
!    integer, parameter :: nc=3
!
!    real(pr) :: n(3), RTtau(nc, nc), alpha(nc, nc), b(nc)
!    real(pr) :: Ge_val, Gen(nc)
!    integer :: i
!
!    real(pr) :: T
!
!
!
!    b = [1, 2, 3]
!    n = [0.5_pr, 0.3_pr, 0.2_pr]
!    T = 250
!
!    RTtau = 0
!    RTtau = reshape( [ &
!       0.0_pr, 0.2_pr, 0.3_pr, &
!       0.4_pr, 0.0_pr, 0.6_pr, &
!       0.7_pr, 0.8_pr, 0.0_pr &
!       ], [nc, nc])
!
!    alpha = reshape( [ &
!       0.0_pr, 0.1_pr, 0.2_pr, &
!       0.3_pr, 0.0_pr, 0.4_pr, &
!       0.5_pr, 0.6_pr, 0.0_pr &
!       ], [nc, nc])
!
!    call excess_gibbs(n, T, Ge_val, Gen)
!
!    print *, Ge_val, Ge(n, T)
!
!    do i=1,nc
!       print *, dGe(i), Gen(i)
!    end do
!
! contains
!    real(pr) function Ge(n, T)
!       real(pr), intent(in) :: n(:)
!       real(pr), intent(in) :: T
!
!       real(pr) :: tau(nc, nc)
!       integer :: i, j
!
!       real(pr) :: U, D, E(nc, nc)
!
!       tau = RTtau / (R*T)
!       E = exp(alpha * tau)
!
!       U = 0
!       D = 0
!       Ge = 0
!       do i=1,nc
!          U = 0
!          D = 0
!          do j=1,nc
!             U = U + n(j) * b(j) * tau(j, i) * E(j, i)
!             D = D + n(j) * b(j) * E(j, i)
!          end do
!          Ge = Ge + n(i) * U / D
!       end do
!
!       Ge = Ge * R * T
!
!    end function Ge
!
!    subroutine excess_gibbs(n, T, Ge, Gen)
!       real(pr), intent(in) :: n(:)
!       real(pr), intent(in) :: T
!       real(pr), intent(out) :: Ge
!       real(pr), intent(out) :: Gen(:)
!
!       real(pr) :: tau(nc, nc)
!       integer :: i, j
!
!       real(pr) :: U, D, E(nc, nc)
!
!       real(pr) :: dUDdn(nc), dRdn(nc)
!
!       tau = RTtau / (R*T)
!       E = exp(alpha * tau)
!
!       Ge = 0
!       dRdn = 0
!       do i=1,nc
!          U = 0
!          D = 0
!          do j=1,nc
!             U = U + n(j) * b(j) * tau(j, i) * E(j, i)
!             D = D + n(j) * b(j) * E(j, i)
!
!          end do
!          dRdn(i) = -sum(n * b * tau(:, i)* E(:, i)) / D
!
!          Ge = Ge + n(i) * U / D
!       end do
!
!       Ge = Ge * R * T
!
!    end subroutine excess_gibbs
!
!    real(pr) function dGe(i)
!       real(pr) :: f1, f2, f3
!       real(pr) :: X(nc+1), dX(nc+1)
!
!       integer :: i
!
!       X = [n, T]
!
!       dX = 0
!
!       f1 = Ge(X(:nc), X(nc+1))
!
!       dX(i) = 1e-8
!       X = X + dX
!       f2 = Ge(X(:nc), X(nc+1))
!
!
!       dGe = (f2 - f1) / (dX(i))
!
!    end function dGe
! end program main
!
!
!
! ! program main
! !    use yaeos, only: pr, CubicEoS, EquilibriumState, saturation_pressure, HV, PurePsat, pure_saturation_line
! !    use yaeos__equilibria_boundaries_pure_saturation, only: solve_point
! !    use yaeos__models_ar_cubic_mixing_base, only: lamdba_hv
! !    use yaeos__models_ge_nrtlhv, only: NRTLHV
! !    use yaeos, only: PXEnvel2, px_envelope_2ph
! !    use testing_aux, only: test_title, assert
! !
! !
! !    integer, parameter :: nc=2
! !    type(NRTLHV) :: model
! !    type(CubicEoS) :: cubic
! !    type(PurePsat) :: psat
! !    type(EquilibriumState) :: sat
! !    type(HV) :: mr
! !    type(PXEnvel2) :: px
! !
! !    real(pr) :: C(nc, nc), alpha(nc, nc)
! !    real(pr) :: n(nc), T
! !    real(pr) :: Ge, GeT
! !    real(pr) :: a(nc), kij(nc, nc)
! !
! !    integer :: i
! !
! !    write(*, *) test_title("NRTL-HV model")
! !
! !    call setup_cubic
! !
! !
! !    ! ==========================================================================
! !    ! Valuation of Ge for a mixture of methane and water, using the reduction
! !    ! to classical mixing rules
! !    ! --------------------------------------------------------------------------
! !
! !    n = [0.5, 0.5]
! !    T = 200
! !    call nrtl_params_with_cubic(alpha, C)
! !    model%alpha = alpha
! !    model%b = cubic%b
! !    model%C = C
! !
! !    call model%excess_gibbs(n, T, Ge)
! !
! !    ! call assert(abs(Ge - 46.02789) < 1e-4, "NRTL-HV excess Gibbs energy")
! !
! !    ! ==========================================================================
! !    ! Valuation of Ge for a mixture of acetone and water, using the parameters
! !    ! from the original HV paper
! !    ! --------------------------------------------------------------------------
! !    T = 200 + 273.15
! !    model = acetone_water_paper()
! !    mr = HV(ge=model, del1=cubic%del1, bi=cubic%b)
! !
! !    block
! !       real(pr) :: z0(2), zi(2), a
! !       z0 = [0, 1]
! !       zi = [1, 0]
! !
! !       a = 0.9999
! !       n = a * zi + (1-a)*z0
! !
! !       psat = pure_saturation_line(cubic, 2, 1._pr, T-10)
! !       sat = saturation_pressure(cubic, n, T, kind="bubble", P0=psat%get_P(T))
! !       px = px_envelope_2ph(cubic, z0, a, zi, sat, delta_0=-0.0001_pr)
! !       do i=1,size(px%points)
! !          write(1, *) px%points(i)
! !       end do
! !
! !       call cubic%set_mixrule(mr)
! !       psat = pure_saturation_line(cubic, 2, 1._pr, T-10)
! !       sat = saturation_pressure(cubic, n, T, kind="bubble", P0=psat%get_P(T))
! !       px = px_envelope_2ph(cubic, z0, a, zi, sat, delta_0=-0.0001_pr)
! !       do i=1,size(px%points)
! !          write(2, *) px%points(i)
! !       end do
! !    end block
! !
! !
! ! contains
! !
! !    subroutine setup_cubic
! !       use yaeos, only: SoaveRedlichKwong
! !       use forsus, only: forsus_default_dir, forsus_dir, Substance
! !       type(Substance) :: sus(nc)
! !
! !       forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
! !       kij(1, 2) = 0.1
! !       kij(2, 1) = 0.1
! !
! !       sus(1) = Substance("methane")
! !       sus(2) = Substance("water")
! !
! !       cubic = SoaveRedlichKwong(&
! !          tc=sus%critical%critical_temperature%value, &
! !          pc=sus%critical%critical_pressure%value/1e5, &
! !          w=sus%critical%acentric_factor%value, &
! !          kij=kij)
! !    end subroutine setup_cubic
! !
! !    subroutine nrtl_params_with_cubic(alpha, C)
! !       real(pr), intent(out) :: alpha(nc, nc), C(nc, nc)
! !
! !       real(pr) :: g_ji(nc, nc), g_ii(nc), b(nc), a(nc), dadt(nc), dadt2(nc)
! !       real(pr) :: lambda, lambdadn(nc), dlambdadn2(nc, nc)
! !       real(pr) :: d1, dd1(nc), dd12(nc, nc)
! !
! !       real(pr) :: Tr(nc), Tc(nc)
! !       integer :: i, j
! !       g_ji = 0
! !       alpha = 0
! !
! !       ! ==============================================================
! !       ! Obtain parameters
! !       ! --------------------------------------------------------------
! !       b = cubic%b
! !       Tc = cubic%components%Tc
! !       Tr = T / Tc
! !       call cubic%alpha%alpha(Tr, a, dadt, dadt2)
! !       a = cubic%ac * a
! !
! !       call cubic%mixrule%D1mix(n, cubic%del1, d1, dd1, dd12)
! !       call lamdba_hv(d1, dd1, dd12, lambda, lambdadn, dlambdadn2)
! !
! !       g_ii = -a/b * lambda
! !       do i=1,nc
! !          do j=1,nc
! !             g_ji(i, j) = -2 * sqrt(b(i)*b(j))/(b(i) + b(j)) * sqrt(g_ii(i)*g_ii(j)) * (1-kij(i, j))
! !          end do
! !       end do
! !
! !       do i=1,nc
! !          C(:, i) = g_ji(:, i) - g_ii(i)
! !       end do
! !
! !       alpha = 0
! !    end subroutine nrtl_params_with_cubic
! !
! !    subroutine C_from_kij(cubic, T, C)
! !       use yaeos, only: size, QMR, R
! !       type(CubicEoS) :: cubic
! !       real(pr), intent(in) :: T
! !       real(pr), intent(out) :: C(:, :)
! !
! !       real(pr) :: a(size(cubic)), b(size(cubic)), kij(size(cubic), size(cubic))
! !       real(pr) :: dadt(size(cubic)), dadt2(size(cubic))
! !       real(pr) :: lambda, lambdadn(size(cubic)), dlambdadn2(size(cubic), size(cubic))
! !
! !       real(pr) :: d1, dd1(size(cubic)), dd12(size(cubic), size(cubic))
! !       real(pr) :: g_ji(size(cubic), size(cubic)), g_ii(size(cubic))
! !
! !       real(pr) :: Tr(size(cubic)), Tc(size(cubic))
! !       integer :: i, j, k
! !
! !       real(pr) :: D, dDT, dDT2, dDi(size(cubic)), dDiT(size(cubic)), dDij(size(cubic), size(cubic))
! !       real(pr) :: aux
! !       integer :: f
! !
! !       open(newunit=f, file="D")
! !
! !       b = cubic%b
! !       Tc = cubic%components%Tc
! !       Tr = T / Tc
! !
! !       ! ==============================================================
! !       ! Obtain the a_i values
! !       ! --------------------------------------------------------------
! !       call cubic%alpha%alpha(Tr, a, dadt, dadt2)
! !       a = cubic%ac * a
! !       dadt = cubic%ac * dadt / Tc
! !       dadt2 = cubic%ac * dadt2 / Tc**2
! !
! !       associate(mr => cubic%mixrule)
! !          select type(mr)
! !           class is (QMR)
! !             kij = mr%k
! !          end select
! !       end associate
! !
! !       call cubic%mixrule%D1mix(n, cubic%del1, d1, dd1, dd12)
! !       call lamdba_hv(d1, dd1, dd12, lambda, lambdadn, dlambdadn2)
! !       print *, "D1", d1
! !
! !       do i=1,99
! !          n = [0.01*i, 1-0.01*i]
! !          call cubic%mixrule%Dmix(n, T, a, dadt, dadt2, D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDiT, dDij=dDij)
! !          write(f, *) n(1), D, dDT, dDT2, dDi(1)
! !       end do
! !       write(f, *)
! !       write(f, *)
! !
! !       g_ii = -a/b * lambda
! !       do i=1,size(cubic)
! !          do j=1,size(cubic)
! !             g_ji(j, i) = -2*lambda * sqrt(a(i)*a(j)) * 1/(b(i)+b(j)) * (1-kij(i, j))
! !          end do
! !       end do
! !
! !       do i=1,size(cubic)
! !          do j=1,size(cubic)
! !             C(j, i) = (g_ji(j, i) - g_ji(i,i))
! !          end do
! !       end do
! !
! !       alpha = 0.0
! !       model%alpha = alpha
! !       model%C = C
! !       model%b = cubic%b
! !
! !       mr = HV(ge=model, del1=cubic%del1, bi=cubic%b)
! !       call cubic%set_mixrule(mr)
! !
! !       call cubic%mixrule%Dmix(&
! !          n, T, a, dadt, dadt2, &
! !          D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDiT, dDij=dDij)
! !
! !
! !       print *, nc
! !       do i=1,99
! !          n = [0.01*i, 1-0.01*i]
! !
! !          call model%excess_gibbs(n, T, Ge)
! !          aux = 0
! !          do j=1,nc
! !             do k=1,nc
! !                aux = aux + n(j)*n(k)*sqrt(a(j)*a(k)) * b(k)/(b(j) + b(k)) * (1-kij(j, k))
! !             end do
! !          end do
! !
! !          aux = -2*lambda/(sum(n*b)) * aux
! !          aux = aux + lambda * sum(n*a/b)
! !
! !          call cubic%mixrule%Dmix(n, T, a, dadt, dadt2, D=D, dDdT=dDT, dDdT2=dDT2, dDi=dDi, dDidT=dDiT, dDij=dDij)
! !          write(f, *) n(1), D, dDT, dDT2, dDi(1)
! !       end do
! !
! !       call exit
! !
! !
! !    end subroutine C_from_kij
! !
! !    type(NRTLHV) function acetone_water_paper() result(model)
! !       !! Parameters from the original HV paper for the system acetone-water
! !       use yaeos, only: SoaveRedlichKwong
! !       use forsus, only: forsus_default_dir, forsus_dir, Substance
! !
! !       type(Substance) :: sus(2)
! !       real(pr) :: alpha(2, 2), C(2, 2), b(2), kij(2,2)
! !
! !       forsus_dir = "build/dependencies/forsus/" // forsus_default_dir
! !       sus(1) = Substance("acetone")
! !       sus(2) = Substance("water")
! !
! !       kij = 0
! !       kij(1, 2) = -0.087
! !       kij(2, 1) = kij(1, 2)
! !       cubic = SoaveRedlichKwong(&
! !          tc=sus%critical%critical_temperature%value, &
! !          pc=sus%critical%critical_pressure%value/1e5, &
! !          w=sus%critical%acentric_factor%value, kij=kij)
! !
! !       alpha = 0
! !       alpha(1, 2) = 0.461
! !       alpha(2, 1) = 0.461
! !
! !       C = 0
! !       C(1, 2) = 73011._pr * 100000._pr/101325._pr / 1000
! !       C(2, 1) = 58849._pr * 100000._pr/101325._pr / 1000
! !
! !       ! model%alpha = alpha
! !       ! model%C = C
! !       ! model%b = cubic%b
! !
! !       call C_from_kij(cubic, T, C)
! !    end function acetone_water_paper
! ! end program main
! !
!
