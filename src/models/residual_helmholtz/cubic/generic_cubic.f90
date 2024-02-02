module yaeos_models_ar_genericcubic
    use yaeos_constants, only: pr
    use yaeos_models_ar, only: ArModel
    use yaeos_substance, only: Substances
    implicit none
    
    type, abstract :: AlphaFunction
    contains
        procedure(abs_alpha), deferred :: alpha
    end type

    type, abstract :: CubicMixRule
    contains
        procedure(abs_Dmix), deferred :: Dmix
        procedure(abs_Bmix), deferred :: Bmix
    end type

    type, extends(ArModel) :: CubicEoS
        !! Cubic Equation of State.
        class(Substances), allocatable :: components
        class(CubicMixRule), allocatable :: mixrule
        class(AlphaFunction), allocatable :: alpha
        real(pr), allocatable :: ac(:), b(:), del1(:), del2(:)
    contains
        procedure :: residual_helmholtz => GenericCubic_Ar
        procedure :: get_v0 => v0
    end type

    abstract interface
        subroutine abs_alpha(self, Tr, a, dadt, dadt2)
            import AlphaFunction, pr
            class(AlphaFunction), intent(in) :: self
            real(pr), intent(in) :: Tr(:)
            real(pr), intent(out) :: a(:), dadt(:), dadt2(:)
        end subroutine

        subroutine abs_Dmix(self, n, T, &
            ai, daidt, daidt2, &
            D, dDdT, dDdT2, dDi, dDidT, dDij&
            )
            import CubicMixRule, pr
            class(CubicMixRule), intent(in) :: self
            real(pr), intent(in) :: T, n(:)
            real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
            real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)
        end subroutine
        
        subroutine abs_Bmix(self, n, bi, B, dBi, dBij)
            import CubicMixRule, pr
            class(CubicMixRule), intent(in) :: self
            real(pr), intent(in) :: n(:)
            real(pr), intent(in) :: bi(:)
            real(pr), intent(out) :: B, dBi(:), dBij(:, :)
        end subroutine
    end interface
contains

   subroutine GenericCubic_Ar(&
        self, n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
    )
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[P = \frac{RT}{V-b} - \frac{a_c\alpha(T_r)}{(V+b\delta_1)(V+b\delta_2)}\]
      use yaeos_constants, only: R
      class(CubicEoS), intent(in) :: self
      real(pr), intent(in) :: n(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), optional, intent(out) :: ar !! Residual Helmholtz
      real(pr), optional, intent(out) :: arv !! dAr/dV
      real(pr), optional, intent(out) :: ArT !! dAr/dT
      real(pr), optional, intent(out) :: artv !! dAr2/dTV
      real(pr), optional, intent(out) :: arv2 !! dAr2/dV2
      real(pr), optional, intent(out) :: ArT2 !! dAr2/dT2
      real(pr), optional, intent(out) :: Arn(size(n)) !! dAr/dn
      real(pr), optional, intent(out) :: ArVn(size(n)) !! dAr2/dVn
      real(pr), optional, intent(out) :: ArTn(size(n)) !! dAr2/dTn
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! dAr2/dn2


      real(pr) :: Bmix, dBi(size(n)), dBij(size(n), size(n))
      real(pr) :: D, dDi(size(n)), dDij(size(n), size(n)), dDidT(size(n)), dDdT, dDdT2

      real(pr) :: totn, d1, d2

      real(pr) :: f, g, fv, fB, gv, fv2, gv2, AUX, FFB, FFBV, FFBB

      real(pr) :: Tr(size(n)), a(size(n)), dadt(size(n)), dadt2(size(n))


      integer :: i, j, nc

      nc = size(n)
      TOTN = sum(n)

      D1 = self%del1(1)
      D2 = (1._pr - D1)/(1._pr + D1)

      Tr = T/self%components%Tc
      call self%alpha%alpha(Tr, a, dadt, dadt2)
      call self%mixrule%Bmix(n, self%b, Bmix, dBi, dBij)

      
      a = self%ac * a
      dadt = self%ac * dadt / self%components%Tc
      dadt2 = self%ac * dadt2 / self%components%Tc**2

      call self%mixrule%Dmix(&
            n, T, a, dadt, dadt2, D, dDdT, dDdT2, dDi, dDidT, dDij&
      )

      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f and its derivatives as defined by Mollerup
      f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/(D1 - D2)
      g = R*log(1 - Bmix/V)
      fv = -1/((V + D1*Bmix)*(V + D2*Bmix))
      fB = -(f + V*fv)/Bmix
      gv = R*Bmix/(V*(V - Bmix))
      fv2 = (-1/(V + D1*Bmix)**2 + 1/(V + D2*Bmix)**2)/Bmix/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - Bmix)**2)

      ! Reduced Helmholtz Energy and derivatives
      if (present(Ar)) Ar = -TOTN*g*T - D*f
      if (present(ArV)) ArV = -TOTN*gv*T - D*fv
      if (present(ArV2)) ArV2 = -TOTN*gv2*T - D*fv2

      AUX = R*T/(V - Bmix)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
      FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2

      if (present(Arn))  Arn(:)  = -g*T + FFB*dBi(:) - f*dDi(:)
      if (present(ArVn)) ArVn(:) = -gv*T + FFBV*dBi(:) - fv*dDi(:)
      if (present(ArTn)) ArTn(:) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(:) - f*dDidT(:)

      if (present(Arn2)) then
      do i = 1, nc
         do j = 1, i
            Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                        + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
            Arn2(j, i) = Arn2(i, j)
         end do
      end do
      end if

      ! TEMPERATURE DERIVATIVES
      if (present(ArT))  ArT = -TOTN*g - dDdT*f
      if (present(ArTV)) ArTV = -TOTN*gv - dDdT*fV
      if (present(ArT2)) ArT2 = -dDdT2*f
   end subroutine GenericCubic_Ar

   function v0(self, n, p, t)
    class(CubicEoS), intent(in) :: self
    real(pr), intent(in) :: n(:), p, t
    real(pr) :: v0

    real(pr) :: dbi(size(n)), dbij(size(n), size(n))
    call self%mixrule%Bmix(n, self%b, v0, dbi, dbij)
   end function

end module