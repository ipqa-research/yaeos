module yaeos__models_ar_genericcubic_base
   use yaeos__constants, only: pr, R
   implicit none

contains

   subroutine GenericCubic_Ar(&
      n, V, T, &
      B, dBi, dBij, &
      D, dDi, dDij, dDidT, dDdT, dDdT2, &
      D1, dD1i, dD1ij, &
      Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
      )
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[
      !!   P = \frac{RT}{V-B}
      !!       - \frac{D(T)}{(V+B \delta_1)(V+B\delta_2)}
      !! \]
      !!
      !! This routine assumes that the \(\delta_1\) is not a constant parameter
      !! (as it uses to be in classical Cubic EoS) to be compatible with the
      !! three parameter EoS RKPR where \(delta_1\) is not a constant and
      !! has its own mixing rule.
      use yaeos__constants, only: R
      real(pr), intent(in) :: n(:) !! Number of moles
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]

      real(pr), intent(in) :: B !! Repulsive parameter [L]
      real(pr), intent(in) :: dBi(size(n)) !! \(dB/dn_i\)
      real(pr), intent(in) :: dBij(size(n), size(n)) !! \(d^2B/dn_{ij}\)

      real(pr), intent(in) :: D !! Attractive parameter
      real(pr), intent(in) :: dDi(size(n)) !! \(dD/dn_i\)
      real(pr), intent(in) :: dDij(size(n), size(n)) !! \(d^2D/dn_{ij}\)
      real(pr), intent(in) :: dDidT(size(n)) !! \(d^2D/dTdn_i\)
      real(pr), intent(in) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(in) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)

      real(pr), intent(in) :: D1 !! \(\delta_1\) parameter
      real(pr), intent(in) :: dD1i(size(n)) !! \(d\delta_1/dn_i\)
      real(pr), intent(in) :: dD1ij(size(n), size(n)) !! \(d^2\delta_1/dn_{ij}\)

      real(pr), optional, intent(out) :: ar !! Residual Helmholtz
      real(pr), optional, intent(out) :: arv !! \(\frac{dAr}{dV}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
      real(pr), optional, intent(out) :: artv !! \(\frac{d^2Ar}{dTdV}\)
      real(pr), optional, intent(out) :: arv2 !! \(\frac{d^2Ar}{dV^2}\)
      real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
      real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVdn_i}\)
      real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTdn_i}\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\frac{d^2Ar}{dn_{ij}}\)


      real(pr) :: totn
      real(pr) :: auxD2, fD1, fBD1, fVD1, fD1D1
      real(pr) d2

      real(pr) :: f, g, fv, fB, gv, fv2, gv2, AUX, FFB, FFBV, FFBB

      integer :: i, j, nc

      nc = size(n)
      TOTN = sum(n)

      ! Delta 2 parameter is calculated from Delta 1
      D2 = (1._pr - D1)/(1._pr + D1)

      ! ========================================================================
      ! Main functions defined by Møllerup
      ! The f's and g's used here are for Ar, not F (reduced Ar)
      ! This requires to multiply by R all g, f
      ! ------------------------------------------------------------------------
      f = log((V + D1*B)/(V + D2*B))/B/(D1 - D2)
      g = R*log(1 - B/V)
      fv = -1/((V + D1*B)*(V + D2*B))
      fB = -(f + V*fv)/B
      gv = R*B/(V*(V - B))
      fv2 = (-1/(V + D1*B)**2 + 1/(V + D2*B)**2)/B/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - B)**2)

      ! DERIVATIVES OF f WITH RESPECT TO DELTA1
      auxD2 = (1 + 2/(1 + D1)**2)
      fD1 = (1/(V + D1*B) + 2/(V + D2*B)/(1 + D1)**2) - f*auxD2
      fD1 = fD1/(D1 - D2)
      fBD1 = -(fB*auxD2 + D1/(V + D1*B)**2 + 2*D2/(V + D2*B)**2/(1 + D1)**2)
      fBD1 = fBD1/(D1 - D2)
      fVD1 = -(fV*auxD2 + 1/(V + D1*B)**2 + 2/(V + D2*B)**2/(1 + D1)**2)/(D1 - D2)
      fD1D1 = 4*(f - 1/(V + D2*B))/(1 + D1)**3 + B*(-1/(V + D1*B)**2 &
         + 4/(V + D2*B)**2/(1 + D1)**4) - 2*fD1*(1 + 2/(1 + D1)**2)
      fD1D1 = fD1D1/(D1 - D2)

      AUX = R*T/(V - B)
      FFB = TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - B) + D*(2*fv + V*fv2)/B
      FFBB = TOTN*AUX/(V - B) - D*(2*f + 4*V*fv + V**2*fv2)/B**2

      ! ========================================================================
      ! Reduced Helmholtz Energy and derivatives
      ! ------------------------------------------------------------------------
      if (present(Ar)) Ar = -TOTN*g*T - D*f
      if (present(ArV)) ArV = -TOTN*gv*T - D*fv
      if (present(ArV2)) ArV2 = -TOTN*gv2*T - D*fv2

      if (present(Arn))  Arn(:)  = -g*T + FFB*dBi(:) - f*dDi(:) - D*fD1 * dD1i(:)
      if (present(ArVn)) ArVn(:) = -gv*T + FFBV*dBi(:) - fv*dDi(:) - D*fVD1*dD1i(:)
      if (present(ArTn)) ArTn(:) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(:) - f*dDidT(:) - dDdT*fD1*dD1i(:)

      if (present(Arn2)) then
         do i = 1, nc
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                  + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) - f*dDij(i, j)
               Arn2(i, j) = Arn2(i, j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                  - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                  - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end do
      end if

      ! TEMPERATURE DERIVATIVES
      if (present(ArT))  ArT = -TOTN*g - dDdT*f
      if (present(ArTV)) ArTV = -TOTN*gv - dDdT*fV
      if (present(ArT2)) ArT2 = -dDdT2*f
   end subroutine GenericCubic_Ar

   subroutine GenericCubic_Ar_dddlc(&
      n, V, T, &
      B, dBi, dBij, &
      D, &
      dDdV, dDdV2, dDdT, dDdT2, &
      dDdTV, dDi, dDidV, dDidT, dDij, &
      D1, dD1i, dD1ij, &
      Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz Energy for a generic Cubic Equation of State.
      !!
      !! Calculates the residual Helmholtz Energy for a generic Cubic EoS as
      !! defined by Michelsen and Møllerup:
      !!
      !! \[
      !!   P = \frac{RT}{V-B}
      !!       - \frac{D(T,V,n)}{(V+B \delta_1)(V+B\delta_2)}
      !! \]
      !!
      !! When the mixing rule produces a D that depends on volume (e.g. the
      !! s-DDLC local composition rule), the non-zero volume derivatives
      !! `dDdV`, `dDdV2`, `dDdTV`, and `dDidV` drive additional correction
      !! terms that are added to `ArV`, `ArV2`, `ArTV`, and `ArVn`.
      !! For standard mixing rules all four inputs should be passed as zero,
      !! and the correction block is skipped entirely at run time.
      use yaeos__constants, only: R
      real(pr), intent(in) :: n(:)  !! Mole numbers
      real(pr), intent(in) :: V     !! Volume [L]
      real(pr), intent(in) :: T     !! Temperature [K]

      real(pr), intent(in) :: B             !! Repulsive parameter [L]
      real(pr), intent(in) :: dBi(size(n)) !! \(dB/dn_i\)
      real(pr), intent(in) :: dBij(size(n), size(n)) !! \(d^2B/dn_{ij}\)

      real(pr), intent(in) :: D                       !! Attractive parameter
      real(pr), intent(in) :: dDi(size(n))            !! \(dD/dn_i\)
      real(pr), intent(in) :: dDij(size(n), size(n))  !! \(d^2D/dn_{ij}\)
      real(pr), intent(in) :: dDidT(size(n))          !! \(d^2D/dT\,dn_i\)
      real(pr), intent(in) :: dDidV(size(n))          !! \(d^2D/dV\,dn_i\)
      real(pr), intent(in) :: dDdT  !! \(dD/dT\)
      real(pr), intent(in) :: dDdT2 !! \(d^2D/dT^2\)
      real(pr), intent(in) :: dDdV  !! \(dD/dV\) — zero for V-independent mixing rules
      real(pr), intent(in) :: dDdV2 !! \(d^2D/dV^2\) — zero for V-independent mixing rules
      real(pr), intent(in) :: dDdTV !! \(d^2D/dT\,dV\) — zero for V-independent mixing rules

      real(pr), intent(in) :: D1              !! \(\delta_1\) parameter
      real(pr), intent(in) :: dD1i(size(n))   !! \(d\delta_1/dn_i\)
      real(pr), intent(in) :: dD1ij(size(n), size(n)) !! \(d^2\delta_1/dn_{ij}\)

      real(pr), optional, intent(out) :: Ar              !! Residual Helmholtz
      real(pr), optional, intent(out) :: ArV             !! \(dAr/dV\)
      real(pr), optional, intent(out) :: ArT             !! \(dAr/dT\)
      real(pr), optional, intent(out) :: ArTV            !! \(d^2Ar/dT\,dV\)
      real(pr), optional, intent(out) :: ArV2            !! \(d^2Ar/dV^2\)
      real(pr), optional, intent(out) :: ArT2            !! \(d^2Ar/dT^2\)
      real(pr), optional, intent(out) :: Arn(size(n))    !! \(dAr/dn_i\)
      real(pr), optional, intent(out) :: ArVn(size(n))   !! \(d^2Ar/dV\,dn_i\)
      real(pr), optional, intent(out) :: ArTn(size(n))   !! \(d^2Ar/dT\,dn_i\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(d^2Ar/dn_{ij}\)

      real(pr) :: totn
      real(pr) :: auxD2, fD1, fBD1, fVD1, fD1D1
      real(pr) :: D2
      real(pr) :: f, g, fv, fB, gv, fv2, gv2, AUX, FFB, FFBV, FFBB

      integer :: i, j, nc

      nc   = size(n)
      TOTN = sum(n)

      ! Delta 2 from Delta 1
      ! -----------------------------------------------------------------------
      D2 = (1._pr - D1)/(1._pr + D1)

      ! =======================================================================
      ! Møllerup auxiliary functions  (for Ar, not reduced F, so all g/f × R)
      ! -----------------------------------------------------------------------
      f   = log((V + D1*B)/(V + D2*B))/B/(D1 - D2)
      g   = R*log(1 - B/V)
      fv  = -1/((V + D1*B)*(V + D2*B))
      fB  = -(f + V*fv)/B
      gv  = R*B/(V*(V - B))
      fv2 = (-1/(V + D1*B)**2 + 1/(V + D2*B)**2)/B/(D1 - D2)
      gv2 = R*(1/V**2 - 1/(V - B)**2)

      ! Derivatives of f with respect to delta1
      auxD2 = (1 + 2/(1 + D1)**2)
      fD1   = (1/(V + D1*B) + 2/(V + D2*B)/(1 + D1)**2) - f*auxD2
      fD1   = fD1/(D1 - D2)
      fBD1  = -(fB*auxD2 + D1/(V + D1*B)**2 + 2*D2/(V + D2*B)**2/(1 + D1)**2)
      fBD1  = fBD1/(D1 - D2)
      fVD1  = -(fV*auxD2 + 1/(V + D1*B)**2 + 2/(V + D2*B)**2/(1 + D1)**2)/(D1 - D2)
      fD1D1 = 4*(f - 1/(V + D2*B))/(1 + D1)**3 &
         + B*(-1/(V + D1*B)**2 + 4/(V + D2*B)**2/(1 + D1)**4) &
         - 2*fD1*(1 + 2/(1 + D1)**2)
      fD1D1 = fD1D1/(D1 - D2)

      AUX  = R*T/(V - B)
      FFB  =  TOTN*AUX - D*fB
      FFBV = -TOTN*AUX/(V - B) + D*(2*fv + V*fv2)/B
      FFBB =  TOTN*AUX/(V - B) - D*(2*f + 4*V*fv + V**2*fv2)/B**2

      ! ========================================================================
      ! Reduced Helmholtz Energy and derivatives
      ! ------------------------------------------------------------------------
      if (present(Ar))   Ar   = -TOTN*g*T - D*f

      if (present(ArV)) then
         ArV = -TOTN*gv*T - D*fv
         ! s-DDLC correction: ArV += (-f)*dDdV
         ArV = ArV - f*dDdV
      end if

      if (present(ArV2)) then
         ArV2 = -TOTN*gv2*T - D*fv2
         ! s-DDLC correction: ArV2 += 2*(-fv)*dDdV + (-f)*dDdV2
         ArV2 = ArV2 - 2*fv*dDdV - f*dDdV2
      end if

      ! =======================================================================
      ! Composition derivatives
      ! -----------------------------------------------------------------------
      if (present(Arn)) &
         Arn(:) = -g*T + FFB*dBi(:) - f*dDi(:) - D*fD1*dD1i(:)

      if (present(ArVn)) then
         ArVn(:) = -gv*T + FFBV*dBi(:) - fv*dDi(:) - D*fVD1*dD1i(:)
         ArVn(:) = ArVn(:) - fB*dBi(:)*dDdV - f*dDidV(:)
      end if

      if (present(Arn2)) then
         do i = 1, nc
            do j = 1, i
               Arn2(i, j) = AUX*(dBi(i) + dBi(j)) &
                  - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i)) &
                  + FFB*dBij(i, j) + FFBB*dBi(i)*dBi(j) &
                  - f*dDij(i, j) &
                  - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i)) &
                  - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) &
                  - D*fD1*dD1ij(i, j) - D*fD1D1*dD1i(i)*dD1i(j)
               Arn2(j, i) = Arn2(i, j)
            end do
         end do
      end if

      ! =======================================================================
      ! Temperature derivatives
      ! -----------------------------------------------------------------------
      if (present(ArT))  ArT  = -TOTN*g - dDdT*f
      if (present(ArT2)) ArT2 = -dDdT2*f

      if (present(ArTV)) then
         ArTV = -TOTN*gv - dDdT*fV
         ! s-DDLC correction: ArTV += dDdV*f/T + (-f)*dDdTV
         ArTV = ArTV + dDdV*f/T - f*dDdTV
      end if

      if (present(ArTn)) &
         ArTn(:) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(:) - f*dDidT(:) - dDdT*fD1*dD1i(:)
   end subroutine GenericCubic_Ar_dddlc

end module yaeos__models_ar_genericcubic_base
