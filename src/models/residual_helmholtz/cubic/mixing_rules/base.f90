module yaeos__models_ar_cubic_mixing_base
   !! # Mixing rules core math
   !! Procedures of the core calculations of CubicEoS mixing rules.
   !!
   !! # Description
   !! This module holds all the basic math to use mixing rules in other codes.
   !! Keeping it simple and accessible.
   !!
   !! # Examples
   !!
   !! ```fortran
   !! bi = [0.2_pr, 0.3_pr]
   !! lij = reshape([0.0_pr, 0.2_pr, 0.2_pr, 0.0_pr], [2,2])
   !!
   !! ! Calculate B parameter with Quadratic Mixing Rules.
   !! call bmix_qmr(n, bi, lij, b, dbi, dbij)
   !!
   !! ```
   !!
   !! # References
   use yaeos__constants, only: pr, solving_volume
   implicit none
contains

   pure subroutine bmix_linear(n, bi, b, dbi, dbij)
      !! # bmix_linear
      !! Calculate the linear mixture covolume and composition derivatives.
      !!
      !! # Description
      !! Uses \(B=\sum_i n_i b_i\), so its mole-number gradient is `bi` and
      !! its Hessian is zero.
      !!
      !! # Examples
      !! ```fortran
      !! call bmix_linear(n, bi, B, dBdn, d2Bdn2)
      !! ```
      real(pr), intent(in) :: n(:) !! Component mole amounts.
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes.
      real(pr), intent(out) :: b !! Mixture covolume.
      real(pr), intent(out) :: dbi(:) !! Mole-number gradient of the mixture covolume.
      real(pr), intent(out) :: dbij(:, :) !! Mole-number Hessian of the mixture covolume.

      b = sum(n*bi)
      dbi = bi
      dbij = 0
   end subroutine bmix_linear

   pure subroutine bmix_qmr(n, bi, lij, b, dbi, dbij)
      !! # bmix_qmr
      !! Calculate the QMR covolume and its first and second mole derivatives.
      !!
      !! # Description
      !! Forms pair covolumes \(b_{ij}=\frac{b_i+b_j}{2}(1-l_{ij})\) and
      !! evaluates \(B=\sum_{ij}n_i n_j b_{ij}/\sum_i n_i\). The derivatives
      !! are returned with respect to the component mole amounts.
      !!
      !! # Examples
      !! ```fortran
      !! call bmix_qmr(n, bi, lij, B, dBdn, d2Bdn2)
      !! ```
      real(pr), intent(in) :: n(:) !! Component mole amounts.
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes.
      real(pr), intent(in) :: lij(:, :) !! Repulsive binary interaction matrix.
      real(pr), intent(out) :: b !! Mixture covolume.
      real(pr), intent(out) :: dbi(:) !! Mole-number gradient of the mixture covolume.
      real(pr), intent(out) :: dbij(:, :) !! Mole-number Hessian of the mixture covolume.

      real(pr) :: bij(size(n), size(n))

      real(pr) :: totn !! Total number of moles.
      real(pr) :: aux(size(n)) !! Intermediate accumulator.

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.

      nc = size(n)
      totn = sum(n)
      B = 0
      dBi = 0
      dBij = 0
      aux = 0

      do i = 1, nc
         do j = 1, nc
            bij(i, j) = 0.5_pr * (bi(i) + bi(j)) * (1.0_pr - lij(i,j))
            aux(i) = aux(i) + n(j) * bij(i, j)
         end do
         B = B + n(i)*aux(i)
      end do

      B = B/totn

      if (solving_volume) return

      do i = 1, nc
         dBi(i) = (2*aux(i) - B)/totn
         do j = 1, i
            dBij(i, j) = (2*bij(i, j) - dBi(i) - dBi(j))/totn
            dBij(j, i) = dBij(i, j)
         end do
      end do
   end subroutine bmix_qmr

   pure subroutine d1mix_rkpr(n, d1i, d1, dd1i, dd1ij)
      !! # d1mix_rkpr
      !!
      !! # Description
      !! RKPR \(\delta_1\) parameter mixing rule.
      !!
      !! The RKPR EoS doesn't have a constant \(\delta_1\) value for each
      !! component, so a proper mixing rule should be provided. A linear
      !! combination is used.
      !!
      !! \[
      !!     \Delta_1 = \frac{\sum_i^N n_i \delta_{1i}}{\sum_i^N n_i}
      !! \]
      !!
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: d1i(:) !! Pure-component \(\delta_{1i}\) values
      real(pr), intent(out) :: D1 !! Mixture \(\Delta_1\)
      real(pr), intent(out) :: dD1i(:) !! \(\partial \Delta_1/\partial n_i\)
      real(pr), intent(out) :: dD1ij(:, :) !! \(\partial^2 \Delta_1/\partial n_i \partial n_j\)

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.
      real(pr) :: totn !! Total number of moles.

      nc = size(n)
      totn = sum(n)

      D1 = sum(n * d1i)/totn

      if (solving_volume) return

      do i = 1, nc
         dD1i(i) = (d1i(i) - D1)/totn
         do j = 1, nc
            dD1ij(i, j) = (2 * D1 - d1i(i) - d1i(j))/totn**2
         end do
      end do
   end subroutine d1mix_rkpr

   subroutine lamdba_hv(nc, d1, dd1i, dd1ij, L, dLi, dLij)
      !! # lamdba_hv
      !!
      !! # Description
      !! Infinite pressure limit parameter \(\Lambda\)
      !!
      !! \[
      !! \Lambda = \frac{1}{\delta_1 - \delta_2} \ln \frac{1 + \delta_1}{1 + \delta_2}
      !! \]
      !!
      !! with \(\delta_2 = (1-\delta_1)/(1+\delta_1)\). The mole derivatives
      !! are only computed when `dLij` is present and no volume solve is active.
      integer, intent(in) :: nc !! Number of components
      real(pr), intent(in) :: d1 !! Mixture \(\delta_1\)
      real(pr), optional, intent(in) :: dd1i(nc) !! \(\partial \delta_1/\partial n_i\)
      real(pr), optional, intent(in) :: dd1ij(nc, nc) !! \(\partial^2 \delta_1/\partial n_i \partial n_j\)
      real(pr), intent(out) :: L !! \(\Lambda\)
      real(pr), optional, intent(out) :: dLi(nc) !! \(\partial \Lambda/\partial n_i\)
      real(pr), optional, intent(out) :: dLij(nc, nc) !! \(\partial^2 \Lambda/\partial n_i \partial n_j\)

      real(pr) :: f !! Residual or scalar equation value.
      real(pr) :: g
      real(pr) :: h
      real(pr), dimension(nc) :: df
      real(pr), dimension(nc) :: dg
      real(pr), dimension(nc) :: dh
      real(pr), dimension(nc, nc) :: d2f
      real(pr), dimension(nc, nc) :: d2g
      real(pr), dimension(nc, nc) :: d2h

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      f = d1 + 1
      g = (d1 + 1)*d1 + d1 - 1
      h = log((d1+1)**2 / 2)

      L = f/g * h

      if (solving_volume .or. .not. present(dLij)) return

      df = dd1i
      dg = 2*(d1 + 1)*dd1i
      dh = 2 * dd1i/(d1 + 1)

      dLi = f/g * dh - f*h*dg/g**2 + h * df/g

      do concurrent (i=1:nc, j=1:nc)
         d2f(i, j) = dd1ij(i, j)
         d2g(i, j) = 2*dd1ij(i, j)*(d1 + 1) + 2*dd1i(i)*dd1i(j)
         d2h(i, j) = 2*(dd1ij(i, j)/(d1 + 1) - dd1i(i)*dd1i(j)/(d1 + 1)**2)
      end do

      ! This derivative probably could be simplifyied
      do concurrent (i=1:nc, j=1:nc)
         dLij(i, j) = &
            f * d2h(i,j)/g - &
            f * h *d2g(i, j)/g**2 - &
            f * dg(i) * dh(j)/g**2 - &
            f * dg(j) * dh(i)/g**2 + &
            2 * f * h * dg(i) * dg(j)/g**3 + &
            h * d2f(i, j)/g + &
            df(i)*dh(j)/g + &
            df(j)*dh(i)/g - &
            h * df(i)*dg(j)/g**2 - &
            h * df(j) * dg(i)/g**2
      end do
   end subroutine lamdba_hv

   subroutine DmixHV(n, T, &
      bi, B, dBi, dBij, &
      D1, dD1i, dD1ij, &
      ai, daidt, daidt2, &
      Ge, GeT, GeT2, Gen, GeTn, Gen2, &
      D, dDdT, dDdT2, dDi, dDidT, dDij &
      )
      !! # `DmixHV`
      !! Attractive parameter calculation for the Huron-Vidal mixing rule.
      !!
      !! # Description
      !! This subroutine calculates the attractive parameter \(D\) and its
      !! derivatives for a mixture using the Huron-Vidal mixing rule.
      !! The Huron-Vidal mixing rule combines the pure component parameters
      !! using an excess Gibbs energy model.
      !! The expression of the attractive parameter is:
      !!
      !! \[
      !!   D(n, T) =
      !!     B\left(\sum_i n_i\frac{a_i}{b_i}
      !!     - \frac{G^E}{\Lambda}\right)
      !! \]
      !!
      !! where \(\Lambda\) is evaluated by `lamdba_hv` from the mixture
      !! \(\delta_1\). The excess Gibbs energy and its derivatives are supplied
      !! by the caller, normally from a `GeModel` evaluated at \(T\).
      !!
      !! # References
      !! Huron, M.-J., and J. Vidal, *Fluid Phase Equilibria* 3 (1979) 255-271.
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: bi(:) !! Pure-component covolumes [L/mol]
      real(pr), intent(in) :: B !! Mixture covolume \(B\) [L]
      real(pr), intent(in) :: dBi(:) !! \(\partial B/\partial n_i\)
      real(pr), intent(in) :: dBij(:, :) !! \(\partial^2 B/\partial n_i \partial n_j\)
      real(pr), intent(in) :: D1 !! Mixture \(\delta_1\)
      real(pr), intent(in) :: dD1i(:) !! \(\partial \delta_1/\partial n_i\)
      real(pr), intent(in) :: dD1ij(:, :) !! \(\partial^2 \delta_1/\partial n_i \partial n_j\)
      real(pr), intent(in) :: ai(:) !! Pure-component attractive parameters \(a_i(T)\)
      real(pr), intent(in) :: daidt(:) !! \(da_i/dT\)
      real(pr), intent(in) :: daidt2(:) !! \(d^2a_i/dT^2\)
      real(pr), intent(in) :: Ge !! Excess Gibbs energy \(G^E\) [bar L]
      real(pr), intent(in) :: GeT !! \(\partial G^E/\partial T\)
      real(pr), intent(in) :: GeT2 !! \(\partial^2 G^E/\partial T^2\)
      real(pr), intent(in) :: Gen(:) !! \(\partial G^E/\partial n_i\)
      real(pr), intent(in) :: GeTn(:) !! \(\partial^2 G^E/\partial T \partial n_i\)
      real(pr), intent(in) :: Gen2(:, :) !! \(\partial^2 G^E/\partial n_i \partial n_j\)
      real(pr), intent(out) :: D !! Mixture attractive parameter \(D\) [bar L^2]
      real(pr), intent(out) :: dDdT !! \(\partial D/\partial T\)
      real(pr), intent(out) :: dDdT2 !! \(\partial^2 D/\partial T^2\)
      real(pr), intent(out) :: dDi(:) !! \(\partial D/\partial n_i\)
      real(pr), intent(out) :: dDidT(:) !! \(\partial^2 D/\partial T \partial n_i\)
      real(pr), intent(out) :: dDij(:, :) !! \(\partial^2 D/\partial n_i \partial n_j\)

      real(pr) :: f !! Residual or scalar equation value.
      real(pr) :: fdt
      real(pr) :: fdt2 !! Second-order intermediate quantity used by this procedure.
      real(pr) :: fdi(size(n))
      real(pr) :: fdit(size(n))
      real(pr) :: fdij(size(n), size(n))
      real(pr) :: totn !! Total number of moles

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: nc !! Number of components.
      real(pr) :: L !! Infinite-pressure parameter \(\Lambda\)
      real(pr) :: dL(size(n))
      real(pr) :: dL2(size(n), size(n))

      nc = size(n)
      totn = sum(n)

      call lamdba_hv(nc, D1, dD1i, dD1ij, L, dL, dL2)

      f    = sum(n*ai/bi) - Ge/L
      fdt  = sum(n*daidt/bi) - GeT/L
      fdt2 = sum(n*daidt2/bi) - GeT2/L

      fdi = ai/bi - (Gen/L - dL * Ge/L**2)
      fdiT = daidt/bi - (GeTn/L - dL * GeT/L**2)

      do concurrent(i=1:nc, j=1:nc)
         fdij(i, j) = &
            Ge * dL2(i, j) / L**2 &
            - 2 * Ge * dL(i) * dL(j) / L**3 &
            - Gen2(i, j) / L &
            + Gen(i) * dL(j) / L**2 &
            + Gen(j) * dL(i) / L**2
      end do

      dDi = B*fdi + f*dBi
      dDidT = B*fdiT + fdT*dBi

      D = f*B
      dDdT = fdT*B
      dDdT2 = fdT2*B
      dDij = fdij

      do i=1,nc
         do j=1,nc
            dDij(i, j) = dBi(j)*fdi(i) + B*fdij(j, i) + f*dBij(i, j) + fdi(j)*dBi(i)
         end do
      end do

   end subroutine DmixHV

   ! ===========================================================================
   ! Cubic Mixing Rules
   ! ---------------------------------------------------------------------------

   pure subroutine CMR_aijk(&
      ai, daidt, daidt2, k, dkdt, dkdt2, &
      a, dadt, dadt2 &
      )
      !! # `CMR_aijk`
      !! Calculate the \(a_{ijk}\) tensor of the CMR.
      !!
      !! # Description
      !! \(a_{ijk} = (a_i a_j a_k)^{1/3}(1 - k_{ijk})\), symmetric in its
      !! three indices, together with its first two temperature derivatives.
      real(pr), intent(in) :: ai(:) !! \(a_i\)
      real(pr), intent(in) :: daidt(:) !! \( \frac{da_i}{dT} \)
      real(pr), intent(in) :: daidt2(:) !! \( \frac{d^2a_i}{dT^2} \)
      real(pr), intent(in) :: k(:, :, :) !! \(k_{ijk}\) tensor
      real(pr), intent(in) :: dkdt(:, :, :) !! \(  \frac{d k_{ijk}}{dT} \)
      real(pr), intent(in) :: dkdt2(:, :, :) !! \(  \frac{d^2k_{ijk}}{dT^2} \)
      real(pr), intent(out) :: a(:, :, :) !! \(a_{ijk}\) matrix
      real(pr), intent(out) :: dadt(:, :, :) !! \(\frac{da_{ijk}}{dT}\) matrix
      real(pr), intent(out) :: dadt2(:, :, :) !! \(\frac{d^2a_{ijk}}{dT^2}\) matrix

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: l !! Loop index `l`.
      integer :: nc !! Number of components.

      real(pr) :: aijk_3
      real(pr) :: aaa
      real(pr) :: daaadt

      real(pr), parameter :: third = 1._pr/3._pr

      a = 0
      dadt = 0
      dadt2 = 0

      nc = size(ai)

      do i = 1,nc
         aaa = ai(i) * ai(i) * ai(i)
         aijk_3 = aaa ** third

         a(i, i, i) =  aijk_3

         daaadt = (&
            daidt(i) * ai(i) * ai(i) &
            + ai(i) * daidt(i) * ai(i) &
            + ai(i) * ai(i) * daidt(i) &
            ) * third

         dadt(i, i, i) = daidt(i)
         dadt2(i, i, i) = daidt2(i)

         do j = i,nc
            do l = i,nc
               aaa = ai(i) * ai(j) * ai(l)
               aijk_3 = aaa ** third

               a(i, j, l) =  aijk_3 * (1 - k(i, j, l))

               a(j, i, l) = a(i, j, l)
               a(j, l, i) = a(i, j, l)

               daaadt = (&
                  daidt(i) * ai(j) * ai(l) &
                  + ai(i) * daidt(j) * ai(l) &
                  + ai(i) * ai(j) * daidt(l) &
                  ) * third

               dadt(i, j, l) = - dkdt(i, j, l) * aijk_3 + daaadt / (aaa) * a(i, j, l)
               dadt(j, i, l) = dadt(i, j, l)
               dadt(j, l, i) = dadt(i, j, l)

               dadt2(i, j, l) = &
                  - 2 * dkdt(i, j, l) * daaadt * aijk_3 / aaa &
                  - daidt(i) * daaadt * a(i, j, l) / (aaa * ai(i)) &
                  - daidt(j) * daaadt * a(i, j, l) / (aaa * ai(j)) &
                  - daidt(l) * daaadt * a(i, j, l) / (aaa * ai(l)) &
                  - dkdt2(i, j, l) * aijk_3 &
                  + ((&
                  2 * daidt(i) * daidt(j) * ai(l) &
                  + 2 * daidt(i) * daidt(l) * ai(j) &
                  + 2 * daidt(j) * daidt(l) * ai(i) &
                  + daidt2(i) * ai(j) * ai(l) &
                  + daidt2(j) * ai(i) * ai(l) &
                  + daidt2(l) * ai(j) * ai(i) &
                  )/3._pr) * (a(i, j, l) / aaa) &
                  + (a(i, j, l)*daaadt**2/aaa**2)
               dadt2(j, i, l) = dadt2(i, j, l)
               dadt2(j, l, i) = dadt2(i, j, l)

            end do
         end do
      end do
   end subroutine CMR_aijk

   subroutine CMR_Dmix(n, V, T, &
      a, dadt, dadt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! # CMR_Dmix
      !! Evaluate the ternary cubic attractive mixture parameter and derivatives.
      !!
      !! # Description
      !! \(D = \sum_{ijk} n_i n_j n_k a_{ijk}/n_T\), with its temperature and
      !! composition derivatives. Volume derivatives are returned as zero.
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: a(:, :, :) !! \(a_{ijk}\) tensor (see `CMR_aijk`)
      real(pr), intent(in) :: dadt(:, :, :) !! \(da_{ijk}/dT\)
      real(pr), intent(in) :: dadt2(:, :, :) !! \(d^2a_{ijk}/dT^2\)

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dV}\) (zero)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV}\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :) !! \(\frac{d^2D}{dn_{ij}}\)

      real(pr) :: aux(size(n)) !! Intermediate accumulator.
      real(pr) :: auxij(size(n),size(n))
      real(pr) :: auxT(size(n))
      real(pr) :: auxTij(size(n),size(n))
      real(pr) :: auxT2(size(n)) !! Second-order intermediate quantity used by this procedure.
      real(pr) :: auxT2ij(size(n),size(n))
      real(pr) :: aijk(size(n),size(n),size(n))
      real(pr) :: daijkdT(size(n),size(n),size(n))
      real(pr) :: daijkdT2(size(n),size(n),size(n))

      real(pr) :: totn !! Total number of moles.
      real(pr) :: cum

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: k !! Loop index `k`.
      integer :: nc !! Number of components.

      nc = size(n)

      totn = sum(n)

      ! Initialize derivatives as 0
      D = 0
      dDdT = 0
      dDdT2 = 0
      dDi = 0
      dDidT = 0
      dDij = 0

      ! Not density-dependent
      dDdV = 0
      dDdTV = 0
      dDdV2 = 0

      aux = 0
      auxij = 0
      auxT = 0
      auxT2 = 0
      auxTij = 0

      aijk = a
      daijkdT = dadt
      daijkdT2 = dadt2

      do i=1,nc
         do j=1,nc
            auxT2ij(i, j) = 0
            do k=1,nc
               auxij(i, j) = auxij(i, j) + n(k) * aijk(i, j, k)

               auxTij(i, j) = auxTij(i, j) + n(k) * daijkdT(i, j, k)
               auxT2ij(i, j) = auxT2ij(i, j) + n(k) * daijkdT2(i, j, k)
            end do

            aux(i) = aux(i) + n(j) * auxij(i,j)

            auxT(i) = auxT(i) + n(j)*auxTij(i,j)
            auxT2(i) = auxT2(i) + n(j)*auxT2ij(i,j)
         end do

         D = D + n(i)*aux(i)

         dDdT = dDdT + n(i)*auxT(i)
         dDdT2 = dDdT2 + n(i)*auxT2(i)
      end do

      D = D / totn
      dDdT = dDdT / totn
      dDdT2 = dDdT2 / totn

      do i=1,nc
         dDi(i) = (3*aux(i) - D) / totn
         dDidT(i) = (3 * auxT(i) - dDdT) / totn

         do j=1,i
            dDij(i, j) = (6*auxij(i, j) - dDi(i) - dDi(j)) / totn
            dDij(j, i) = dDij(i, j)
         end do
      end do
   end subroutine CMR_Dmix

   pure subroutine CMR_Bmix(n, bijk, Bmix, dBi, dBij)
      !! # CMR_Bmix
      !! Evaluate the ternary cubic repulsive mixture parameter and derivatives.
      !!
      !! # Description
      !! Contracts the ternary covolume tensor with component mole amounts and
      !! returns the mixture value, gradient, and Hessian.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol]
      real(pr), intent(in) :: bijk(:, :, :) !! Ternary covolume tensor \(b_{ijk}\) [L/mol]
      real(pr), intent(out) :: Bmix !! \(B = \sum_{ijk} n_i n_j n_k b_{ijk}/n_T^2\) [L]
      real(pr), intent(out) :: dBi(:) !! \(\partial B/\partial n_i\)
      real(pr), intent(out) :: dBij(:,:) !! \(\partial^2 B/\partial n_i \partial n_j\)
      real(pr) :: aux(size(n)) !! Intermediate accumulator.
      real(pr) :: auxij(size(n), size(n))

      real(pr) :: totn !! Total number of moles
      real(pr) :: sqn !! totn^2

      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.
      integer :: k !! Loop index `k`.
      integer :: nc !! Number of components.

      nc = size(n)
      totn = sum(n)
      sqn = totn*totn
      Bmix = 0
      aux = 0
      auxij = 0

      do i=1,nc
         do j=1,nc
            do k=1,nc
               auxij(i,j) = auxij(i,j) + n(k)*bijk(i,j,k)
            end do
            aux(i) = aux(i) + n(j)*auxij(i,j)
         end do
         Bmix = Bmix + n(i)*aux(i)
      end do

      Bmix=Bmix/sqn

      do i=1,nc
         dBi(i)=(3*aux(i)-2*totn*Bmix)/sqn
         do j=1,i
            dBij(i,j)=(6*auxij(i,j)-2*(Bmix+totn*dBi(i)+totn*dBi(j)))/sqn
            dBij(j,i)=dBij(i,j)
         end do
      end do
   end subroutine CMR_Bmix
end module yaeos__models_ar_cubic_mixing_base
