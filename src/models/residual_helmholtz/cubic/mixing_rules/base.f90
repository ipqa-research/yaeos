module yaeos__models_ar_cubic_mixing_base
   !! # Mixing rules core math
   !! Procedures of the core calculations of CubicEoS mixing rules.
   !!
   !! # Description
   !! This module holds all the basic math to use mixing rules in other codes.
   !! Keeping it simple and accesible.
   !!
   !! # Examples
   !!
   !! ```fortran
   !! bi = [0.2, 0.3]
   !! lij = reshape([0.0, 0.2, 0.2, 0], [2,2])
   !!
   !! ! Calculate B parameter with Quadratric Mixing Rules.
   !! call bmix_qmr(n, bi, lij, b, dbi, dbij)
   !!
   !! ```
   !!
   !! # References
   use yaeos__constants, only: pr, solving_volume
   implicit none
contains

   pure subroutine bmix_linear(n, bi, b, dbi, dbij)
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(out) :: b, dbi(:), dbij(:, :)

      b = sum(n*bi)
      dbi = bi
      dbij = 0
   end subroutine bmix_linear

   pure subroutine bmix_qmr(n, bi, lij, b, dbi, dbij)
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bi(:)
      real(pr), intent(in) :: lij(:, :)
      real(pr), intent(out) :: b, dbi(:), dbij(:, :)

      real(pr) :: bij(size(n), size(n))

      real(pr) :: totn, aux(size(n))

      integer :: i, j, nc

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
      !! RKPR \(\delta_1\) parameter mixing rule.
      !!
      !! The RKPR EoS doesn't have a constant \(\delta_1\) value for each
      !! component, so a proper mixing rule should be provided. A linear
      !! combination is used.
      !!
      !! \[
      !!     \Delta_1 = \sum_i^N n_i \delta_{1i}
      !! \]
      !!
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: d1i(:)
      real(pr), intent(out) :: D1
      real(pr), intent(out) :: dD1i(:)
      real(pr), intent(out) :: dD1ij(:, :)

      integer :: i, j, nc
      real(pr) :: totn

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
      !! Infinite pressure limit parameter \(\Lambda\)
      !!
      !! \[
      !! \Lambda = \frac{1}{\delta_1 - \delta_2} \ln \frac{1 + \delta_1}{1 + \delta_2}
      !! \]
      integer, intent(in) :: nc
      real(pr), intent(in) :: d1
      real(pr), optional, intent(in) :: dd1i(nc)
      real(pr), optional, intent(in) :: dd1ij(nc, nc)
      real(pr), intent(out) :: L
      real(pr), optional, intent(out) :: dLi(nc)
      real(pr), optional, intent(out) :: dLij(nc, nc)

      real(pr) :: f, g, h
      real(pr), dimension(nc) :: df, dg, dh
      real(pr), dimension(nc, nc) :: d2f, d2g, d2h

      integer :: i, j

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
      !! # Examples
      !!
      !! # References
      !!
      real(pr), intent(in) :: T, n(:)
      real(pr), intent(in) :: bi(:) !! Covolume parameter
      real(pr), intent(in) :: B !! mixture covolume parameter
      real(pr), intent(in) :: dBi(:), dBij(:, :)
      real(pr), intent(in) :: D1, dD1i(:), dD1ij(:, :)
      real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
      real(pr), intent(in) :: Ge, GeT, GeT2
      real(pr), intent(in) :: Gen(:), GeTn(:), Gen2(:, :)
      real(pr), intent(out) :: D, dDdT, dDdT2, dDi(:), dDidT(:), dDij(:, :)

      real(pr) :: f, fdt, fdt2, fdi(size(n)), fdit(size(n)), fdij(size(n), size(n))
      real(pr) :: totn !! Total number of moles

      integer :: i, j, nc
      real(pr) :: L, dL(size(n)), dL2(size(n), size(n))

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
      real(pr), intent(in) :: ai(:) !! \(\a_i\)
      real(pr), intent(in) :: daidt(:) !! \( \frac{d\a_i}{dT} \)
      real(pr), intent(in) :: daidt2(:)!! \( \frac{d^2\a_i}{dT^2} \)
      real(pr), intent(in) :: k(:, :, :) !! k_{ijk} matrix
      real(pr), intent(in) :: dkdt(:, :, :) !! \(  \frac{d k_{ijk}}{dT} \)
      real(pr), intent(in) :: dkdt2(:, :, :) !! \(  \frac{d^2k_{ijk}}{dT^2} \)
      real(pr), intent(out) :: a(:, :, :) !! \(a_{ijk}\) matrix
      real(pr), intent(out) :: dadt(:, :, :) !! \(\frac{da_{ijk}}{dT}\) matrix
      real(pr), intent(out) :: dadt2(:, :, :) !! \(\frac{d^2a_{ijk}{dT^2}\) matrix

      integer :: i, j, l, nc

      real(pr) :: aijk_3
      real(pr) :: aaa, daaadt

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
      real(pr), intent(in) :: V !! Volume [L] (unused)
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(in) :: n(:) !! Moles vector [mol]
      real(pr), intent(in) :: a(:, :, :) !!
      real(pr), intent(in) :: dadt(:, :, :) !!
      real(pr), intent(in) :: dadt2(:, :, :) !!

      real(pr), intent(out) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
      real(pr), intent(out) :: dDdV !! \(\frac{dD}{dT}\)
      real(pr), intent(out) :: dDdT !! \(\frac{dD}{dV}\)
      real(pr), intent(out) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
      real(pr), intent(out) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
      real(pr), intent(out) :: dDdTV !! \(\frac{d^2D}{dTV\)
      real(pr), intent(out) :: dDi(:) !! \(\frac{dD}{dn_i}\)
      real(pr), intent(out) :: dDidV(:) !! \(\frac{d^2D}{dVn_i}\)
      real(pr), intent(out) :: dDidT(:) !! \(\frac{d^2D}{dTn_i}\)
      real(pr), intent(out) :: dDij(:, :)!! \(\frac{d^2D}{dn_{ij}}\)

      real(pr) :: aux(size(n)),auxij(size(n),size(n))
      real(pr) :: auxT(size(n)),auxTij(size(n),size(n)),auxT2(size(n)),auxT2ij(size(n),size(n))
      real(pr) :: aijk(size(n),size(n),size(n)),daijkdT(size(n),size(n),size(n)),daijkdT2(size(n),size(n),size(n))

      real(pr) :: totn
      real(pr) :: cum

      integer :: i, j, k, nc

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
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: bijk(:, :, :)
      real(pr), intent(out) :: Bmix
      real(pr), intent(out) :: dBi(:), dBij(:,:)
      real(pr) :: aux(size(n)), auxij(size(n), size(n))

      real(pr) :: totn !! Total number of moles
      real(pr) :: sqn !! totn^2

      integer :: i, j, k, nc

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
