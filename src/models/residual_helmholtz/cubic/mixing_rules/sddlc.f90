module yaeos__models_ar_cubic_mixing_sddlc
   use yaeos__constants, only: pr
   use yaeos__models_ar_cubic_quadratic_mixing, only: QMR

   type, extends(QMR) :: sDDLC
      real(pr), allocatable :: q(:) !! Segment size
   contains
      procedure :: Dmix => ddlc_Dmix
   end type sDDLC


contains
   subroutine ddlc_Dmix(&
      self, n, V, T, &
      ai, daidt, daidt2, &
      D, &
      dDdV, dDdT, dDdV2, dDdT2, dDi, dDdTV, dDidV, dDidT, dDij &
      )
      !! s-DDLC D mixing rule including V, n, and T derivatives.
      !!
      !! Corresponds to `DLCandTnVder` in the legacy code.
      !!
      !! The local-composition D is:
      !! \[
      !!   D = \frac{\sum_{ij} n_i n_j E_{ij}(V,T)\,a_{ij}(T)}{\Lambda(V,T,n)}
      !! \]
      !! where \(E_{ij} = \exp\!\left(\frac{a_{ij}\sum_k n_k q_k}{q_i q_j RTV}\right)\)
      !! and \(\Lambda = \sum_{ij} s_i s_j E_{ij}\), with segment fractions
      !! \(s_i = n_i q_i / \sum_k n_k q_k\).
      class(DDLC_MixRule), intent(in)  :: self
      real(pr),            intent(in)  :: n(:), V, T
      real(pr),            intent(in)  :: ai(:), daidt(:), daidt2(:)
      real(pr),            intent(out) :: D, dDdT, dDdT2, dDdV, dDdV2, dDdTV
      real(pr),            intent(out) :: dDi(:), dDidT(:), dDidT(:), dDij(:,:)

      integer  :: i, j, k, m, nc
      real(pr) :: RTV, sumnq

      ! -- aij and its T-derivatives --
      real(pr) :: aij(size(n),size(n))
      real(pr) :: daijdT(size(n),size(n)), daijdT2(size(n),size(n))

      ! -- Segment fractions and composition derivatives --
      real(pr) :: sfrac(size(n))
      real(pr) :: dsfdn(size(n),size(n))
      real(pr) :: d2sfdn2(size(n),size(n),size(n))

      ! -- E matrix and all its derivatives --
      real(pr) :: rata(size(n),size(n))
      real(pr) :: E(size(n),size(n))
      real(pr) :: dEdn(size(n),size(n),size(n))
      real(pr) :: d2Edn2(size(n),size(n),size(n),size(n))
      real(pr) :: dEdT(size(n),size(n)), d2EdT2(size(n),size(n))
      real(pr) :: dEdV(size(n),size(n)), d2EdV2(size(n),size(n))
      real(pr) :: d2EdnT(size(n),size(n),size(n))
      real(pr) :: d2EdnV(size(n),size(n),size(n))
      real(pr) :: d2EdVT(size(n),size(n))

      ! -- Denominator Den and its derivatives --
      real(pr) :: Den
      real(pr) :: dDendn(size(n)), d2Dendn2(size(n),size(n))
      real(pr) :: d2DendnT(size(n)), d2DendnV(size(n))
      real(pr) :: dDendV, d2DendV2, d2DendVT, dDendT, d2DendT2

      ! -- Loop temporaries --
      real(pr) :: aux, auxV, auxVV, auxT, auxTT, auxVT
      real(pr) :: auxjmi, auxjmiT, auxjmiV
      real(pr) :: auxdfE, auxd2fE, auxdfdE, auxfd2E, auxd2Ek
      real(pr) :: auxTij, auxdE, auxd2E

      nc  = size(n)
      RTV = R*T*V

      ! -----------------------------------------------------------------------
      ! Build aij matrix (reusing the shared helper)
      ! -----------------------------------------------------------------------
      call self%aij(self, T, ai, daidt, daidt2, aij, daijdt, daijdt2)

      ! -----------------------------------------------------------------------
      ! Segment fractions and their n-derivatives
      ! -----------------------------------------------------------------------
      sumnq = sum(n*self%q)

      do i = 1, nc
         sfrac(i) = n(i)*self%q(i)/sumnq
         do j = 1, nc
            dsfdn(i,j) = -sfrac(i)*self%q(j)/sumnq
            if (i == j) dsfdn(i,i) = dsfdn(i,i) + self%q(i)/sumnq
            do k = j, nc
               d2sfdn2(i,j,k) = 2.0_pr*sfrac(i)*self%q(j)*self%q(k)/sumnq**2
               if (i == j) d2sfdn2(i,i,k) = d2sfdn2(i,i,k) &
                  - self%q(i)*self%q(k)/sumnq**2
               if (i == k) d2sfdn2(i,j,i) = d2sfdn2(i,j,i) &
                  - self%q(i)*self%q(j)/sumnq**2
               d2sfdn2(i,k,j) = d2sfdn2(i,j,k)
            end do
         end do
      end do

      ! -----------------------------------------------------------------------
      ! E matrix, V/n/T-derivatives, and accumulation of Den
      ! -----------------------------------------------------------------------
      Den        = 0.0_pr
      dDendn     = 0.0_pr
      d2DendnT   = 0.0_pr
      d2DendnV   = 0.0_pr
      dDendV     = 0.0_pr
      d2DendV2   = 0.0_pr
      d2DendVT   = 0.0_pr
      dDendT     = 0.0_pr
      d2DendT2   = 0.0_pr
      d2Edn2     = 0.0_pr
      dEdT       = 0.0_pr
      d2EdT2     = 0.0_pr
      d2EdVT     = 0.0_pr
      d2EdnT     = 0.0_pr

      do k = 1, nc
         do m = k, nc
            rata(k,m) = aij(k,m)/(self%q(k)*self%q(m)*RTV)
            E(k,m)    = exp(rata(k,m)*sumnq)
            E(m,k)    = E(k,m)

            ! -- n-derivatives of E --
            do i = 1, nc
               dEdn(k,m,i)     = E(k,m)*rata(k,m)*self%q(i)
               dEdn(m,k,i)     = dEdn(k,m,i)
               d2Edn2(k,m,i,i) = dEdn(k,m,i)**2/E(k,m)
               d2Edn2(m,k,i,i) = d2Edn2(k,m,i,i)
               do j = 1, i-1
                  d2Edn2(k,m,i,j) = dEdn(k,m,i)*dEdn(k,m,j)/E(k,m)
                  d2Edn2(k,m,j,i) = d2Edn2(k,m,i,j)
                  d2Edn2(m,k,i,j) = d2Edn2(k,m,i,j)
                  d2Edn2(m,k,j,i) = d2Edn2(k,m,i,j)
               end do
            end do

            ! -- V-derivatives of E --
            dEdV(k,m)     = -E(k,m)*rata(k,m)*sumnq/V
            d2EdV2(k,m)   =  dEdV(k,m)**2/E(k,m) - 2.0_pr*dEdV(k,m)/V
            d2EdnV(k,m,:) =  self%q(:)*dEdV(k,m)/sumnq &
               + dEdV(k,m)*dEdn(k,m,:)/E(k,m)
            dEdV(m,k)     = dEdV(k,m)
            d2EdV2(m,k)   = d2EdV2(k,m)
            d2EdnV(m,k,:) = d2EdnV(k,m,:)

            ! -- T-derivatives of E (computed only when requested) --
            if (calc_T) then
               dEdT(k,m)    = E(k,m)*rata(k,m)*sumnq &
                  *(daijdT(k,m)/aij(k,m) - 1.0_pr/T)
               d2EdT2(k,m)  = dEdT(k,m)**2/E(k,m) - 2.0_pr*dEdT(k,m)/T &
                  + E(k,m)*rata(k,m)*sumnq*daijdT2(k,m)/aij(k,m)
               d2EdVT(k,m)  = dEdV(k,m)* &
                  (dEdT(k,m)/E(k,m) + daijdT(k,m)/aij(k,m) - 1.0_pr/T)
               d2EdnT(k,m,:)= self%q(:)*dEdT(k,m)/sumnq &
                  + dEdT(k,m)*dEdn(k,m,:)/E(k,m)
               dEdT(m,k)    = dEdT(k,m)
               d2EdT2(m,k)  = d2EdT2(k,m)
               d2EdVT(m,k)  = d2EdVT(k,m)
               d2EdnT(m,k,:)= d2EdnT(k,m,:)
            end if

            ! -- Accumulate Den and its derivatives --
            if (k == m) then
               Den      = Den      + sfrac(k)**2*E(k,k)
               dDendV   = dDendV   + sfrac(k)**2*dEdV(k,k)
               d2DendV2 = d2DendV2 + sfrac(k)**2*d2EdV2(k,k)
               if (calc_T) then
                  dDendT   = dDendT   + sfrac(k)**2*dEdT(k,k)
                  d2DendT2 = d2DendT2 + sfrac(k)**2*d2EdT2(k,k)
                  d2DendVT = d2DendVT + sfrac(k)**2*d2EdVT(k,k)
               end if
               do i = 1, nc
                  dDendn(i)   = dDendn(i) &
                     + 2.0_pr*sfrac(k)*dsfdn(k,i)*E(k,k) &
                     + sfrac(k)**2*dEdn(k,k,i)
                  d2DendnV(i) = d2DendnV(i) &
                     + 2.0_pr*sfrac(k)*dsfdn(k,i)*dEdV(k,k) &
                     + sfrac(k)**2*d2EdnV(k,k,i)
                  if (calc_T) then
                     d2DendnT(i) = d2DendnT(i) &
                        + 2.0_pr*sfrac(k)*dsfdn(k,i)*dEdT(k,k) &
                        + sfrac(k)**2*d2EdnT(k,k,i)
                  end if
               end do
            else
               Den      = Den      + 2.0_pr*sfrac(k)*sfrac(m)*E(k,m)
               dDendV   = dDendV   + 2.0_pr*sfrac(k)*sfrac(m)*dEdV(k,m)
               d2DendV2 = d2DendV2 + 2.0_pr*sfrac(k)*sfrac(m)*d2EdV2(k,m)
               if (calc_T) then
                  dDendT   = dDendT   + 2.0_pr*sfrac(k)*sfrac(m)*dEdT(k,m)
                  d2DendT2 = d2DendT2 + 2.0_pr*sfrac(k)*sfrac(m)*d2EdT2(k,m)
                  d2DendVT = d2DendVT + 2.0_pr*sfrac(k)*sfrac(m)*d2EdVT(k,m)
               end if
               do i = 1, nc
                  dDendn(i)   = dDendn(i) &
                     + 2.0_pr*sfrac(k)*sfrac(m)*dEdn(k,m,i) &
                     + 2.0_pr*(sfrac(k)*dsfdn(m,i) + sfrac(m)*dsfdn(k,i))*E(k,m)
                  d2DendnV(i) = d2DendnV(i) &
                     + 2.0_pr*sfrac(k)*sfrac(m)*d2EdnV(k,m,i) &
                     + 2.0_pr*(sfrac(k)*dsfdn(m,i)+sfrac(m)*dsfdn(k,i))*dEdV(k,m)
                  if (calc_T) then
                     d2DendnT(i) = d2DendnT(i) &
                        + 2.0_pr*sfrac(k)*sfrac(m)*d2EdnT(k,m,i) &
                        + 2.0_pr*(sfrac(k)*dsfdn(m,i)+sfrac(m)*dsfdn(k,i))*dEdT(k,m)
                  end if
               end do
            end if
         end do
      end do

      ! -----------------------------------------------------------------------
      ! Second composition derivative of Den
      ! -----------------------------------------------------------------------
      d2Dendn2 = 0.0_pr
      do i = 1, nc
         do j = i, nc
            do k = 1, nc
               auxdfE  = 0.0_pr
               auxd2fE = 0.0_pr
               auxdfdE = 0.0_pr
               auxfd2E = 0.0_pr
               do m = 1, nc
                  auxdfE  = auxdfE  + dsfdn(m,i)*E(k,m)
                  auxd2fE = auxd2fE + d2sfdn2(m,i,j)*E(k,m)
                  auxdfdE = auxdfdE &
                     + dsfdn(m,i)*dEdn(k,m,j) + dsfdn(m,j)*dEdn(k,m,i)
                  auxfd2E = auxfd2E + sfrac(m)*d2Edn2(k,m,i,j)
               end do
               d2Dendn2(i,j) = d2Dendn2(i,j) &
                  + 2.0_pr*dsfdn(k,j)*auxdfE &
                  + sfrac(k)*(2.0_pr*auxd2fE + 2.0_pr*auxdfdE + auxfd2E)
            end do
            d2Dendn2(j,i) = d2Dendn2(i,j)
         end do
      end do

      ! -----------------------------------------------------------------------
      ! Accumulate D and its n/V/T derivatives (numerator terms, before /Den)
      ! -----------------------------------------------------------------------
      D      = 0.0_pr
      dDdT   = 0.0_pr
      dDdT2  = 0.0_pr
      dDdV     = 0.0_pr
      dDdV2    = 0.0_pr
      dDdTV    = 0.0_pr
      dDi    = 0.0_pr
      dDidT    = 0.0_pr
      dDidT  = 0.0_pr

      do i = 1, nc
         aux   = 0.0_pr; auxV  = 0.0_pr; auxVV = 0.0_pr
         auxT  = 0.0_pr; auxTT = 0.0_pr; auxVT = 0.0_pr
         do j = 1, nc
            auxjmi  = 0.0_pr; auxjmiT = 0.0_pr; auxjmiV = 0.0_pr
            do m = 1, nc
               auxjmi  = auxjmi  + n(m)*dEdn(m,j,i)*aij(m,j)
               auxjmiV = auxjmiV + n(m)*d2EdnV(m,j,i)*aij(m,j)
               if (calc_T) then
                  auxjmiT = auxjmiT + n(m)*( &
                     dEdn(m,j,i)*daijdT(m,j) + d2EdnT(m,j,i)*aij(m,j))
               end if
            end do

            dDi(i) = dDi(i) + 2.0_pr*n(j)*E(i,j)*aij(i,j) + n(j)*auxjmi
            dDidT(i) = dDidT(i) + 2.0_pr*n(j)*dEdV(i,j)*aij(i,j) + n(j)*auxjmiV

            if (calc_T) then
               auxTij  = E(i,j)*daijdT(i,j) + dEdT(i,j)*aij(i,j)
               auxT    = auxT + n(j)*auxTij
               dDidT(i)= dDidT(i) + 2.0_pr*n(j)*auxTij + n(j)*auxjmiT
               auxTT   = auxTT + n(j)*( &
                  E(i,j)*daijdT2(i,j) + 2.0_pr*dEdT(i,j)*daijdT(i,j) &
                  + d2EdT2(i,j)*aij(i,j))
               auxVT   = auxVT + n(j)*( &
                  d2EdVT(i,j)*aij(i,j) + dEdV(i,j)*daijdT(i,j))
            end if

            aux   = aux   + n(j)*E(i,j)*aij(i,j)
            auxV  = auxV  + n(j)*dEdV(i,j)*aij(i,j)
            auxVV = auxVV + n(j)*d2EdV2(i,j)*aij(i,j)
         end do

         D   = D   + n(i)*aux
         dDdV  = dDdV  + n(i)*auxV
         dDdV2 = dDdV2 + n(i)*auxVV
         if (calc_T) then
            dDdT  = dDdT  + n(i)*auxT
            dDdT2 = dDdT2 + n(i)*auxTT
            dDdTV   = dDdTV   + n(i)*auxVT
         end if
      end do

      ! -----------------------------------------------------------------------
      ! Complete D and its derivatives (divide by Den, apply quotient rule)
      ! -----------------------------------------------------------------------
      D   =  D/Den
      dDdV  = (dDdV  - D*dDendV)/Den
      dDdV2 = (dDdV2 - 2.0_pr*dDdV*dDendV - D*d2DendV2)/Den

      dDi(1:nc) = (dDi(1:nc)  - D*dDendn(1:nc))/Den
      dDidT(1:nc) = (dDidT(1:nc)  - dDdV*dDendn(1:nc) - D*d2DendnV(1:nc) &
         - dDi(1:nc)*dDendV)/Den

      dDdT  = (dDdT  - D*dDendT)/Den
      dDdT2 = (dDdT2 - 2.0_pr*dDdT*dDendT - D*d2DendT2)/Den
      dDdTV   = (dDdTV   - dDdT*dDendV - dDdV*dDendT - D*d2DendVT)/Den
      dDidT(1:nc) = (dDidT(1:nc) - dDdT*dDendn(1:nc) - D*d2DendnT(1:nc) &
         - dDi(1:nc)*dDendT)/Den

      ! -----------------------------------------------------------------------
      ! Second composition derivative of D (= dDij)
      ! -----------------------------------------------------------------------
      do i = 1, nc
         do j = i, nc
            auxdE  = 0.0_pr
            auxd2E = 0.0_pr
            do k = 1, nc
               auxd2Ek = 0.0_pr
               do m = 1, nc
                  auxd2Ek = auxd2Ek + n(m)*d2Edn2(k,m,i,j)*aij(k,m)
               end do
               auxdE  = auxdE + n(k)*( &
                  dEdn(i,k,j)*aij(i,k) + dEdn(k,j,i)*aij(k,j))
               auxd2E = auxd2E + n(k)*auxd2Ek
            end do
            dDij(i,j) = (2.0_pr*E(i,j)*aij(i,j) + 2.0_pr*auxdE + auxd2E &
               - dDi(j)*dDendn(i) - dDi(i)*dDendn(j) &
               - D*d2Dendn2(i,j)) / Den
            dDij(j,i) = dDij(i,j)
         end do
      end do
   end subroutine ddlc_Dmix
end module yaeos__models_ar_cubic_mixing_sddlc

