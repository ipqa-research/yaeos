module yaeos__equilibria_stability
   !! # Phase Stability module
   !! Phase stability related calculations.
   !!
   !! # Description
   !! Contains the basics rotuines to make phase stability analysis for
   !! phase-equilibria detection.
   !!
   !! - `tpd(model, z, w, P, T)`: reduced Tangent-Plane-Distance
   !! - `min_tpd(model, z, P, T, mintpd, w)`: Find minimal tpd for a multicomponent mixture
   !!
   !! # Examples
   !!
   !! ```fortran
   !!   ! Obtain the minimal tpd for a binary mixture at \(z_1 = 0.13\)
   !!   model = PengRobinson76(tc, pc, ac, kij, lij)
   !!
   !!   z = [0.13, 1-0.13]
   !!   w = [0.1, 0.9]
   !!
   !!   P = 45.6_pr
   !!   T = 190._pr
   !!
   !!   z = z/sum(z)
   !! -----------------------------------------------
   !! ```
   !!
   !! # References
   !! 1. Thermodynamic Models: Fundamental and Computational Aspects, Michael L.
   !! Michelsen, Jørgen M. Mollerup. Tie-Line Publications, Denmark (2004)
   !! [doi](http://dx.doi.org/10.1016/j.fluid.2005.11.032)
   use yaeos__constants, only: pr, r
   use yaeos__models_base, only: BaseModel
   use yaeos__models_ar, only: ArModel
   use yaeos__models_ge, only: GeModel
   implicit none

contains

   real(pr) function tm(model, z, w, P, T, d, dtpd)
      !! # tm
      !! Evaluate Michelsen's modified reduced tangent-plane distance.
      !!
      !! # Description
      !! Alternative formulation of the reduced tangent plane \(tpd\) function,
      !! where the test phase is defined in moles, which enables for unconstrained
      !! minimization.
      !! \[
      !!   tm(W) = 1 + \sum_i W_i (\ln W_i + \ln \phi_i(W) - d_i - 1)
      !! \]
      !!
      !! # Examples
      !!
      !! ## Calculation of `tm`
      !! ```fortran
      !!  tpd_value = tm(model, z, w, P, T)
      !! ```
      !!
      !! ## Using precalculated feed data
      !! It is possible to calculate externally the feed `d_i` vector and use it for
      !! later calculations.
      !! ```fortran
      !! call model%lnphi_pt(z, T=T, P=P, V=Vz, root_type="stable", lnPhi=lnphi_z)
      !! d = log(z) + lnphi_z
      !! tpd_value = tm(model, z, w, P, T, d=d)
      !! ```
      !!
      !! # References
      !! 1. Thermodynamic Models: Fundamental and Computational Aspects, Michael L.
      !! Michelsen, Jørgen M. Mollerup. Tie-Line Publications, Denmark (2004)
      !! [doi](http://dx.doi.org/10.1016/j.fluid.2005.11.032)
      class(BaseModel), intent(in) :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Feed composition
      real(pr), intent(in) :: w(:) !! Test-phase mole numbers vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(in) :: d(:) !! \(d_i\) vector
      real(pr), optional, intent(out) :: dtpd(:) !! Gradient of `tm` with respect to trial-phase mole amounts.

      real(pr) :: di(size(z))
      real(pr) :: vz
      real(pr) :: vw
      real(pr) :: lnphi_z(size(z))
      real(pr) :: lnphi_w(size(z))

      select type (model)
       class is (ArModel)
         call model%lnphi_pt(&
            w, T=T, P=P, V=Vw, root_type="stable", lnPhi=lnPhi_w &
            )
         if (.not. present(d)) then
            call model%lnphi_pt(&
               z, T=T, P=P, V=Vz, root_type="stable", lnPhi=lnPhi_z&
               )
            di = log(z) + lnphi_z
         else
            di = d
         end if

       class is (GeModel)
         call model%ln_activity_coefficient(w, T=T, lngamma=lnPhi_w)

         if (.not. present(d)) then
            call model%ln_activity_coefficient(z, T=T, lngamma=lnPhi_z)
            di = log(z) + lnphi_z
         else
            di = d
         end if
      end select

      ! tpd = sum(w * (log(w) + lnphi_w - di))
      tm = 1 + sum(w * (log(w) + lnPhi_w - di - 1))

      if (present(dtpd)) then
         dtpd = log(w) + lnPhi_w - di
      end if
   end function tm

   subroutine min_tpd(model, z, P, T, mintpd, w, all_minima)
      !! # min_tpd
      !!
      !! Find the minimum of the reduced tangent-plane distance.
      !!
      !! ## Description
      !! Minimizes Michelsen's modified tangent-plane distance with respect to
      !! the trial-phase mole numbers. A negative minimum indicates that the
      !! feed is unstable. The routine can optionally return all the minima
      !! found from different starting points.
      class(BaseModel), target :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Feed composition
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: w(:) !! Trial composition
      real(pr), intent(out) :: mintpd !! Minimal value of \(tm\)
      real(pr), optional, intent(out) :: all_minima(:, :)
      !! All the found minima

      real(pr) :: dx(size(w))
      real(pr) :: lnphi_z(size(z))
      real(pr) :: di(size(z))

      real(pr) :: lnphi_w(size(w))
      real(pr) :: dw(size(w))
      real(pr) :: mins(size(w))
      real(pr) :: ws(size(w), size(w))
      real(pr) :: V
      integer :: i !! Loop index `i`.
      integer :: j !! Loop index `j`.

      integer :: nc !! Number of components.
      integer :: stat
      integer :: max_iters
      integer :: iters

      nc = size(z)

      dx = 0.001_pr

      ! Calculate feed di
      select type (model)
       class is (ArModel)
         call model%lnphi_pt(z, T=T, P=P, V=V, root_type="stable", lnPhi=lnPhi_z)
       class is (GeModel)
         call model%ln_activity_coefficient(z, T=T, lngamma=lnPhi_z)
      end select
      di = log(z) + lnphi_z

      ! ==============================================================
      ! Minimize for each component using each quasi-pure component
      ! as initialization.
      ! --------------------------------------------------------------
      max_iters = 1000
      mins = 10
      do i=1,nc
         iters = 0
         ! w = 1e-15
         ! w(i) = 1 - 1e-15
         ! w = w/sum(w)
         w = 1
         w(i) = 1000
         dw = 100
         do while(maxval(abs(dw)) > 1e-8 .and. abs(mins(i)) > 1e-4 .and. iters < max_iters)
            iters = iters + 1
            select type (model)
             class is (ArModel)
               call model%lnphi_pt(w, T=T, P=P, V=V, root_type="stable", lnPhi=lnPhi_w)
             class is (GeModel)
               call model%ln_activity_coefficient(w, T=T, lngamma=lnPhi_w)
            end select

            dw = exp(di - lnphi_w) - w
            do while(any(dw + w < 0))
               dw = dw/2
            end do

            mins(i) = 1 + sum(w * (log(w) + lnPhi_w - di - 1))
            w = w + dw
         end do
         w = w/sum(w)
         mins(i) = tm(model, z, w, P, T, d=di)
         ws(i, :) = w
      end do

      i = minloc(mins, dim=1)
      mintpd = mins(i)
      w = ws(i, :)

      if(present(all_minima)) then
         do i=1,nc
            all_minima(i, :nc) = ws(i, :)
            all_minima(i, nc+1) = mins(i)
         end do
      end if
   end subroutine min_tpd
end module yaeos__equilibria_stability
