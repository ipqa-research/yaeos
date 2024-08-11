module yaeos__phase_equilibria_stability
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
   use yaeos__models_ar, only: ArModel
   implicit none

contains

   real(pr) function tm(model, z, w, P, T, d, dtpd)
      !! # Alternative formulation of tangent-plane-distance
      !! Michelsen's modified \(tpd\) function, \(tm\).
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
      !!  tm = tpd(model, z, w, P, T)
      !!  ---------------------------
      !! ```
      !!
      !! ## Using precalculated trial-phase data
      !! It is possible to calculate externaly the `d_i` vector and use it for
      !! later calculations.
      !! ```fortran
      !! call fugacity_tp(&
      !!   model, z, T=T, P=P, V=Vz, root_type="stable", lnphip=lnphi_z&
      !! )
      !! lnphi_z = lnphi_z - log(P)
      !! di = log(z) + lnphi_z
      !! tm = tpd(model, z, w, P, T, d=di)
      !!  ---------------------------
      !! ```
      !!
      !! # References
      !! 1. Thermodynamic Models: Fundamental and Computational Aspects, Michael L.
      !! Michelsen, Jørgen M. Mollerup. Tie-Line Publications, Denmark (2004)
      !! [doi](http://dx.doi.org/10.1016/j.fluid.2005.11.032)
      class(ArModel), intent(in) :: model !! Thermodynamic model
      real(pr), intent(in) :: z(:) !! Feed composition
      real(pr), intent(in) :: w(:) !! Test-phase mole numbers vector
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), optional, intent(in) :: d(:) !! \(d_i\) vector
      real(pr), optional, intent(out) :: dtpd(:)

      real(pr) :: di(size(z)), vz, vw
      real(pr) :: lnphi_z(size(z)), lnphi_w(size(z))

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


      ! tpd = sum(w * (log(w) + lnphi_w - di))
      tm = 1 + sum(w * (log(w) + lnPhi_w - di - 1))

      if (present(dtpd)) then
         dtpd = log(w) + lnPhi_w - di
      end if
   end function tm

   subroutine min_tpd(model, z, P, T, mintpd, w, all_minima)
      use yaoes__optimizers_nlopt_wrap, only: NLOPTWrapper
      class(ArModel) :: model
      real(pr), intent(in) :: z(:) !! Feed composition
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: w(:) !! Trial composition
      real(pr), intent(out) :: mintpd !! Minimal value of \(tm\)
      real(pr), optional, intent(out) :: all_minima(:, :) 
         !! All the found minima
      type(NLOPTWrapper) :: opt

      real(pr) :: dx(size(w))
      real(pr) :: lnphi_z(size(z)), di(size(z))

      real(pr) :: mins(size(w)), ws(size(w), size(w)), V
      integer :: i

      integer :: stat

      dx = 0.001_pr

      ! Calculate feed di
      call model%lnphi_pt(z, T=T, P=P, V=V, root_type="stable", lnPhi=lnPhi_z)
      di = log(z) + lnphi_z

      ! ==============================================================
      ! Minimize for each component using each quasi-pure component
      ! as initialization.
      ! --------------------------------------------------------------
      opt%parameter_step = dx
      !$OMP PARALLEL DO PRIVATE(i, w, mintpd, stat) SHARED(opt, ws, mins)
      do i=1,size(w)
         w = 0.001_pr
         w(i) = 0.999_pr
         call opt%optimize(obj_func, w, mintpd)
         mins(i) = mintpd
         ws(i, :) = w
      end do
      !$OMP END PARALLEL DO

      i = minloc(mins, dim=1)
      mintpd = mins(i)
      w = ws(i, :)

      if(present(all_minima)) all_minima = ws
   contains
      subroutine obj_func(X, F, dF, data)
         real(pr), intent(in) :: X(:)
         real(pr), intent(out) :: F
         real(pr), optional, intent(out) :: dF(:)
         class(*), optional, intent(in out) :: data
         F = tm(model, z, x, P, T, d=di)
      end subroutine
   end subroutine min_tpd
end module yaeos__phase_equilibria_stability
