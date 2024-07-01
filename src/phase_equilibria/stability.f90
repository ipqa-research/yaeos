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
   use yaeos__thermoprops, only: fugacity_vt, fugacity_tp
   use yaeos__models_ar, only: ArModel
   implicit none

contains

   real(pr) function tpd(model, z, w, P, T, d, dtpd, lnphiw, outd)
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
      real(pr), optional, intent(out) :: lnphiw(:)
      real(pr), optional, intent(out) :: outd(:)

      real(pr) :: di(size(z)), vz, vw
      real(pr) :: lnphi_z(size(z)), lnphi_w(size(z))

      call fugacity_tp(model, w, T=T, P=P, V=Vw, root_type="stable", lnphip=lnphi_w)
      lnphi_w = lnphi_w - log(P)

      if (.not. present(d)) then
         call fugacity_tp(&
            model, z, T=T, P=P, V=Vz, root_type="stable", lnphip=lnphi_z&
         )
         lnphi_z = lnphi_z - log(P)
         di = log(z) + lnphi_z
      end if

      tpd = 1 + sum(w * (log(w) + lnphi_w - di - 1))

      if (present(dtpd)) then
         dtpd = log(w) + lnphi_w - di
      end if

      if (present(lnphiw)) lnphiw = lnphi_w
      if (present(outd)) outd = di
   end function tpd

   subroutine min_tpd(model, z, P, T, mintpd, w, all_minima)
      use nlopt_wrap, only: create, destroy, nlopt_opt, nlopt_algorithm_enum
      use nlopt_callback, only: nlopt_func, create_nlopt_func
      class(ArModel) :: model
      real(pr), intent(in) :: z(:) !! Feed composition
      real(pr), intent(in) :: P !! Pressure [bar]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr), intent(out) :: w(:) !! Trial composition
      real(pr), intent(out) :: mintpd !! Minimal value of \(tm\)
      real(pr), optional, intent(out) :: all_minima(:, :) 
         !! All the found minima

      real(pr) :: dx(size(w))
      real(pr) :: lnphi_z(size(z)), di(size(z))

      real(pr) :: mins(size(w)), ws(size(w), size(w))
      integer :: i

      type(nlopt_opt) :: opt !! Optimizer
      type(nlopt_func) :: f !! Function to optimize

      integer :: stat

      f = create_nlopt_func(foo)
      dx = 0.001_pr

      
      ! Calculate feed di
      call fugacity_tp(&
         model, z, T=T, P=P, root_type="stable", lnphip=lnphi_z&
      )
      di = log(z) + lnphi_z - log(P)


      ! ==============================================================
      ! Setup optimizer
      ! --------------------------------------------------------------
      ! opt = nlopt_opt(nlopt_algorithm_enum%LN_NELDERMEAD, size(w))
      ! opt = nlopt_opt(nlopt_algorithm_enum%LD_TNEWTON, size(w))
      opt = nlopt_opt(nlopt_algorithm_enum%LN_NELDERMEAD, size(w))
      call opt%set_ftol_rel(0.001_pr)
      call opt%set_ftol_abs(0.00001_pr)
      call opt%set_min_objective(f)
      call opt%set_initial_step(dx)

      ! ==============================================================
      ! Minimize for each component using each quasi-pure component
      ! as initialization.
      ! --------------------------------------------------------------
      do i=1,size(w)
         w = 0.001_pr
         w(i) = 0.999_pr
         call opt%optimize(w, mintpd, stat)
         mins(i) = mintpd
         ws(i, :) = w
      end do

      i = minloc(mins, dim=1)
      mintpd = mins(i)
      w = ws(i, :)

      if(present(all_minima)) all_minima = ws

      call destroy(opt)
   contains
      real(pr) function foo(x, gradient, func_data)
         real(pr), intent(in) :: x(:)
         real(pr), optional, intent(in out) :: gradient(:)
         class(*), optional, intent(in) :: func_data
         foo = tpd(model, z, x, P, T, d=di, dtpd=gradient)
      end function foo
   end subroutine min_tpd
end module yaeos__phase_equilibria_stability
