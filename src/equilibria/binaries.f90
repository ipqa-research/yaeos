module yaeos__equilibria_binaries
   !! Module with routines particular to binary mixtures.
   use yaeos__constants, only: pr
   use yaeos__models, only: ArModel
   use yaeos__equilibria_equilibrium_state, only: EquilibriumState
   use yaeos__math, only: solve_system

   private
   public :: find_llcl, three_phase_line_F_solve, three_phase_line_F
   public :: BinaryThreePhase, binary_llv_from_cep

   type :: BinaryThreePhase
      !! # `BinaryThreePhase`
      !! Structure to hold the results of a binary LLV line calculation.
      !! # Description
      !! This structure holds the results of a binary LLV line calculation.
      !! Pressures are calculated using the composition of the liquid phase `x`.
      real(pr), allocatable :: x1(:) !! Mole fraction of component 1 in liquid phase x
      real(pr), allocatable :: y1(:) !! Mole fraction of component 1 in liquid phase y
      real(pr), allocatable :: w1(:) !! Mole fraction of component 1 in vapor phase w
      real(pr), allocatable :: Vx(:) !! Volume of liquid phase x [L/mol]
      real(pr), allocatable :: Vy(:) !! Volume of liquid phase y [L/mol]
      real(pr), allocatable :: Vw(:) !! Volume of vapor phase w [L/mol]
      real(pr), allocatable :: T(:) !! Temperature [K]
      real(pr), allocatable :: P(:) !! Pressure [bar]
   end type BinaryThreePhase

   type :: LLVVarEnum
      integer :: x1=1
      integer :: x2=2
      integer :: y1=3
      integer :: y2=4
      integer :: w1=5
      integer :: w2=6
      integer :: Vx=7
      integer :: Vy=8
      integer :: Vw=9
      integer :: T=10
   end type LLVVarEnum

   type(LLVVarEnum), parameter :: llv_vars = LLVVarEnum()

   integer, parameter :: n_llv_vars = 10

contains

   subroutine find_llcl(model, z0, zi, P, a, V, T)
      !! # `find_llcl`
      !! Find an initial guess for the critical L-L line of a binary mixture.
      !!
      !! # Description
      !!
      !!
      !! # Examples
      !!
      !!
      !! # References
      !! [1] M. Cismondi, M.L. Michelsen, Global phase equilibrium
      !! calculations:
      !! Critical lines, critical end points and liquid–liquid–vapour
      !! equilibrium in binary mixtures, The Journal of Supercritical Fluids 39
      !!  (2007) 287–295. https://doi.org/10.1016/j.supflu.2006.03.011.
      implicit none
      class(ArModel), intent(in) :: model !! Thermodynamic model to use
      real(kind=pr), intent(in) :: P !! Pressure [bar]
      real(kind=pr), intent(in) :: z0(2) !! Mole fractions of original fluid
      real(kind=pr), intent(in) :: zi(2) !! Mole fractions of new fluid
      real(kind=pr), intent(out) :: a !! Mole fraction of new fluid
      real(kind=pr), intent(out) :: V !! Volume [L/mol]
      real(kind=pr), intent(in out) :: T !! Temperature [K]

      real(kind=pr) :: z(2), lnphi(2), dlnphidn(2, 2)

      integer :: i
      integer :: tries
      real(kind=pr) :: as(50)
      real(kind=pr) :: lambdas(50)

      do tries = 1, 2
         do i = 1, 50
            a = real(i, pr) / 51.
            z = a * zi + (1 - a) * z0
            call model%lnphi_pt(&
               z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=dlnphidn, root_type=&
               "liquid")
            lambdas(i) = 1 - dlnphidn(1, 2)

            as(i) = a
         end do

         i = minloc(lambdas, dim=1)
         a = as(i)
         z = a * zi + (1 - a) * z0

         if (lambdas(i) > 0) then
            do while (lambdas(i) > 0)
               T = T - 1
               call model%lnphi_pt(z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=&
                  dlnphidn, root_type="liquid")
               lambdas(i) = 1 - dlnphidn(1, 2)
            end do
         else
            do while (lambdas(i) < 0)
               T = T + 1
               call model%lnphi_pt(z, P=P, T=T, V=V, lnPhi=lnPhi, dlnphidn=&
                  dlnphidn, root_type="liquid")
               lambdas(i) = 1 - dlnphidn(1, 2)
            end do
         end if
      end do
   end subroutine find_llcl

   type(BinaryThreePhase) function binary_llv_from_cep(model, cep) result(llv)
      !! # `binary_llv_from_cep`
      !! Calculate the LLV line from a converged critical end point (CEP).
      !!
      !! # Description
      !! From a converged critical end point (CEP) of a binary mixture, this
      !! function calculates the three-phase line (LLV) by solving the
      !! corresponding system of equations (defined at [[three_phase_line_F]])
      !! at each point.
      !! To trace the whole line a continuation method is used to obtain good
      !! initial guesses for each point.
      !! The specification used to trace the line is initially the difference
      !! between the mole fractions of the two liquid phases, and then it is
      !! switched to temperature.
      class(ArModel), intent(in) :: model !! Thermodynamic model to use
      type(EquilibriumState), intent(in) :: cep !! Converged critical end point.

      integer, parameter :: nvars = 10

      real(pr) :: X(nvars)
      real(pr) :: dFdS(nvars), F(nvars), dF(nvars, nvars)
      integer :: ns
      real(pr) ::  dS, S, dXdS(nvars)

      real(pr) :: T, P, Vx, Vy, Vw
      real(pr) :: xx(2), y(2), w(2)

      real(pr), allocatable :: Ts(:), Ps(:)

      integer :: points
      integer :: iters

      real(pr) :: delta

      dFdS = 0
      dFdS(10) = -1

      delta = 0.0001 * cep%x(1)
      xx(1) = cep%x(1) - delta
      y(1) = cep%x(1) + delta
      w(1) = cep%y(1)

      xx(2) = 1 - xx(1)

      y(2) = 1 - y(1)
      w(2) = 1 - w(1)

      T = cep%T
      P = cep%P

      call model%volume(xx, P, T, Vx, root_type="liquid")
      call model%volume(y, P, T, Vy, root_type="liquid")
      Vw = cep%Vy

      X = log([&
         xx, &
         y, &
         w, &
         Vx, &
         Vy, &
         Vw, &
         T &
         ])
      T = HUGE(1._pr)
      ns = 0
      ! if (ns == 0) then
      !    S = exp(X(2)) - exp(X(1))
      ! else if (ns == -1) then
      !    S = X(4) - X(5)
      ! else
      !    S = X(ns)
      ! end if
      S = xx(1) - y(1) - delta
      dS = -0.0001

      points = 0
      F = 0
      P = 10
      allocate(&
         llv%T(0), llv%P(0), llv%x1(0), llv%y1(0), llv%w1(0), &
         llv%Vx(0), llv%Vy(0), llv%Vw(0))
      do while((T > 100 .or. P > 1e-8_pr) .and. maxval(abs(F)) < 1e-9)
         points = points + 1
         call three_phase_line_F_solve(model, X, ns, S, F, dF, iters)
         if (iters >= 100) then
            X = X - 0.9 * dXdS * dS 
            S = S - 0.9 * dS
            call three_phase_line_F_solve(model, X, ns, S, F, dF, iters)
         end if

         xx = exp(X(llv_vars%x1:llv_vars%x2))
         y = exp(X(llv_vars%y1:llv_vars%y2))
         w = exp(X(llv_vars%w1:llv_vars%w2))
         Vx = exp(X(llv_vars%Vx))
         Vy = exp(X(llv_vars%Vy))
         Vw = exp(X(llv_vars%Vw))
         T = exp(X(llv_vars%T))
         call model%pressure(xx, Vx, T, P)

         if (maxval(abs(F)) < 1e-9 .and. P > 0) then
            llv%T = [llv%T, T]
            llv%P = [llv%P, P]
            llv%x1 = [llv%x1, xx(1)]
            llv%y1 = [llv%y1, y(1)]
            llv%w1 = [llv%w1, w(1)]
            llv%Vx = [llv%Vx, Vx]
            llv%Vy = [llv%Vy, Vy]
            llv%Vw = [llv%Vw, Vw]
         end if


         dXdS = solve_system(dF, -dFdS)

         ns = llv_vars%y1
         dS = dXdS(ns) * dS
         dXdS = dXdS / dXdS(ns)

         ! if (T < llv%T(1) - 10) then
         !    ns = 7
         !    dS = dXdS(ns) * dS
         !    dS = -0.001 * 3. / real(iters, pr)
         !    ! dS = sign(max(abs(dS), 0.001_pr), dS)
         !    dXdS = dXdS / dXdS(ns)
         ! end if

         do while(&
            abs(Vx - exp(X(llv_vars%Vx) + dS*dXdS(llv_vars%Vx))) > 0.5*Vx&
            .or. abs(Vy - exp(X(llv_vars%Vy) + dS*dXdS(llv_vars%Vy))) > 0.5*Vy&
            .or. abs(Vw - exp(X(llv_vars%Vw) + dS*dXdS(llv_vars%Vw))) > 0.5*Vw&
            )
            dS = dS/2
         end do

         X = X + dXdS * dS

         if (ns == 0) then
            S = exp(X(llv_vars%x1)) - exp(X(llv_vars%y1))
         else if (ns == -1) then
            S = X(llv_vars%Vx) - X(llv_vars%Vy)
         else
            S = X(ns)
         end if
      end do
   end function binary_llv_from_cep

   subroutine three_phase_line_F(model, Xvars, ns, S, F, dF)
      !! # `three_phase_line_F`
      !!
      !! # Description
      !! Calculate the function vector and Jacobian for the three-phase
      !! line (LLV) of a binary mixture. Phases are defined as `x`, `y` and `w`.
      !! Which are two liquid phases and one vapor phase respectively.
      !!
      !! The system of equations is defined as:
      !! \[
      !! \begin{align*}
      !! & \ln \phi_i^{(x)} - \ln \phi_i^{(w)} = 0 \\
      !! & \ln \phi_i^{(y)} - \ln \phi_i^{(w)} = 0 \\
      !! & P^{(x)} - P^{(w)} = 0 \\
      !! & P^{(y)} - P^{(w)} = 0 \\
      !! & \text{specification}
      !! \begin{cases}
      !! y_1 - x_1 - S = 0 & \text{if } ns = 0 \\
      !! \ln v_x - \ln v_y - S = 0 & \text{if } ns = -1 \\
      !! X_{ns} - S = 0 & \text{otherwise}
      !! \end{cases}
      !! \end{align*}
      !! \]
      !!
      !! # References
      !! [1] M. Cismondi, M.L. Michelsen, Global phase equilibrium
      !! calculations:
      !! Critical lines, critical end points and liquid–liquid–vapour
      !! equilibrium in binary mixtures, The Journal of Supercritical Fluids 39
      !!  (2007) 287–295. https://doi.org/10.1016/j.supflu.2006.03.011.
      use yaeos__math, only: derivative_dxk_dni
      class(ArModel), intent(in) :: model !! Thermodynamic model to use
      real(kind=pr), intent(in) :: Xvars(:) !! Input vector
      integer, intent(in) :: ns !! Specified variable index
      real(pr), intent(in) :: S !! Specified variable value
      real(kind=pr), intent(out) :: F(:) !! Function vector
      real(kind=pr), intent(out) :: dF(:, :) !! Jacobian

      ! Variables describing the three phases
      real(pr) :: x1, y1, w1
      real(pr) :: vx, vy, vw, T
      real(pr) :: x(2), y(2), w(2)

      real(pr) :: Px, Py, Pw
      real(pr) :: isofug_1(2), isofug_2(2)
      real(pr) :: Peq_1, Peq_2

      real(pr) :: lnf_x(2), lnf_y(2), lnf_w(2)
      real(pr) :: dlnfxdn(2, 2), dlnfydn(2, 2), dlnfwdn(2, 2)
      real(pr) :: dlnfxdv(2), dlnfydv(2), dlnfwdv(2)
      real(pr) :: dlnfxdP(2), dlnfydP(2), dlnfwdP(2)
      real(pr) :: dlnfxdT(2), dlnfydT(2), dlnfwdT(2)

      real(pr) :: dPxdN(2), dPydN(2), dPwdN(2)
      real(pr) :: dPxdV, dPydV, dPwdV
      real(pr) :: dPxdT, dPydT, dPwdT
      real(pr) :: dxdn(2, 2)

      x = exp(XVars(llv_vars%x1:llv_vars%x2))
      y = exp(XVars(llv_vars%y1:llv_vars%y2))
      w = exp(XVars(llv_vars%w1:llv_vars%w2))
      vx = exp(Xvars(llv_vars%Vx))
      vy = exp(Xvars(llv_vars%Vy))
      vw = exp(Xvars(llv_vars%Vw))
      T =  exp(Xvars(llv_vars%T))
      x1 = x(1)
      y1 = y(1)

      ! x = [x1, 1-x1]
      ! y = [y1, 1-y1]
      ! w = [w1, 1-w1]

      call model%lnfug_vt(&
         x, Vx, T, Px, lnf_x, &
         dlnfxdV, dlnfxdT, dlnfxdn, &
         dPxdV, dPxdT, dPxdn &
         )
      call model%lnfug_vt(&
         y, Vy, T, Py, lnf_y, &
         dlnfydV, dlnfydT, dlnfydn, &
         dPydV, dPydT, dPydn &
         )
      call model%lnfug_vt(&
         w, Vw, T, Pw, lnf_w, &
         dlnfwdV, dlnfwdT, dlnfwdn, &
         dPwdV, dPwdT, dPwdn &
         )


      ! Calculate isofugacity coefficients
      isofug_1 = lnf_x - lnf_w
      isofug_2 = lnf_y - lnf_w
      Peq_1 = Px - Pw
      Peq_2 = Py - Pw

      F(1:2) = isofug_1
      F(3:4) = isofug_2
      F(5) = Peq_1
      F(6) = Peq_2
      F(7) = sum(x) - 1
      F(8) = sum(y) - 1
      F(9) = sum(w) - 1

      if (ns == 0) then
         F(10) = y1 - x1 - S
      else if (ns == -1) then
         F(10) = log(vx / vy) - S
      else
         F(10) = Xvars(ns) - S
      end if

      df = 0

      ! Derivatives wrt x1
      df(1, llv_vars%x1) = dlnfxdn(1, 1)
      df(2, llv_vars%x1) = dlnfxdn(2, 1)
      df(3:4, llv_vars%x1) = 0
      df(5, llv_vars%x1) = dPxdN(1)
      df(6, llv_vars%x1) = 0
      df(7, llv_vars%x1) = 1

      ! Derivatives wrt x2
      df(1, llv_vars%x2) = dlnfxdn(1, 2)
      df(2, llv_vars%x2) = dlnfxdn(2, 2)
      df(3:4, llv_vars%x2) = 0
      df(5, llv_vars%x2) = dPxdN(2)
      df(6, llv_vars%x2) = 0
      df(7, llv_vars%x2) = 1

      ! Derivatives wrt y1
      df(1:2, llv_vars%y1) = 0
      df(3, llv_vars%y1) = dlnfydn(1, 1)
      df(4, llv_vars%y1) = dlnfydn(2, 1)
      df(5, llv_vars%y1) = 0
      df(6, llv_vars%y1) = dPydN(1)
      df(8, llv_vars%y1) = 1

      ! Derivatives wrt y2
      df(1:2, llv_vars%y2) = 0
      df(3, llv_vars%y2) = dlnfydn(1, 2)
      df(4, llv_vars%y2) = dlnfydn(2, 2)
      df(5, llv_vars%y2) = 0
      df(6, llv_vars%y2) = dPydN(2)
      df(8, llv_vars%y2) = 1

      ! Derivatives wrt w1
      df(1, llv_vars%w1) =  - dlnfwdn(1, 1)
      df(2, llv_vars%w1) =  - dlnfwdn(2, 1)
      df(3, llv_vars%w1) =  - dlnfwdn(1, 1)
      df(4, llv_vars%w1) =  - dlnfwdn(2, 1)
      df(5, llv_vars%w1) = - dPwdN(1)
      df(6, llv_vars%w1) = - dPwdN(1)
      df(9, llv_vars%w1) = 1

      ! Derivatives wrt w2
      df(1, llv_vars%w2) =  - dlnfwdn(1, 2)
      df(2, llv_vars%w2) =  - dlnfwdn(2, 2)
      df(3, llv_vars%w2) =  - dlnfwdn(1, 2)
      df(4, llv_vars%w2) =  - dlnfwdn(2, 2)
      df(5, llv_vars%w2) = - dPwdN(2)
      df(6, llv_vars%w2) = - dPwdN(2)
      df(9, llv_vars%w2) = 1

      ! Derivatives wrt Vx
      df(1, llv_vars%Vx) = dlnfxdv(1)
      df(2, llv_vars%Vx) = dlnfxdv(2)
      df(3:4, llv_vars%Vx) = 0
      df(5, llv_vars%Vx) = dPxdV
      df(6, llv_vars%Vx) = 0

      ! Derivatives wrt Vy
      df(1:2, llv_vars%Vy) = 0
      df(3, llv_vars%Vy) = dlnfydv(1)
      df(4, llv_vars%Vy) = dlnfydv(2)
      df(5, llv_vars%Vy) = 0
      df(6, llv_vars%Vy) = dPydV

      ! Derivatives wrt vw
      df(1, llv_vars%Vw) = -dlnfwdv(1)
      df(2, llv_vars%Vw) = -dlnfwdv(2)
      df(3, llv_vars%Vw) = -dlnfwdv(1)
      df(4, llv_vars%Vw) = -dlnfwdv(2)
      df(5, llv_vars%Vw) = -dPwdV
      df(6, llv_vars%Vw) = -dPwdV

      ! Derivatives wrt T
      df(1:2, llv_vars%T) = dlnfxdT - dlnfwdT
      df(3:4, llv_vars%T) = dlnfydT - dlnfwdT
      df(5, llv_vars%T) = dPxdT - dPwdT
      df(6, llv_vars%T) = dPydT - dPwdT

      if (ns == 0) then
         df(10, llv_vars%x1) = -1
         df(10, llv_vars%y1) = 1
      else if (ns == -1) then
         df(10, llv_vars%Vx) = 1/vx
         df(10, llv_vars%Vy) = -1/vy
      else
         df(10, ns) = 1
      end if

      ! Multiply by variable values due to log-transform
      df(:, llv_vars%x1) = df(:, llv_vars%x1) * x(1)
      df(:, llv_vars%x2) = df(:, llv_vars%x2) * x(2)
      df(:, llv_vars%y1) = df(:, llv_vars%y1) * y(1)
      df(:, llv_vars%y2) = df(:, llv_vars%y2) * y(2)
      df(:, llv_vars%w1) = df(:, llv_vars%w1) * w(1)
      df(:, llv_vars%w2) = df(:, llv_vars%w2) * w(2)
      df(:, llv_vars%Vx) = df(:, llv_vars%Vx) * vx
      df(:, llv_vars%Vy) = df(:, llv_vars%Vy) * vy
      df(:, llv_vars%Vw) = df(:, llv_vars%Vw) * vw
      df(:, llv_vars%T) = df(:, llv_vars%T) * T
   end subroutine three_phase_line_F

   subroutine three_phase_line_F_solve(model, X, ns, S, F, dF, iters)
      !! # `three_phase_line_F_solve`
      !!
      !! # Description
      !! Solve the system of equations defined in `three_phase_line_F` using
      !! a Newton-Raphson method.
      !! It will make a maximum of 50 iterations to converge. With a tolerance
      !! of 1e-9 in the maximum absolute value of the function vector or the
      !! step vector.
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model !! Thermodynamic model to use
      real(kind=pr), intent(in out) :: X(:) !! Input/output vector
      integer, intent(in) :: ns !! Specified variable index
      real(pr), intent(in) :: S !! Specified variable value
      real(kind=pr), intent(out) :: F(:) !! Function vector
      real(kind=pr), intent(out) :: dF(:, :) !! Jacobian
      integer, optional, intent(out) :: iters !! Number of iterations performed

      integer :: i
      integer :: max_tries
      real(kind=pr) :: tol
      real(kind=pr) :: res_norm
      real(kind=pr) :: dX(size(X))
      real(kind=pr) :: Xold(size(X))

      tol = 1e-9_pr
      max_tries = 100

      dX = 10
      F = 10
      do i = 1, max_tries
         call three_phase_line_F(model, X, ns, S, F, dF)
         res_norm = maxval(abs(F))
         if (res_norm < tol .or. maxval(abs(dX)) < 1e-9) exit

         Xold = X

         ! Solve the linear system dF * dX = -F
         dX = solve_system(dF, -F)

         do while (maxval(abs(dX)) > 0.1)
            dX = dX / 2
         end do

         ! if (exp(X(1) + dX(1)) < 0 .or. X(1) + dX(1) > 0) then
         !    dX(1) = log(1e-15_pr) - X(1)
         ! end if
         ! if (exp(X(2) + dX(2)) < 0 .or. X(2) + dX(2) > 0) then
         !    dX(2) = log(1e-15_pr) - X(2)
         ! end if

         X = Xold + dX


         if (present(iters)) iters = i
      end do
   end subroutine three_phase_line_F_solve
end module yaeos__equilibria_binaries
