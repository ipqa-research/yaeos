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

      real(pr) :: X(7)
      real(pr) :: dFdS(7), F(7), dF(7, 7)
      integer :: ns
      real(pr) ::  dS, S, dXdS(7)

      real(pr) :: T, P, Vx, Vy, Vw
      real(pr) :: x1, y1, w1

      real(pr), allocatable :: Ts(:), Ps(:)

      integer :: points
      integer :: iters

      dFdS = 0
      dFdS(7) = -1

      X = log([&
         CEP%x(1)+1e-9, &
         cep%x(1)-1e-9, &
         cep%y(1), &
         cep%Vx+1e-9, &
         cep%Vx-1e-9, &
         cep%Vy, &
         cep%T &
         ])
      T = HUGE(1._pr)
      ns = 0
      if (ns == 0) then
         S = exp(X(2)) - exp(X(1))
      else if (ns == -1) then
         S = X(4) - X(5)
      else
         S = X(ns)
      end if

      dS = 0.01

      points = 0
      F = 0
      P = 10
      allocate(llv%T(0), llv%P(0), llv%x1(0), llv%y1(0), llv%w1(0), llv%Vx(0), llv%Vy(0), llv%Vw(0))
      do while(maxval(abs(F)) < 1e-9 .and. P > 0)
         points = points + 1
         call three_phase_line_F_solve(model, X, ns, S, F, dF, iters)

         x1 = exp(X(1))
         y1 = exp(X(2))
         w1 = exp(X(3))
         Vx = exp(X(4))
         Vy = exp(X(5))
         Vw = exp(X(6))
         T = exp(X(7))
         call model%pressure([x1, 1-x1], Vx, T, P)
         
         if (maxval(abs(F)) < 1e-9 .and. P > 0) then
            llv%T = [llv%T, T]
            llv%P = [llv%P, P]
            llv%x1 = [llv%x1, x1]
            llv%y1 = [llv%y1, y1]
            llv%w1 = [llv%w1, w1]
            llv%Vx = [llv%Vx, Vx]
            llv%Vy = [llv%Vy, Vy]
            llv%Vw = [llv%Vw, Vw]
         end if


         ns = -1
         dS = -0.01 * 3. / real(iters, pr)

         dXdS = solve_system(dF, -dFdS)

         if (T < llv%T(1) - 10) then
            ns = 7
            dS = -0.01 * 3. / real(iters, pr)
            dXdS = dXdS / dXdS(ns)
         end if

         ! Avoid large steps in volume, since they can lead to indeterminations
         ! Probably a better way could be avoiding reaching to the co-volume.
         ! But equations like GERG do not have a co-volume value.
         do while(&
              abs(Vx - exp(X(4) + dS*dXdS(4))) > 0.5*Vx&
              .or. abs(Vy - exp(X(5) + dS*dXdS(5))) > 0.5*Vy&
              .or. abs(Vw - exp(X(6) + dS*dXdS(6))) > 0.5*Vw&
            )
            dS = dS/2
         end do

         X = X + dXdS * dS

         if (ns == 0) then
            S = exp(X(2)) - exp(X(1))
         else if (ns == -1) then
            S = X(4) - X(5)
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

      x1 = exp(XVars(1))
      y1 = exp(XVars(2))
      w1 = exp(XVars(3))
      vx = exp(Xvars(4))
      vy = exp(Xvars(5))
      vw = exp(Xvars(6))
      T =  exp(Xvars(7))

      x = [x1, 1-x1]
      y = [y1, 1-y1]
      w = [w1, 1-w1]

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

      if (ns == 0) then
         F(7) = y1 - x1 - S
      else if (ns == -1) then
         F(7) = log(vx / vy) - S
      else
         F(7) = Xvars(ns) - S
      end if

      df = 0

      df(1, 1) = dlnfxdn(1, 1) - dlnfxdn(1, 2)
      df(2, 1) = dlnfxdn(2, 1) - dlnfxdn(2, 2)
      df(3:4, 1) = 0
      df(5, 1) = dPxdN(1) - dPxdN(2)
      df(6, 1) = 0

      ! Derivatives wrt y1
      df(1:2, 2) = 0
      df(3, 2) = dlnfydn(1, 1) - dlnfydn(1, 2)
      df(4, 2) = dlnfydn(2, 1) - dlnfydn(2, 2)
      df(5, 2) = 0
      df(6, 2) = dPydN(1) - dPydN(2)

      ! Derivatives wrt w1
      df(1, 3) =  - (dlnfwdn(1, 1) - dlnfwdn(1, 2))
      df(2, 3) =  - (dlnfwdn(2, 1) - dlnfwdn(2, 2))
      df(3, 3) =  - (dlnfwdn(1, 1) - dlnfwdn(1, 2))
      df(4, 3) =  - (dlnfwdn(2, 1) - dlnfwdn(2, 2))
      df(5, 3) = - (dPwdN(1) - dPwdN(2))
      df(6, 3) = - (dPwdN(1) - dPwdN(2))

      ! Derivatives wrt vx
      df(1, 4) = dlnfxdv(1)
      df(2, 4) = dlnfxdv(2)
      df(3:4, 4) = 0
      df(5, 4) = dPxdV
      df(6, 4) = 0

      ! Derivatives wrt vy
      df(1:2, 5) = 0
      df(3, 5) = dlnfydv(1)
      df(4, 5) = dlnfydv(2)
      df(5, 5) = 0
      df(6, 5) = dPydV

      ! Derivatives wrt vw
      df(1, 6) = -dlnfwdv(1)
      df(2, 6) = -dlnfwdv(2)
      df(3, 6) = -dlnfwdv(1)
      df(4, 6) = -dlnfwdv(2)
      df(5, 6) = -dPwdV
      df(6, 6) = -dPwdV

      ! Derivatives wrt T
      df(1:2, 7) = dlnfxdT - dlnfwdT
      df(3:4, 7) = dlnfydT - dlnfwdT
      df(5, 7) = dPxdT - dPwdT
      df(6, 7) = dPydT - dPwdT

      if (ns == 0) then
         df(7, 1) = -1
         df(7, 2) = 1
      else if (ns == -1) then
         df(7, 4) = 1/vx
         df(7, 5) = -1/vy
      else
         df(7, ns) = 1
      end if

      ! Multiply by variable values due to log-transform
      df(:, 1) = df(:, 1) * x1
      df(:, 2) = df(:, 2) * y1
      df(:, 3) = df(:, 3) * w1
      df(:, 4) = df(:, 4) * vx
      df(:, 5) = df(:, 5) * vy
      df(:, 6) = df(:, 6) * vw
      df(:, 7) = df(:, 7) * T
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
      max_tries = 50

      dX = 10
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

         X = Xold + dX

         if (present(iters)) iters = i

         if (any(isnan(F))) error stop
      end do
   end subroutine three_phase_line_F_solve
end module yaeos__equilibria_binaries
