module yaeos__m_s_sp
   !! Module to calculate saturation points
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel, size
   implicit none

contains

   subroutine saturation_F(model, z, X, ns, S, F, dF)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(out) :: F(:)
      real(pr), optional, intent(out) :: dF(:, :)

      ! Variables
      real(pr) :: T, Vz, Vy
      real(pr) :: z(size(model))

      ! Main phase variables
      real(pr) :: lnfug_z(size(model)), dlnfug_dn_z(size(model), size(model))
      real(pr) :: dlnfug_dT_z(size(model)), dlnfug_dV_z(size(model))
      real(pr) :: dlnfug_dP_z(size(model))
      real(pr) :: Pz, dPdTz, dPdVz, dPdn_z(size(z))

      ! incipient phase variables
      real(pr) :: y(size(z))
      real(pr) :: lnfug_y(size(model)), dlnfug_dn_y(size(model), size(model))
      real(pr) :: dlnfug_dT_y(size(model)), dlnfug_dV_y(size(model))
      real(pr) :: dlnfug_dP_y(size(model))
      real(pr) :: Py, dPdTy, dPdVy, dPdn_y(size(z))

      integer :: j, nc

      nc = size(z)

      y  = z * exp(X(:nc))
      Vz = exp(X(nc+1))
      Vy = exp(X(nc+2))
      T  = X(nc+3)

      if (present(df)) then
         call model%lnfug_vt(&
            n=z, V=Vz, T=T, P=Pz, dPdT=dPdTz, dPdV=dPdVz, dPdn=dPdn_z, &
            lnf=lnfug_z, &
            dlnfdV=dlnfug_dV_z, dlnfdT=dlnfug_dT_z, dlnfdn=dlnfug_dn_z &
            )

         call model%lnfug_vt(&
            n=y, V=Vy, T=T, P=Py, dPdT=dPdTy, dPdV=dPdVy, dPdn=dPdn_y, &
            lnf=lnfug_y, &
            dlnfdV=dlnfug_dV_y, dlnfdT=dlnfug_dT_y, dlnfdn=dlnfug_dn_y &
            )
      else
         call model%lnfug_vt(n=z, V=Vz, T=T, P=Pz, lnf=lnfug_z)
         call model%lnfug_vt(n=y, V=Vy, T=T, P=Py, lnf=lnfug_y)
      end if

      F = 0

      F(:nc) = lnfug_y - lnfug_z
      F(nc + 1) = sum(y - z)
      F(nc + 2) = Py - Pz
      F(nc + 3) = X(ns) - S

      if (present(dF)) then
         dF = 0

         ! isofugacity
         do j=1,nc
            df(:nc, j) = y(j) * dlnfug_dn_y(:, j)
         end do

         dF(:nc, nc+1) = -dlnfug_dV_z * Vz
         dF(:nc, nc+2) =  dlnfug_dV_y * Vy
         dF(:nc, nc+3) =  dlnfug_dT_y - dlnfug_dT_z

         ! mass balance
         df(nc+1, :nc) = y

         ! pressure equality
         df(nc+2, :nc)  = y * dPdn_y
         df(nc+2, nc+1) = -dPdVz * Vz
         df(nc+2, nc+2) = dPdVy * Vy
         df(nc+2, nc+3) = dPdTy - dPdTz

         df(nc+3, ns) = 1
      end if
   end subroutine saturation_F

   subroutine saturation_TP(model, kind, z, X, ns, S, F, dF, dFdS)
      class(ArModel), intent(in) :: model
      character(len=*), intent(in) :: kind
      real(pr), intent(in) :: z(size(model))

      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S

      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:, :)
      real(pr), intent(out) :: dFdS(:)

      character(len=14) :: kind_z, kind_y

      real(pr) :: y(size(X)-2)
      real(pr) :: Vz, Vy
      real(pr) :: lnPhi_z(size(X)-2), lnPhi_y(size(X)-2)
      real(pr) :: dlnphi_dt_z(size(X)-2), dlnphi_dt_y(size(X)-2)
      real(pr) :: dlnphi_dp_z(size(X)-2), dlnphi_dp_y(size(X)-2)
      real(pr) :: dlnphi_dn_z(size(X)-2, size(X)-2), dlnphi_dn_y(size(X)-2, size(X)-2)

      real(pr) :: T, P, K(size(X)-2)

      integer :: i, j, nc

      nc = size(X)-2

      F = 0
      dF = 0

      K = exp(X(:nc))
      T = exp(X(nc+1))
      P = exp(X(nc+2))

      y = K*z

      select case(kind)
       case ("bubble")
         kind_z = "liquid"
         kind_y = "vapor"
       case ("dew")
         kind_z = "vapor"
         kind_y = "liquid"
       case ("liquid-liquid")
         kind_z = "liquid"
         kind_y = "liquid"
       case default
         kind_z = "stable"
         kind_y = "stable"
      end select

      call model%lnphi_pt(&
         z, P, T, V=Vz, root_type=kind_z, &
         lnPhi=lnphi_z, dlnPhidt=dlnphi_dt_z, &
         dlnPhidp=dlnphi_dp_z, dlnphidn=dlnphi_dn_z &
         )
      call model%lnphi_pt(&
         y, P, T, V=Vy, root_type=kind_y, &
         lnPhi=lnphi_y, dlnPhidt=dlnphi_dt_y, &
         dlnPhidp=dlnphi_dp_y, dlnphidn=dlnphi_dn_y &
         )

      F(:nc) = X(:nc) + lnPhi_y - lnPhi_z
      F(nc + 1) = sum(y - z)
      F(nc + 2) = X(ns) - S

      ! Jacobian Matrix
      do j=1,nc
         df(:nc, j) = dlnphi_dn_y(:, j) * y(j)
         df(j, j) = dF(j, j) + 1
      end do

      df(:nc, nc + 1) = T * (dlnphi_dt_y - dlnphi_dt_z)
      df(:nc, nc + 2) = P * (dlnphi_dp_y - dlnphi_dp_z)

      df(nc + 1, :nc) = y

      df(nc + 2, :) = 0
      df(nc + 2, ns) = 1

      dFdS = 0
      dFdS(nc+2) = -1
   end subroutine saturation_TP

   subroutine solve_TP(model, kind, z, X, ns, S, tol, max_iterations, its)
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      character(len=*), intent(in) :: kind
      real(pr), intent(in) :: z(:)
      real(pr), intent(inout) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_iterations
      integer, intent(out) :: its

      real(pr) :: F(size(X))
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: dx(size(X))

      its = 0
      do while (its < max_iterations)
         call saturation_TP(model, kind, z, X, ns, S, F, dF, dFdS)
         if (all(abs(F) < tol)) exit

         dX = solve_system(dF, -F)
         X = X + dX

         do while (any(isnan(F)))
            X = X - dx*0.9
            call saturation_TP(model, kind, z, X, ns, S, F, dF, dFdS)
         end do
         its = its + 1
      end do

   end subroutine solve_TP

   subroutine solve_VxVyT(model, z, X, ns, S, tol, max_iterations, its)
      use yaeos__math, only: solve_system
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(inout) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(in) :: tol
      integer, intent(in) :: max_iterations
      integer, intent(out) :: its

      real(pr) :: F(size(X))
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: dx(size(X))

      its = 0
      do while (its < max_iterations)
         call saturation_F(model, z, X, ns, S, F, dF)
         if (all(abs(F) < tol)) exit

         dX = solve_system(dF, -F)

         X = X + dX

         do while (any(isnan(F)))
            X = X - 0.9_pr * dx
            call saturation_F(model, z, X, ns, S, F, dF)
         end do
         its = its + 1
      end do
   end subroutine solve_VxVyT
end module yaeos__m_s_sp
