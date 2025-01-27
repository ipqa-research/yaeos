module yaeos__m_s_sp
   !! Module to calculate saturation points
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel, size
   implicit none

contains

   subroutine saturation_F(model, z, X, ns, S, F, dF, dPdVz, dPdVy)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      real(pr), intent(out) :: F(:)
      real(pr), optional, intent(out) :: dF(:, :)
      real(pr), intent(out) :: dPdVz, dPdVy

      ! Variables
      real(pr) :: T, Vz, Vy
      real(pr) :: z(size(model))

      ! Main phase variables
      real(pr) :: lnfug_z(size(model)), dlnfug_dn_z(size(model), size(model))
      real(pr) :: dlnfug_dT_z(size(model)), dlnfug_dV_z(size(model))
      real(pr) :: dlnfug_dP_z(size(model))
      real(pr) :: Pz, dPdTz, dPdn_z(size(z))

      ! incipient phase variables
      real(pr) :: y(size(z))
      real(pr) :: lnfug_y(size(model)), dlnfug_dn_y(size(model), size(model))
      real(pr) :: dlnfug_dT_y(size(model)), dlnfug_dV_y(size(model))
      real(pr) :: dlnfug_dP_y(size(model))
      real(pr) :: Py, dPdTy, dPdn_y(size(z))

      real(pr) :: lnPspec

      integer :: j, nc

      nc = size(z)

      y  = z * exp(X(:nc))
      Vz = exp(X(nc+1))
      Vy = exp(X(nc+2))
      T  = exp(X(nc+3))
      lnPspec = X(nc+4)

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
      F(nc + 3) = lnPspec - log(Py)
      F(nc + 4) = X(ns) - S

      if (present(dF)) then
         dF = 0

         ! isofugacity
         do j=1,nc
            df(:nc, j) = y(j) * dlnfug_dn_y(:, j)
         end do

         dF(:nc, nc+1) = -dlnfug_dV_z * Vz
         dF(:nc, nc+2) =  dlnfug_dV_y * Vy
         dF(:nc, nc+3) =  T * (dlnfug_dT_y - dlnfug_dT_z)

         ! mass balance
         df(nc+1, :nc) = y

         ! pressure equality
         df(nc+2, :nc)  = y * dPdn_y
         df(nc+2, nc+1) = -dPdVz * Vz
         df(nc+2, nc+2) = dPdVy * Vy
         df(nc+2, nc+3) = T*(dPdTy - dPdTz)

         df(nc+3, :nc) = -y * dPdn_y/Py
         df(nc+3, nc+1) = 0
         df(nc+3, nc+2) = -dPdVy * Vy / Py
         df(nc+3, nc+3) = -dPdTy * T / Py
         df(nc+3, nc+4) = 1

         df(nc+4, ns) = 1
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

      integer :: nc

      real(pr) :: F(size(X))
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: dx(size(X))

      nc = size(X) - 2

      its = 0
      dX = 1
      do while (its < max_iterations)
         its = its + 1
         call saturation_TP(model=model, z=z, kind=kind, X=X, ns=ns, S=S, F=F, dF=dF, dFdS=dFdS)

         dX = solve_system(dF, -F)
         
         do while (abs(dX(nc+1)) > 0.1)
            dX(nc+1) = dX(nc+1) / 2
         end do
         
         do while (abs(dX(nc+2)) > 0.1)
            dX(nc+2) = dX(nc+2) / 2
         end do

         X = X + dX
         if (all(abs(F) < tol)) exit
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

      real(pr) :: dPdVz, dPdVy
      
      integer :: nc
      real(pr) :: F(size(X))
      real(pr) :: dF(size(X), size(X))
      real(pr) :: dFdS(size(X))
      real(pr) :: Xold(size(X)), dx(size(X)), dx_old(size(x))

      nc = size(X) - 4

      its = 0
      do while (its < max_iterations)
         call saturation_F(model, z, X, ns, S, F, dF, dPdVz, dPdVy)
         if (all(abs(F) < tol)) exit

         dX = solve_system(dF, -F)

         X = X + dX

         its = its + 1
      end do
   end subroutine solve_VxVyT
end module yaeos__m_s_sp
