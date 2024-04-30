program flash_tv
   use yaeos, only: pr, EquilibriaState, flash, PengRobinson76, ArModel, fugacity_tp
   use yaeos__math_linalg, only: solve_system
   use fixtures_models, only: binary_PR76

   implicit none

   integer, parameter :: nc = 2

   real(pr) :: z(nc) = [0.4, 0.6]
   real(pr) :: x(nc) = [0.32424471950363210, 0.67575531029866709]
   real(pr) :: y(nc) = [0.91683466155334536, 8.3165368249135715E-002]
   real(pr) :: vx = 8.4918883298198036E-002
   real(pr) :: vy = 0.32922132295944545

   real(pr) :: F(nc+3), dF(nc+3, nc+3)
   real(pr) :: Fdx(nc+3),dFdx(nc+3, nc+3)
   real(pr) :: dxold(nc+3)

   real(pr) :: K(nc), T=200, V
   real(pr) :: beta=0.12

   real(pr) :: vars(nc+3), dx(nc+3)

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: flash_result

   integer :: its

   model = binary_PR76()

   V = beta*vy + (1-beta) * vx

   K = y/x
   dfdX=0
   df = 0

   flash_result = flash(model, z, T, p_spec=5._pr, iters=its)

   x = flash_result%x
   y = flash_result%y
   y(1) = y(1) + 0.2
   y = y/sum(y)
   beta = flash_result%beta
   vx = flash_result%vx
   vy = flash_result%vy
   K = y/x
   V = beta *vy + (1-beta)*vx

   print *, "true Results"
   print *, flash_result%vx, flash_result%vy, flash_result%beta, beta * vy + (1-beta)*vx
   print *, flash_result%x
   print *, flash_result%y
   print *, "------------"

   call flash_TV_DF(model, z, T, V, K, beta, vx, vy, F, dF)
   ! print *, "numdiff"
   ! do its = 1,nc+3
   !    print *, df(its, :)
   ! end do

   call flash_TV_F(model, z, T, V, K, beta, vx, vy, F, dF)
   ! print *, "anadiff"
   ! do its = 1,nc+3
   !    print *, df(its, :)
   ! end do

   print *, "solve?"
   do its=1,1000
      call flash_TV_DF(model, z, T, V, K, beta, vx, vy, F, dF)

      dxold = dx
      dx = solve_system(dF, -F)

      if (any(isnan(dx)) .or. any(isnan(F))) then
         print *, "reducing"
         K = K - dxold(:nc)
         beta = beta - dxold(nc+1)
         vx = vx - dxold(nc+2)
         vy = vy - dxold(nc+3)

         dx = dxold/2
      end if

      ! print *, its, F

      do while (beta + dx(nc+1) > 1 .or. beta + dx(nc+1) < 0)
         dx = dx/5
      end do

      K = K + dx(:nc)
      beta = beta + dx(nc+1)
      vx = vx + dx(nc+2)
      vy = vy + dx(nc+3)

      if (maxval(abs(F)) < 1e-3) then
         print *, "conv"
         exit
      end if
   end do

   x = z / (beta * (K - 1) + 1)
   y = K * x
   print *, vx, vy, beta, beta*vy + (1-beta)*vx
   print *, x
   print *, y


contains
   subroutine flash_TV_F(model, z, T, V, K, beta, vx, vy, F, dF)
      use yaeos, only: fugacity_vt
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: T
      real(pr), intent(in) :: V
      real(pr), intent(in) :: K(:)
      real(pr), intent(in) :: beta
      real(pr), intent(in) :: vx
      real(pr), intent(in) :: vy
      real(pr), intent(out) :: F(:)
      real(pr), optional, intent(out) :: dF(:, :)

      real(pr) :: x(nc), Px, dPxdn(nc), lnphip_x(nc), dlnPhidx(nc, nc)
      real(pr) :: y(nc), Py, dPydn(nc), lnphip_y(nc), dlnPhidy(nc, nc)

      real(pr) :: dxdK(nc), dydK(nc)

      integer :: i, j

      F = 0

      x = z / (beta * (K - 1) + 1)
      y = K * x

      dxdK = - beta * z / (beta*(K - 1)+1)**2
      dydK  = x + K * dxdK

      call fugacity_vt(model, x, vx, T, P=Px, dPdN=dPxdn, lnphip=lnphip_x, dlnphidn=dlnPhidx)
      call fugacity_vt(model, y, vy, T, P=Py, dPdN=dPydn, lnphip=lnphip_y, dlnphidn=dlnPhidy)

      F(:nc)    = log(K) - lnphip_x + lnphip_y
      F(nc + 1) = Px - Py
      F(nc + 2) = V - (beta*Vy + (1-beta)*Vx)
      F(nc + 3) = sum((z * (1 - K)) /(1 + beta*(K - 1)))

      if (present(df)) then
         df = 0
         do i=1,nc
            do j=1,nc
               df(i, j) = &
                  - dxdK(i) * (dlnPhidx(i, j) + dPxdn(j)/Px) &
                  + dydK(i) * (dlnphidy(i, j) + dPydn(j)/Py)

            end do
            df(i, i) = df(i, i) + 1/K(i)
         end do

         df(nc+1, :nc) = (dPxdn * dxdK - dPydn * dydK)
         df(nc+2, :nc) = 0
         df(nc+3, :nc) = beta * z * (K-1)/(beta*(K-1)+1)**2 - z/(beta*(K-1)+1)
      end if
   end subroutine flash_TV_F

   subroutine flash_TV_DF(model, z, T, V, K, beta, vx, vy, F, dF)
      class(ArModel), intent(in) :: model
      real(pr), intent(in) :: z(:)
      real(pr), intent(in) :: T
      real(pr), intent(in) :: V
      real(pr), intent(in) :: K(:)
      real(pr), intent(in) :: beta
      real(pr), intent(in) :: vx
      real(pr), intent(in) :: vy
      real(pr), intent(out) :: F(:)
      real(pr), intent(out) :: dF(:, :)

      real(pr) :: dx
      integer :: i, nc
      real(pr) :: Kdx(size(z)), betadx, vxdx, vydx, Fdx(size(F))

      nc = size(z)

      dx = 1e-10

      f = 0
      df = 0
      call flash_TV_F(model, z, T, V, K, beta, vx, vy, F, dF)

      do i=1, nc
         Kdx = K
         Kdx(i) = K(i) + dx
         call flash_TV_F(model, z, T, V, Kdx, beta, vx, vy, Fdx)
         df(:, i) = (Fdx - F)/dx
      end do

      betadx = beta + dx
      call flash_TV_F(model, z, T, V, K, betadx, vx, vy, Fdx)
      df(:, nc+1) = (Fdx - F)/dx

      vxdx = vx + dx
      call flash_TV_F(model, z, T, V, K, beta, vxdx, vy, Fdx)
      df(:, nc+2) = (Fdx - F)/dx

      vydx = vy + dx
      call flash_TV_F(model, z, T, V, K, beta, vydx, vy, Fdx)
      df(:, nc+3) = (Fdx - F)/dx
   end subroutine flash_TV_DF
end program flash_tv