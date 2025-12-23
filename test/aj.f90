module flasher
   use yaeos

   private

   public :: solve_flash

   class(ArModel), allocatable :: model_
   real(pr) :: z_(4)
   real(pr) :: P_, T_
   integer :: nc_


contains

   subroutine solve_flash(model, nc, z, P, T, k0, beta0, K, beta, info, xvar, fvec)
      use yaeos__math, only: levenberg_marquardt
      use yaeos__equilibria_multiphase_flash, only: solve_mp_flash_point
      use yaeos__equilibria_auxiliar, only: k_wilson
      class(ArModel), intent(in) :: model
      integer, intent(in) :: nc
      real(pr), intent(in) :: z(nc)
      real(pr), intent(in) :: P, T
      real(pr), intent(in) :: K0(nc)
      real(pr), intent(in) :: beta0
      real(pr), intent(out) :: K(nc)
      real(pr), intent(out) :: beta
      real(pr) :: xvar(nc+1), fvec(nc+1)
      integer :: info
      integer :: i

      integer, parameter :: np=1

      real(pr) :: F(nc+np+1+2), X(nc+np+1+2)
      integer :: ns1, ns2, max_iters=1000, beta_0_index, iters
      real(pr) :: S1, S2
      character(len=14), parameter :: kinds_x(np) = "liquid"
      character(len=14), parameter :: kind_w = "vapor"
      logical :: less_phases

      ! ns1 = nc+np+1+1
      ! ns2 = nc+np+1+2
      ! S1 = log(P)
      ! S2 = log(T)

      ! X(:nc) = log(1/k_wilson(model, T, P))
      ! X(nc+np) = beta0
      ! X(nc+np+1) = 1-beta0

      ! call solve_mp_flash_point(model, z, np, kinds_x, kind_w, X, ns1, S1, ns2, S2, max_iters, F, less_phases, beta_0_index, iters)
      ! K = 1/exp(X(:nc))
      ! beta = X(nc+np+1)

      model_ = model
      z_ = z
      P_ = P
      T_ = T
      nc_ = nc

      xvar(:nc) = log(K0)
      xvar(nc+1) = beta0
      call levenberg_marquardt(flash_pt_fun, 1e-15_pr, xvar, fvec, info)

      K = exp(xvar(:nc))
      beta = xvar(nc+1)
   end subroutine solve_flash
   
   subroutine flash_pt_fun(m, n, xvar, fvec, iflag)
      integer, intent(in) :: m, n
      real(pr), intent(in) :: xvar(n)
      real(pr), intent(out) :: fvec(n)
      integer, intent(in out) :: iflag

      real(pr) :: y(nc_), x(nc_), lnphi_x(nc_), lnphi_y(nc_), beta, K(nc_)
      real(pr) :: isofug(nc_), rr, denom(nc_)

      K = exp(xvar(:nc_))
      beta = xvar(nc_ + 1)

      y = z_ * K / (1 + beta*(K - 1))
      x = y/K

      call model_%lnphi_pt(n=x, P=P_, T=T_, lnphi=lnphi_x, root_type="liquid")
      call model_%lnphi_pt(n=y, P=P_, T=T_, lnphi=lnphi_y, root_type="vapor")

      isofug = log(x) + lnphi_x - (log(y) + lnphi_y)

      denom = 1 + beta*(K - 1)
      rr = sum(z_*(K - 1)/denom)

      Fvec(:nc_) = isofug
      Fvec(nc_ + 1) = rr
   end subroutine flash_pt_fun
end module flasher

program test_lm
   !! Test the Levenberg-Marquardt algorithm
   use testing_aux, only: test_title, assert
   use yaeos
   use yaeos__math, only: powel_hybrid, levenberg_marquardt
   use yaeos__equilibria_flash, only: flash_no_beta_limits
   use flasher

   integer, parameter :: nc = 4
   integer, parameter :: ncontacts = 20
   real(pr) :: Tc(nc), Pc(nc), w(nc), kij(nc,nc)
   real(pr) :: K0(nc), beta0, z(nc), P, T, F(nc+1), X(nc+1)
   real(pr) :: xx(nc), yy(nc), beta, K(nc)
   real(pr) :: zi(nc), z0(nc), tl, tl_old, mins_tl(ncontacts)
   real(pr) :: min_tl
   type(CubicEoS) :: model
   type(EquilibriumState) :: fr

   real(pr) :: o(nc), g(nc)
   real(pr) :: zs1(2*ncontacts+2, nc)
   real(pr) :: zs2(2*ncontacts+2, nc)
   integer :: info
   integer :: i, j
   real(pr) :: a

   real(pr) :: last_five(5)

   logical :: plateou_1 = .false.
   logical :: plateou_2 = .false.
   logical :: plateou_3 = .false.

   Tc=[190.578, 425.17, 617.65, 304.20]
   Pc=[46.04, 37.96, 21.07, 73.84]
   w= [0.0104, 0.201, 0.49, 0.225]

   kij(1,:)=[0., 0.027, 0.042, 0.1]
   kij(2,:)=[0.027, 0., 0.008, 0.1257]
   kij(3,:)=[0.042, 0.008, 0., 0.0942]
   kij(4,:)=[0.1, 0.1257, 0.0942, 0.]

   model = PengRobinson76(Tc, Pc, w, kij)

   P = 120
   T = 344

   z = [0.69611886, 0.4128533, 0.59425154, 0.41980931]
   K0 = [100.0, 0.01, 0.01, 100.0]
   beta0 = 0.5
   call solve_flash(model, nc, z, P, T, K0, beta0, K, beta, info, X, F)

   O = [0.2, 0.15, 0.65, 0.]
   G = [0.2, 0., 0., 0.8]

   zs1(1, :) = G
   zs1(2, :) = O
   a = 0.5
   min_tl = 100

   last_five = [50, 30, 15, 62, 3]

   min_tl = sqrt(sum((o-g)**2))
   zs2 = 0

   open(1, file="aj.dat")
   do while(min_tl > 0.01)
      mins_tl = 100
      contacts: do i=1,ncontacts
         zs2(1, :) = G
         tl = 10

         !$OMP PARALLEL DO PRIVATE(j, zi, z0, z, info, beta0, beta, K, yy, xx, tl_old, tl, X, F) &
         !$OMP& default(shared)
         do j=1,2*i-1,2
            zi = zs1(j, :)
            z0 = zs1(j+1, :)

            z = a * zi + (1-a) * z0
            beta0 = a

            ! call solve_flash(model, nc, z, P, T, K0, beta0, K, beta, info, X, F)
            fr = flash_no_beta_limits(model, z, T, beta0=beta0, P_SPEC=P, k0=K0, iters=info)
            K = fr%y/fr%x
            beta = fr%beta
            print *, fr
            yy = z * K / (1 + beta*(K - 1))
            xx = yy/K

            zs2(j+1, :) = xx
            zs2(j+2, :) = yy

            ! print *, i, j
            ! print "(A,x,*(F6.4,x))", "z ", z, beta
            ! print "(A,x,*(F6.4,x))", "x ", xx
            ! print "(A,x,*(F6.4,x))", "y ", yy
            ! print *, ""
            ! print *, "zi", zi
            ! print *, "z0", z0
            ! print *, ""
            ! print "(I1,',',I1,',',*(E14.5,','))", i, j, transpose(zs2)


            tl_old = tl
            tl = sqrt(sum((xx-yy)**2))
            if (i == ncontacts) print *, tl
            if (tl < minval(mins_tl)) mins_tl(j) = tl
         end do
         !$OMP END PARALLEL DO

         min_tl = minval(mins_tl)

         zs2(2*i+2, :) = O
         ! print "(A,x,*(4(F6.4,x),'|'))", "z1", transpose(zs1)
         ! print "(A,x,*(4(F6.4,x),'|'))", "z2", transpose(zs2)
         ! print *, ""
         zs1 = zs2
      end do contacts

      print *, P, i, min_tl
      P = P + 10 * min_tl
      call exit
   end do
   close(1)
end program test_lm
