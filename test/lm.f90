program test_lm
   !! Test the Levenberg-Marquardt algorithm
   use testing_aux, only: test_title, assert
   use yaeos
   use yaeos__math, only: levenberg_marquardt

   integer, parameter :: nc = 4
   real(pr) :: Tc(nc), Pc(nc), w(nc), kij(nc,nc)
   real(pr) :: K0(nc), beta0, z(nc), P, T, F(nc+1), X(nc+1)
   type(CubicEoS) :: model
   integer :: info

   write(*, *) test_title("Levenberg-Marquardt")

   Tc=[190.578, 425.17, 617.65, 304.20]
   Pc=[46.04, 37.96, 21.07, 73.84]
   w= [0.0104, 0.201, 0.49, 0.225]

   kij(1,:)=[0., 0.027, 0.042, 0.1]
   kij(2,:)=[0.027, 0., 0.008, 0.1257]
   kij(3,:)=[0.042, 0.008, 0., 0.0942]
   kij(4,:)=[0.1, 0.1257, 0.0942, 0.]

   model = PengRobinson76(Tc, Pc, w, kij)

   P = 165.5
   T = 344

   z = [0.69611886, 0.4128533, 0.59425154, 0.41980931]
   K0 = [100.0, 0.01, 0.01, 100.0]
   X = [log(K0), 0.1_pr]

   call levenberg_marquardt(flash_pt_fun, 1e-5_pr, X, F, info)
   call assert(all(abs(F) < 1e-5_pr), "Should converge")

contains
   subroutine flash_pt_fun(m, n, xvar, fvec, iflag)
      integer, intent(in) :: m, n
      real(pr), intent(in) :: xvar(n)
      real(pr), intent(out) :: fvec(m)
      integer, intent(in out) :: iflag

      real(pr) :: y(nc), x(nc), lnphi_x(nc), lnphi_y(nc), beta, K(nc)
      real(pr) :: isofug(nc), rr, denom(nc)

      K = exp(xvar(:nc))
      beta = xvar(nc + 1)

      y = z * K / (1 + beta*(K - 1))
      x = y/K

      call model%lnphi_pt(n=x, P=P, T=T, lnphi=lnphi_x, root_type="liquid")
      call model%lnphi_pt(n=y, P=P, T=T, lnphi=lnphi_y, root_type="vapor")

      isofug = log(x) + lnphi_x - (log(y) + lnphi_y)

      denom = 1 + beta*(K - 1)
      rr = sum(z*(K - 1)/denom)

      Fvec(:nc) = isofug
      Fvec(nc + 1) = rr
   end subroutine flash_pt_fun
end program test_lm
