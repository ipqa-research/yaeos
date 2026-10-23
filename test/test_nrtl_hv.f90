program main
   use testing_aux, only: test_title, assert
   use yaeos, only: pr, R
   use yaeos__models_ge_base, only: nrtl_hv_ge, nrtl_hv_tdep
   use yaeos__models_ge_nrtlhv, only: NRTLHV
   use hyperdual_mod
   implicit none

   integer, parameter :: nc = 3
   real(pr) :: n(nc), T, Ge, Gen(nc), GenT(nc), dgendn_num(nc), Gen2(nc, nc), dgendn2(nc, nc), GeT, GeT2
   real(pr) :: alpha(nc, nc), b(nc), gij(nc, nc), tau(nc, nc), dtaudt(nc, nc), dtaudt2(nc, nc)
   real(pr) :: adiff_ge, adiff_get, adiff_get2, adiff_gen(nc), adiff_gent(nc), adiff_gen2(nc, nc)
   real(pr) :: kij(nc, nc)

   type(hyperdual) :: t_hd, n_hd(nc), ge_hdv
   type(NRTLHV) :: model

   integer :: i, j

   write(*, *) test_title("NRTL-HV model")

   gij = reshape([ &
      0.0, 0.1, 0.2, &
      0.3, 0.0, 0.4, &
      0.5, 0.6, 0.0], [nc, nc]) * 100 * R
   alpha = reshape([ &
      0.0, 1.1, 1.2, &
      1.3, 0.0, 1.5, &
      1.6, 1.7, 0.0], [nc, nc])

   b = [0.1, 0.2, 0.3]

   model = NRTLHV(gji0=gij, gjiT=0*gij, alpha=alpha, b=b)

   n = [0.1, 0.3, 0.5]
   T = 350.0

   n_hd = n
   t_hd = t
   t_hd%f1 = 1._pr
   t_hd%f2 = 1._pr
   ge_hdv = ge_hd(n_hd, t_hd)
   
   adiff_ge = ge_hdv%f0
   adiff_get = ge_hdv%f1
   adiff_get2 = ge_hdv%f12

   do i=1,nc
      n_hd = n
      t_hd = t
      n_hd(i)%f1 = 1._pr
      t_hd%f2 = 1._pr
      ge_hdv = ge_hd(n_hd, t_hd)
      adiff_gen(i) = ge_hdv%f1
      adiff_gent(i) = ge_hdv%f12
   end do

   do i=1,nc
      do j=1,nc
         n_hd = n
         t_hd = t
         n_hd(i)%f1 = 1._pr
         n_hd(j)%f2 = 1._pr
         ge_hdv = ge_hd(n_hd, t_hd)
         adiff_gen2(i, j) = ge_hdv%f12
      end do
   end do

   call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GenT, Gen2)
   call assert(abs(Ge - adiff_ge) < 1e-6, "Excess Gibbs energy")
   call assert(abs(GeT - adiff_get) < 1e-6, "Excess Gibbs energy derivative wrt T")
   call assert(abs(GeT2 - adiff_get2) < 1e-6, "Excess Gibbs energy 2nd derivative wrt T")
   call assert(all(abs(Gen - adiff_gen) < 1e-6), "Excess Gibbs energy derivative wrt n")
   call assert(all(abs(GenT - adiff_gent) < 1e-6), "Excess Gibbs energy derivative wrt T and n")
   call assert(all(abs(Gen2 - adiff_gen2) < 1e-6), "Excess Gibbs energy second derivative wrt n")
   
   call model%excess_gibbs(n, T, Ge)
   call assert(abs(Ge - adiff_ge) < 1e-6, "Excess Gibbs energy alone: Individual call")
   
   call model%excess_gibbs(n, T, GeT=GeT)
   call assert(abs(GeT - adiff_get) < 1e-6, "Excess Gibbs energy derivative wrt T: Individual call")
   
   call model%excess_gibbs(n, T, GeT2=GeT2)
   call assert(abs(GeT2 - adiff_get2) < 1e-6, "Excess Gibbs energy 2nd derivative wrt T: Individual call")
   
   call model%excess_gibbs(n, T, Gen=Gen)
   call assert(all(abs(Gen - adiff_gen) < 1e-6), "Excess Gibbs energy derivative wrt n: Individual call")
   call model%excess_gibbs(n, T, GeTn=GenT)
   call assert(all(abs(GenT - adiff_gent) < 1e-6), "Excess Gibbs energy derivative wrt T and n: Individual call")
   call model%excess_gibbs(n, T, Gen2=Gen2)
   call assert(all(abs(Gen2 - adiff_gen2) < 1e-6), "Excess Gibbs energy second derivative wrt n: Individual call")

contains
   type(hyperdual) function ge_hd(n, T)
      type(hyperdual), intent(in) :: n(:), T

      type(hyperdual) :: E(nc, nc), U, D, xi(nc), theta(nc), xxi(nc, nc)
      type(hyperdual) :: tau(nc, nc)

      real(pr) :: tin

      integer :: i, j
      tin = t%f0

      ! call tdep(Tin, gij, tau, dtaudt)
      tau = gij/(R*T)
      E = exp(-alpha * tau)

      do i=1,nc
         xxi(:, i)     = E(:, i) * b * tau(:, i) * n
         xi(i)     = sum(E(:, i) * b * tau(:, i) * n)
         theta(i)  = sum(E(:, i) * b * n)
      end do

      ge_hd = 0._pr
      do i=1,nc
         U = 0._pr
         D = 0._pr
         do j=1,nc
            U = U + n(j) * b(j) * E(j, i) * tau(j, i)
            D = D + n(j) * b(j) * E(j, i)
         end do
         ge_hd = ge_hd + n(i) * U/D
      end do

      ge_hd = ge_hd * R * T
   end function ge_hd
end program main