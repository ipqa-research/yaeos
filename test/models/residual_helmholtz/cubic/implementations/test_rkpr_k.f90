program main
   use testing_aux, only: assert, test_title

   write(*, *) test_title("RKPR replicate acentric factor ")
   call n2
   call binary_n2_c36
contains
   subroutine n2
      use yaeos
      integer, parameter :: nc=1

      type(CubicEoS) :: ar_model

      real(pr) :: tc(nc), pc(nc), w(nc), zc(nc)

      real(pr) :: delta_1(nc)

      type(PurePsat) :: psat
      real(pr) :: Vl, Vv

      integer :: i

      tc = [126.2_pr]
      pc = [34.0_pr]
      w = [0.0377215_pr]
      zc = [0.289_pr]

      do i=10,30
         delta_1 = [i]/10._pr
         ar_model = RKPR(tc, pc, w, zc, delta_1=delta_1)
         call assert(abs(w(1) - (-1 - log10(ar_model%Psat_pure(1, 0.7*Tc(1), Vl, Vv)/Pc(1)))) < 1e-6, "Acentric factor")
      end do
   end subroutine n2

   subroutine binary_n2_c36
      use yaeos
      integer, parameter :: nc=2

      type(CubicEoS) :: ar_model

      real(pr) :: tc(nc), pc(nc), w(nc), zc(nc)

      real(pr) :: delta_1(nc)

      type(QMRTD) :: mixrule
      type(PurePsat) :: psat
      real(pr) :: kij_inf(nc, nc), kij_0(nc, nc), lij(nc, nc)
      real(pr) :: t_ref(nc, nc)

      real(pr) :: Vl, Vv

      integer :: i

      tc = [126.2_pr, 874.0_pr]
      pc = [34.0_pr, 6.8_pr]
      w = [0.0377215_pr, 1.52596_pr]
      zc = [0.289_pr, 0.196_pr]

      delta_1 = [3.0_pr, 3.045_pr]

      ar_model = RKPR(tc, pc, w, zc, delta_1=delta_1)

      kij_0(1, :) = [0.0_pr, -0.17720104893220182_pr]
      kij_0(2, :) = [-0.17720104893220182_pr, 0.0_pr]

      kij_inf(1, :) = [0.0_pr, 0.1515709974620788_pr]
      kij_inf(2, :) = [0.1515709974620788_pr, 0.0_pr]

      lij(1, :) = [0_pr, 0_pr]
      lij(2, :) = [0_pr, 0_pr]

      t_ref(1, :) = [0.0_pr, 126.2_pr]
      t_ref(2, :) = [126.2_pr, 0.0_pr]

      mixrule = QMRTD(k=kij_inf, k0=kij_0, Tref=t_ref, l=lij)

      call ar_model%set_mixrule(mixrule)
         
      call assert(abs(w(1) - (-1 - log10(ar_model%Psat_pure(1, 0.7*Tc(1), Vl, Vv)/Pc(1)))) < 1e-6, "Acentric factor")
      call assert(abs(w(2) - (-1 - log10(ar_model%Psat_pure(2, 0.7*Tc(2), Vl, Vv)/Pc(2)))) < 1e-6, "Acentric factor")

   end subroutine binary_n2_c36

   real(pr) function get_k(model, i)
      use yaeos
      type(CubicEoS) :: model
      integer :: i
      associate(a => model%alpha)
         select type(a)
          type is (AlphaRKPR)
            get_k = a%k(i)
         end select
      end associate
   end function get_k


end program main
