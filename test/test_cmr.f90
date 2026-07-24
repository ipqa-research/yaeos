module tests_cmr
   use yaeos
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert
   integer, parameter :: nc = 3
contains
   subroutine test_aijk_numdiff
      use yaeos__models_ar_cubic_mixing_base, only: CMR_aijk
      use yaeos__extra_fluids, only: co2_h2o_isop, FixtureFluid
      real(pr) :: k(nc, nc, nc), dkdt(nc, nc, nc), dkdt2(nc, nc, nc)
      type(FixtureFluid) :: fluid

      integer :: i, j, l
      real(pr) :: T, dt=1e-3
      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
      real(pr) :: df, df2

      fluid = co2_h2o_isop()

      k(1, 1, 2) = 0.1
      k(1, 2, 1) = 0.1
      k(2, 1, 1) = 0.1

      k(1, 2, 3) = 0.3
      k(2, 1, 3) = 0.3
      k(3, 2, 1) = 0.3

      k(2, 2, 1) = 0.2
      k(2, 1, 2) = 0.2
      k(1, 2, 2) = 0.2
      dkdt = 0
      dkdt2 = 0

      T = 200

      associate(model => fluid%ar_model)
         select type(model)
          type is (CubicEoS)
            call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
            call CMR_aijk(&
               ai, daidt, daidt2, k, dkdt, dkdt2, &
               a, dadt, dadt2 &
               )
         end select
      end associate

      do i=1,nc
         do j=1,nc
            do l=1,nc
               df = (f(T+dT, i, j, l) - f(T-dT, i, j, l))/(2*dt)

               df2 = (f(T+dT, i, j, l) - 2*f(T, i, j, l) + f(T-dT, i, j, l))/(dT**2)
               call assert(abs((df - dadt(i, j, l))/df) < 1e-5, "First order analytic and numerical derivatives should be similar")
               call assert(abs((df2 - dadt2(i, j, l))/df) < 1e-5, "Second order analytic and numerical derivatives should be similar")
            end do
         end do
      end do

   contains
      real(pr) function f(T, i, j, l)
         real(pr), intent(in) :: T
         integer, intent(in) :: i, j, l
         real(pr) :: ai(nc), daidt(nc), daidt2(nc)
         real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
         real(pr) :: eps=1e-5_pr

         associate(model => fluid%ar_model)
            select type(model)
             type is (CubicEoS)
               call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
               call CMR_aijk(&
                  ai, daidt, daidt2, k, dkdt, dkdt2, &
                  a, dadt, dadt2 &
                  )
            end select
         end associate
         f = a(i, j, l)
      end function f
   end subroutine test_aijk_numdiff

   subroutine test_aijk_kexpt_numdiff
      use yaeos__models, only: CMRTD
      use yaeos__models_ar_cubic_cubic_mixing, only: kijk_exp_tdep
      use yaeos__extra_fluids, only: co2_h2o_isop, FixtureFluid
      real(pr) :: k(nc, nc, nc)
      real(pr) :: k0(nc, nc, nc)
      real(pr) :: tref(nc, nc, nc)
      type(FixtureFluid) :: fluid
      type(CMRTD) :: mixrule

      integer :: i, j, l
      real(pr) :: n(nc), T, dt=1e-3
      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
      real(pr) :: df, df2

      fluid = co2_h2o_isop()

      k = 0
      k0 = 0
      tref = 0

      k(1, 1, 2) = 0.1
      k(1, 2, 1) = 0.1
      k(2, 1, 1) = 0.1

      k(1, 2, 3) = 0.3
      k(2, 1, 3) = 0.3
      k(3, 2, 1) = 0.3

      k(2, 2, 1) = 0.2
      k(2, 1, 2) = 0.2
      k(1, 2, 2) = 0.2

      k0(1, 1, 2) = 0.01
      k0(1, 2, 1) = 0.01
      k0(2, 1, 1) = 0.01

      k0(1, 2, 3) = 0.03
      k0(2, 1, 3) = 0.03
      k0(3, 2, 1) = 0.03

      k0(2, 2, 1) = 0.02
      k0(2, 1, 2) = 0.02
      k0(1, 2, 2) = 0.02

      tref(1, 1, 2) = 100
      tref(1, 2, 1) = 100
      tref(2, 1, 1) = 100

      tref(1, 2, 3) = 300
      tref(2, 1, 3) = 300
      tref(3, 2, 1) = 300

      tref(2, 2, 1) = 200
      tref(2, 1, 2) = 200
      tref(1, 2, 2) = 200

      mixrule = CMRTD(k=k, k0=k0, tref=tref, l=0*k0)

      T = 200

      associate(model => fluid%ar_model)
         select type(model)
          type is (CubicEoS)
            call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
            call kijk_exp_tdep(mixrule, &
               T, ai, daidt, daidt2, &
               a, dadt, dadt2 &
               )
         end select
      end associate

      do i=1,nc
         do j=1,nc
            do l=1,nc
               df = (f(T+dT, i, j, l) - f(T-dT, i, j, l))/(2*dt)
               df2 = (f(T+dT, i, j, l) - 2*f(T, i, j, l) + f(T-dT, i, j, l))/(dT**2)
               call assert(abs((df - dadt(i, j, l))/df) < 1e-5, "First order analytic and numerical derivatives should be similar")
               call assert(abs((df2 - dadt2(i, j, l))/df) < 1e-5, "Second order analytic and numerical derivatives should be similar")
            end do
         end do
      end do

   contains
      real(pr) function f(T, i, j, l)
         real(pr), intent(in) :: T
         integer, intent(in) :: i, j, l
         real(pr) :: eps=1e-5_pr
         real(pr) :: ai(nc), daidt(nc), daidt2(nc)
         real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)

         associate(model => fluid%ar_model)
            select type(model)
             type is (CubicEoS)
               call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
               call kijk_exp_tdep(mixrule, &
                  T, ai, daidt, daidt2, &
                  a, dadt, dadt2 &
                  )
            end select
         end associate

         f = a(i, j, l)
      end function f
   end subroutine test_aijk_kexpt_numdiff

   subroutine test_D_kexpt_numdiff
      use yaeos__models, only: CMRTD, CMR
      use yaeos__models_ar_cubic_cubic_mixing, only: kijk_exp_tdep
      use yaeos__extra_fluids, only: co2_h2o_isop, FixtureFluid
      real(pr) :: k(nc, nc, nc)
      real(pr) :: k0(nc, nc, nc)
      real(pr) :: tref(nc, nc, nc)
      type(FixtureFluid) :: fluid
      type(CMRTD) :: mixrule

      integer :: i, j, l
      real(pr) :: n(nc), V=0, T, dt=1e-2, dni=1e-9, dn(nc)
      real(pr) :: ai(nc), daidt(nc), daidt2(nc)
      real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)
      real(pr) :: df, df2
      real(pr) :: f1, f2, f3, f4, f0

      real(pr) :: D, D_num
      real(pr) :: dDdV
      real(pr) :: dDdT
      real(pr) :: dDdT2
      real(pr) :: dDdV2
      real(pr) :: dDdTV
      real(pr) :: dDi(nc), dDi_num(nc)
      real(pr) :: dDidV(nc)
      real(pr) :: dDidT(nc), dDidT_num(nc)
      real(pr) :: dDij(nc, nc), dDij_num(nc, nc)

      fluid = co2_h2o_isop()

      k = 0
      k0 = 0
      tref = 0

      k(1, 1, 2) = 0.1
      k(1, 2, 1) = 0.1
      k(2, 1, 1) = 0.1

      k(1, 2, 3) = 0.3
      k(2, 1, 3) = 0.3
      k(3, 2, 1) = 0.3

      k(2, 2, 1) = 0.2
      k(2, 1, 2) = 0.2
      k(1, 2, 2) = 0.2

      k0(1, 1, 2) = 0.01
      k0(1, 2, 1) = 0.01
      k0(2, 1, 1) = 0.01

      k0(1, 2, 3) = 0.03
      k0(2, 1, 3) = 0.03
      k0(3, 2, 1) = 0.03

      k0(2, 2, 1) = 0.02
      k0(2, 1, 2) = 0.02
      k0(1, 2, 2) = 0.02

      tref(1, 1, 2) = 100
      tref(1, 2, 1) = 100
      tref(2, 1, 1) = 100

      tref(1, 2, 3) = 300
      tref(2, 1, 3) = 300
      tref(3, 2, 1) = 300

      tref(2, 2, 1) = 200
      tref(2, 1, 2) = 200
      tref(1, 2, 2) = 200

      n = fluid%z0
      n = [0.9, 0.5, 0.2]

      mixrule = CMRTD(k=k, k0=k0, tref=tref, l=0*k0)

      T = 200

      print *, "AAAAAAAAAAAA"
      associate(model => fluid%ar_model)
         select type(model)
          type is (CubicEoS)
            call model%set_mixrule(mixrule)
            call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
            call model%mixrule%Dmix(&
               n, V, T, &
               ai, daidt, daidt2, &
               D=D, &
               dDdV=dDdV, &
               dDdT=dDdT, &
               dDdT2=dDdT2, &
               dDdV2=dDdV2, &
               dDdTV=dDdTV, &
               dDi=dDi, &
               dDidV=dDidV, &
               dDidT=dDidT, &
               dDij=dDij &
               )
         end select
      end associate
      print *, "AAAAAAAAAAAAfin"

      call autodiff(n, T, D_num, df, df2, dDi_num, dDidT_num, dDij_num)

      call assert(abs(D - D_num) < 1e-10, "D")
      call assert(abs(dDdT - df) < 1e-8, "dDdT")
      call assert(abs(dDdT2 - df2) < 1e-8, "d2DdT2")
      call assert(allclose(dDi, dDi_num, 1e-8_pr), "dDi")

      call assert(allclose(dDidT, dDidT_num, 1e-8_pr), "dDidT")
      call assert(allclose([dDij], [dDij_num], 1e-8_pr), "dDij")
   contains
      real(pr) function f(n, T)
         real(pr), intent(in) :: n(nc)
         real(pr), intent(in) :: T
         real(pr) :: ai(nc), daidt(nc), daidt2(nc)
         real(pr) :: a(nc, nc, nc), dadt(nc, nc, nc), dadt2(nc, nc, nc)

         real(pr) :: D !! Mixture attractive parameter \(n^2a_{mix}\)
         real(pr) :: dDdV !! \(\frac{dD}{dT}\)
         real(pr) :: dDdT !! \(\frac{dD}{dV}\)
         real(pr) :: dDdT2 !! \(\frac{d^2D}{dT^2}\)
         real(pr) :: dDdV2 !! \(\frac{d^2D}{dV^2}\)
         real(pr) :: dDdTV !! \(\frac{d^2D}{dTV\)
         real(pr) :: dDi(nc) !! \(\frac{dD}{dn_i}\)
         real(pr) :: dDidV(nc) !! \(\frac{d^2D}{dVn_i}\)
         real(pr) :: dDidT(nc) !! \(\frac{d^2D}{dTn_i}\)
         real(pr) :: dDij(nc, nc)!! \(\frac{d^2D}{dn_{ij}}\)

         associate(model => fluid%ar_model)
            select type(model)
             type is (CubicEoS)
               call model%calculate_attractive_parameters(T, ai, daidt, daidt2)
               call model%mixrule%Dmix(&
                  n, V, T, &
                  ai, daidt, daidt2, &
                  D=f, &
                  dDdV=dDdV, &
                  dDdT=dDdT, &
                  dDdT2=dDdT2, &
                  dDdV2=dDdV2, &
                  dDdTV=dDdTV, &
                  dDi=dDi, &
                  dDidV=dDidV, &
                  dDidT=dDidT, &
                  dDij=dDij &
                  )
            end select
         end associate
      end function f

      function hd_d(n, T, ai, daidt, daidt2)
         use hyperdual_mod
         type(hyperdual), intent(in) :: n(:)
         type(hyperdual), intent(in) :: T
         real(pr), intent(in) :: ai(:), daidt(:), daidt2(:)
         type(hyperdual) :: hd_D

         type(hyperdual) :: hd_a(size(n)), hd_k, c

         integer :: i, j, l, nc

         hd_D = 0._pr
         nc = size(n)

         do i=1,nc
            do j=1,nc
               do l=1,nc

                  if (Tref(i, j, l) /= 0) then
                     c = k0(i, j, l) * exp(-T/Tref(i, j, l))
                     hd_k = k(i, j, l) + c
                  else
                     hd_k = 0._pr
                  end if

                  hd_a = ai

                  if (T%f1 == 1 .and. T%f2 == 1) then
                     hd_a%f1 = daidt
                     hd_a%f2 = daidt
                     hd_a%f12 = daidt2
                  else if(T%f1 == 1) then
                     hd_a%f1 = daidt
                  else if(T%f2 == 1) then
                     hd_a%f2 = daidt
                  end if

                  hd_D = hd_d + (n(i) * n(j) * n(l)) * (hd_a(i) * hd_a(j) * hd_a(l))**(1._pr/3._pr) * (1._pr - hd_k)
               end do
            end do
         end do

         hd_D = hd_D / sum(n)
      end function hd_d

      subroutine autodiff(n, T, D, dDdT, dDdT2, dDi, dDiT, dDij)
         use hyperdual_mod
         real(pr), intent(in) :: n(:), T
         real(pr), intent(out) :: D, dDdT, dDdT2
         real(pr), intent(out) :: dDi(:), dDiT(:), dDij(:, :)

         type(hyperdual) :: hd_d_i, hd_t, hd_n(size(n))
         integer :: i, nc

         nc = size(n)

         hd_n = n
         hd_t = t

         hd_t%f1 = 1
         hd_t%f2 = 1
         hd_d_i = hd_d(hd_n, hd_t, ai, daidt, daidt2)

         D = hd_d_i%f0

         dDdT = hd_d_i%f1
         dDdT2 = hd_d_i%f12


         do i=1,nc

            hd_n%f1 = 0
            hd_n%f2 = 0
            hd_n%f12 = 0
            hd_t%f1 = 0
            hd_t%f2 = 0
            hd_t%f12 = 0


            hd_n = n
            hd_t = t
            hd_n(i)%f1 = 1
            hd_t%f2 = 1
            hd_d_i = hd_d(hd_n, hd_t, ai, daidt, daidt2)
            dDi(i) = hd_d_i%f1
            dDidT(i) = hd_d_i%f12
         end do

      end subroutine autodiff
   end subroutine test_D_kexpt_numdiff

   subroutine test_compare_with_qmr
      use yaeos, only: QMR, CMR
      use yaeos__extra_fluids, only: co2_h2o_isop, FixtureFluid
      type(FixtureFluid) :: fluid

      real(pr) :: k(nc, nc, nc), l(nc, nc, nc)
      real(pr) :: n(3), V, T
      real(pr) :: D_qmr
      real(pr) :: dDdV_qmr, dDdT_qmr, dDdV2_qmr, dDdT2_qmr, dDi_qmr, dDdTV_qmr
      real(pr) :: dDidV_qmr(3), dDidT_qmr(3), dDij_qmr(3,3)

      real(pr) :: dDdV_cmr, dDdT_cmr, dDdV2_cmr, dDdT2_cmr, dDi_cmr, dDdTV_cmr
      real(pr) :: dDidV_cmr(3), dDidT_cmr(3), dDij_cmr(3,3)

      fluid = co2_h2o_isop()

      T = 200
      V = 2
      n = fluid%z0

      k = 0

      associate(model => fluid%ar_model)
         select type(model)
          type is (CubicEoS)
            associate (mr => model%mixrule)
               select type (mr)
                type is (QMR)

                  k(1, 1, 2) = mr%k(1, 2)
                  k(1, 2, 1) = k(1, 1, 2)
                  k(1, 2, 1) = k(1, 1, 2)

                  k(1, 1, 3) = mr%k(1, 3)
                  k(1, 3, 1) = k(1, 1, 3)
                  k(1, 3, 1) = k(1, 1, 3)

                  k(2, 2, 3) = mr%k(2, 3)
                  k(2, 3, 2) = k(2, 2, 3)
                  k(2, 3, 2) = k(2, 2, 3)

                  l(1, 1, 2) = mr%l(1, 2)
                  l(1, 2, 1) = l(1, 1, 2)
                  l(1, 2, 1) = l(1, 1, 2)

                  l(1, 1, 3) = mr%l(1, 3)
                  l(1, 3, 1) = l(1, 1, 3)
                  l(1, 3, 1) = l(1, 1, 3)

                  l(2, 2, 3) = mr%l(2, 3)
                  l(2, 3, 2) = l(2, 2, 3)
                  l(2, 3, 2) = l(2, 2, 3)

               end select
            end associate
         end select
      end associate

   end subroutine test_compare_with_qmr
end module tests_cmr

program test_cmr
   use tests_cmr, only: test_aijk_numdiff
   use tests_cmr, only: test_aijk_kexpt_numdiff
   use tests_cmr, only: test_D_kexpt_numdiff

   use testing_aux, only: test_title

   print *, test_title("Cubic Mixing Rules")
   call test_aijk_numdiff
   call test_aijk_kexpt_numdiff
   call test_D_kexpt_numdiff
end program test_cmr
