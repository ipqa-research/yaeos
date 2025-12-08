program main
   use yaeos, only: pr
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title
   implicit none

   call test_dortmund_cons_mix()
   call test_dortmund_cons_pure()
   call test_against_caleb_thermo()

contains
   subroutine test_dortmund_cons_mix()
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_dortmund, UNIFAC
      use yaeos__consistency_gemodel, only: ge_consistency
      use yaeos__consistency_gemodel, only: numeric_ge_derivatives
      use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

      type(UNIFAC) :: model

      integer, parameter :: nc = 4, ng = 4

      type(Groups) :: molecules(nc)

      type(GeGCModelParameters) :: parameters

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_n, Gen_n(nc), GeT_n, GeT2_n, GeTn_n(nc), Gen2_n(nc, nc)

      real(pr) :: n(nc), T
      real(pr) :: dt, dn

      real(pr) :: eq58, eq59(size(n)), eq60(size(n),size(n)), eq61(size(n))

      integer :: i, j

      print *, test_title("Dortmund model mixture consistency test")

      T = 303.15
      n = [400.0, 100.0, 300.0, 200.0]

      dt = 0.1_pr
      dn = 0.01_pr

      molecules(1)%groups_ids = [16]
      molecules(1)%number_of_groups = [1]

      molecules(2)%groups_ids = [1, 2]
      molecules(2)%number_of_groups = [2, 4]

      molecules(3)%groups_ids = [1, 2, 14]
      molecules(3)%number_of_groups = [1, 1, 1]

      molecules(4)%groups_ids = [1, 7, 8, 78, 79]
      molecules(4)%number_of_groups = [2, 1, 1, 3, 1]

      parameters = DortmundParameters()

      model = setup_dortmund(molecules, parameters)

      ! ========================================================================
      ! Call analytic derivatives
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

      ! ========================================================================
      ! Call numeric derivatives
      ! ------------------------------------------------------------------------
      call numeric_ge_derivatives(model, n, T, dn, 0.01_pr, Ge=Ge_n, GeT=GeT_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen=Gen_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeT2=GeT2_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, GeTn=GeTn_n)
      call numeric_ge_derivatives(model, n, T, dn, dt, Ge=Ge_n, Gen2=Gen2_n)

      ! Derivatives checks
      call assert(abs(Ge - Ge_n) < 1e-10, "Ge numeric")
      call assert(abs(GeT - GeT_n) < 1e-4, "GeT numeric")
      call assert(allclose(Gen, Gen_n, 1e-4_pr), "Gen numeric")
      call assert(abs(GeT2 - GeT2_n) < 1e-4_pr, "GeT2 numeric")
      call assert(allclose(GeTn, GeTn_n, 1e-3_pr), "GeTn numeric")
      call assert(allclose(Gen2(1,:), Gen2_n(1,:), 1e-1_pr), "Gen2 numeric")
      call assert(allclose(Gen2(2,:), Gen2_n(2,:), 1e-1_pr), "Gen2 numeric")
      call assert(allclose(Gen2(3,:), Gen2_n(3,:), 1e-1_pr), "Gen2 numeric")

      ! ========================================================================
      ! Consistency tests
      ! ------------------------------------------------------------------------
      call ge_consistency(model, n, t, eq58=eq58)
      call ge_consistency(model, n, t, eq59=eq59)
      call ge_consistency(model, n, t, eq60=eq60)
      call ge_consistency(model, n, t, eq61=eq61)

      ! Eq 58
      call assert(abs(eq58) < 1e-10_pr, "Consistency Eq 58")

      ! Eq 59
      do i=1,size(n)
         call assert(abs(eq59(i)) < 1e-10_pr, "Consistency Eq 59")
      end do

      ! Eq 60
      do i=1,size(n)
         do j=1,size(n)
            call assert(abs(eq60(i, j)) < 1e-10_pr, "Consistency Eq 60")
         end do
      end do

      ! Eq 61
      do i=1,size(n)
         call assert(abs(eq61(i)) < 1e-10_pr, "Consistency Eq 61")
      end do
   end subroutine test_dortmund_cons_mix

   subroutine test_dortmund_cons_pure()
      use yaeos, only: pr
      use yaeos, only: Groups, UNIFAC, setup_dortmund

      type(UNIFAC) :: model
      type(Groups) :: molecules(1)

      real(pr) :: Ge, Gen(1), GeT, GeT2, GeTn(1), Gen2(1, 1), ln_gammas(1)
      real(pr) :: T, n(1)

      integer :: i, j

      print *, test_title("Dortmund model pure consistency test")

      T = 303.15
      n = [400.0]

      ! Hexane [CH3, CH2]
      molecules(1)%groups_ids = [1, 2]
      molecules(1)%number_of_groups = [2, 4]

      model  = setup_dortmund(molecules)

      ! Evaluate Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! All must be zero for a pure compounds
      call assert(abs(Ge) < 1e-10_pr, "Ge = 0")
      call assert(abs(GeT) < 1e-10_pr, "GeT = 0")
      call assert(abs(GeT2) < 1e-10_pr, "GeT2 = 0")

      do i=1,size(n)
         call assert(abs(Gen(i)) < 1e-10_pr, "Gen = 0")
         call assert(abs(GeTn(i)) < 1e-10_pr, "GeTn = 0")
         call assert(abs(ln_gammas(i)) < 1e-10_pr, "ln_gammas = 0")
      end do

      do i=1,size(n)
         do j=1,size(n)
            call assert(abs(Gen2(i, j)) < 1e-10_pr, "Gen2 = 0")
         end do
      end do
   end subroutine test_dortmund_cons_pure

   subroutine test_against_caleb_thermo()
      ! https://github.com/CalebBell/thermo
      use yaeos, only: pr, R
      use yaeos, only: Groups, setup_dortmund, UNIFAC

      type(UNIFAC) :: model

      integer, parameter :: nc = 3, ng = 4

      type(Groups) :: molecules(nc)

      real(pr) :: Ge, Gen(nc), GeT, GeT2, GeTn(nc), Gen2(nc, nc)
      real(pr) :: Ge_i, Gen_i(nc), GeT_i, GeT2_i, GeTn_i(nc), Gen2_i(nc, nc)
      real(pr) :: ln_gammas(nc)

      real(pr) :: n(nc), T, n_t

      print *, test_title("Dortmund test against Caleb Bell's thermo")

      T = 303.15_pr
      n = [2.0_pr, 5.0_pr, 3.0_pr]
      n_t = sum(n)

      molecules(1)%groups_ids = [1, 2]
      molecules(1)%number_of_groups = [2, 4]

      molecules(2)%groups_ids = [1, 2, 14]
      molecules(2)%number_of_groups = [1, 1, 1]

      molecules(3)%groups_ids = [1, 7, 8, 78, 79]
      molecules(3)%number_of_groups = [2, 1, 1, 3, 1]

      ! setup Dortmund model
      model = setup_dortmund(molecules)

      ! Call all Ge and derivatives
      call model%excess_gibbs(n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)

      ! Call Ge and derivatives individually
      call model%excess_gibbs(n, T, Ge_i)
      call model%excess_gibbs(n, T, GeT=GeT_i)
      call model%excess_gibbs(n, T, GeT2=GeT2_i)
      call model%excess_gibbs(n, T, Gen=Gen_i)
      call model%excess_gibbs(n, T, GeTn=GeTn_i)
      call model%excess_gibbs(n, T, Gen2=Gen2_i)

      ! Call GeModel class method
      call model%ln_activity_coefficient(n, T, ln_gammas)

      ! ========================================================================
      ! Test against Caleb Bell's implementation
      ! ------------------------------------------------------------------------
      ! Ge
      ! print *, Ge/n_t
      call assert(abs(Ge / n_t - (14.432048256607143_pr)) <= 1e-4, "Ge/n_t")

      ! ln_gammas
      ! print *, Gen/R/T, ln_gammas
      call assert(allclose(Gen / R / T, [0.727831865699655_pr, 0.49653097358793297_pr, 0.595827311587003_pr], 1e-5_pr), "Gen/R/T")
      call assert(allclose(ln_gammas, [0.727831865699655_pr, 0.49653097358793297_pr, 0.595827311587003_pr], 1e-5_pr), "ln_gammas")

      ! Gen2
      ! print *, Gen2(1,:)
      ! print *, Gen2(2,:)
      ! print *, Gen2(3,:)
      call assert(allclose(Gen2(1,:), &
         [-27.396545333109035_pr, 22.045031798841627_pr, -18.47735610933041_pr] / n_t, 1e-4_pr), "Gen2(1,:)" &
         )
      call assert(allclose(Gen2(2,:), &
         [22.04503179884173_pr, -22.7202581288098_pr, 23.170409015455117_pr] / n_t, 1e-4_pr), "Gen2(1,:)" &
         )
      call assert(allclose(Gen2(3,:), &
         [-18.47735610933015_pr, 23.17040901545502_pr, -26.299110952871235_pr] / n_t, 1e-4_pr), "Gen2(1,:)" &
         )

      ! GeT
      ! print *, GeT/n_t
      call assert(abs(GeT / n_t - 0.019210544360072263_pr) < 1e-5, "GeT/n_t")

      ! GeT2
      ! print *, GeT2 / n_t
      call assert(abs(GeT2 / n_t - (-0.000614035558615488_pr)) < 1e-5, "GeT2/n_t")

      ! GeTn
      ! print *, GeTn
      call assert(allclose(GeTn, [0.023995407535622692_pr, 0.028259971730944012_pr, 0.0009382566249192064_pr], 1e-3_pr), "GeTn")

      ! ========================================================================
      ! Test individual calls
      ! ------------------------------------------------------------------------
      call assert(abs(Ge - Ge_i) <= 1e-10, "Individual calls 1")
      call assert(abs(GeT - GeT_i) <= 1e-10, "Individual calls 2")
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "Individual calls 3")
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "Individual calls 4")
      call assert(allclose(GeTn, GeTn_i, 1e-10_pr), "Individual calls 5")

      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "Individual calls 6")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "Individual calls 7")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "Individual calls 8")

      ! ========================================================================
      ! Test pair calls
      ! ------------------------------------------------------------------------
      ! Ge
      call model%excess_gibbs(n, T, Ge=Ge_i, GeT=GeT_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "pair calls 1")
      call assert(abs(GeT - GeT_i) <= 1e-10, "pair calls 2")

      call model%excess_gibbs(n, T, Ge=Ge_i, GeT2=GeT2_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "pair calls 3")
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "pair calls 4")

      call model%excess_gibbs(n, T, Ge=Ge_i, Gen=Gen_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "pair calls 5")
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "pair calls 6")

      call model%excess_gibbs(n, T, Ge=Ge_i, GeTn=GeTn_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "pair calls 7")
      call assert(allclose(GeTn, GeTn_i, 1e-10_pr), "pair calls 8")

      call model%excess_gibbs(n, T, Ge=Ge_i, Gen2=Gen2_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "pair calls 9")
      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "pair calls 10")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "pair calls 11")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "pair calls 12")

      ! Ge_T
      call model%excess_gibbs(n, T, GeT=GeT_i, GeT2=GeT2_i)
      call assert(abs(GeT - GeT_i) <= 1e-10, "pair calls 14")
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "pair calls 15")

      call model%excess_gibbs(n, T, GeT=GeT_i, Gen=Gen_i)
      call assert(abs(GeT - GeT_i) <= 1e-10, "pair calls 16")
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "pair calls 17")

      call model%excess_gibbs(n, T, GeT=GeT_i, GeTn=GeTn_i)
      call assert(abs(GeT - GeT_i) <= 1e-10, "pair calls 18")
      call assert(allclose(GeTn, GeTn_i, 1e-10_pr), "pair calls 19")

      call model%excess_gibbs(n, T, GeT=GeT_i, Gen2=Gen2_i)
      call assert(abs(GeT - GeT_i) <= 1e-10, "pair calls 20")
      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "pair calls 21")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "pair calls 22")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "pair calls 23")

      ! Ge_T2
      call model%excess_gibbs(n, T, GeT2=GeT2_i, Gen=Gen_i)
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "pair calls 24")
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "pair calls 25")

      call model%excess_gibbs(n, T, GeT2=GeT2_i, GeTn=GeTn_i)
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "pair calls 26")
      call assert(allclose(GeTn, GeTn_i, 1e-10_pr), "pair calls 27")

      call model%excess_gibbs(n, T, GeT2=GeT2_i, Gen2=Gen2_i)
      call assert(abs(GeT2 - GeT2_i) <= 1e-10, "pair calls 28")
      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "pair calls 29")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "pair calls 30")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "pair calls 31")

      ! Gen_i
      call model%excess_gibbs(n, T, Gen=Gen_i, GeTn=GeTn_i)
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "pair calls 32")
      call assert(allclose(GeTn, GeTn_i, 1e-10_pr), "pair calls 33")

      call model%excess_gibbs(n, T, Gen=Gen_i, Gen2=Gen2_i)
      call assert(allclose(Gen, Gen_i, 1e-10_pr), "pair calls 34")
      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "pair calls 35")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "pair calls 36")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "pair calls 37")

      ! ========================================================================
      ! Just one triplet call test
      ! ------------------------------------------------------------------------
      call model%excess_gibbs(n, T, Ge=Ge_i, GeT=GeT_i, Gen2=Gen2_i)
      call assert(abs(Ge - Ge_i) <= 1e-10, "triplet calls 1")
      call assert(abs(GeT - GeT_i) <= 1e-10, "triplet calls 2")
      call assert(allclose(Gen2(1,:), Gen2_i(1,:), 1e-10_pr), "triplet calls 3")
      call assert(allclose(Gen2(2,:), Gen2_i(2,:), 1e-10_pr), "triplet calls 4")
      call assert(allclose(Gen2(3,:), Gen2_i(3,:), 1e-10_pr), "triplet calls 5")
   end subroutine test_against_caleb_thermo
end program main
