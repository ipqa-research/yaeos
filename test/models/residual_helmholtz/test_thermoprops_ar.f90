program test_thermoprops_ar
   use yaeos, only: pr, R, CubicEoS, ArModel
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title

   implicit none

   logical :: passed

   print *, test_title("Residual Thermodynamic Properties")

   ! VT residual properties test
   call test_fugacity_VT(passed)
   call assert(passed, "Fugacity VT")

   call test_enthalpy_residual_vt(passed)
   call assert(passed, "Enthalpy Residual VT")

   call test_gibss_residual_vt(passed)
   call assert(passed, "Gibbs Residual VT")

   call test_entropy_residual_vt(passed)
   call assert(passed, "Entropy Residual VT")

   call test_internal_energy_vt(passed)
   call assert(passed, "Internal Energy VT")

   call test_cp_and_cv_vt(passed)
   call assert(passed, "Cp and Cv VT")

   ! PT residual properties test
   call test_fugacity_TP(passed)
   call assert(passed, "Fugacity TP")

   call test_helmholtz_pt(passed)
   call assert(passed, "Helmholtz PT")

   call test_entropy_pt(passed)
   call assert(passed, "Entropy Residual PT")

   call test_enthalpy_pt(passed)
   call assert(passed, "Enthalpy Residual PT")

   call test_gibbs_pt(passed)
   call assert(passed, "Gibbs Residual PT")

   call test_internal_energy_pt(passed)
   call assert(passed, "Internal Energy PT")

   call test_cv_and_cp_pt(passed)
   call assert(passed, "Cv and Cp PT")

   ! Excess properties test
   call test_lngamma_pt(passed)
   call assert(passed, "lngamma PT")

   call test_excess_gibbs_pt(passed)
   call assert(passed, "Excess Gibbs PT")

   call test_excess_enthalpy_pt(passed)
   call assert(passed, "Excess Enthalpy PT")

   call test_excess_entropy_pt(passed)
   call assert(passed, "Excess Entropy PT")

   call test_excess_internal_energy_pt(passed)
   call assert(passed, "Excess Internal Energy PT")

   call test_excess_helmholtz_pt(passed)
   call assert(passed, "Excess Helmholtz PT")

   call test_excess_volume_pt(passed)
   call assert(passed, "Excess Volume PT")
contains
   ! ===========================================================================
   ! Fixtures
   ! ---------------------------------------------------------------------------
   type(CubicEoS) function binary_PR76() result(ar_model)
      use yaeos, only: PengRobinson76
      ! Auxiliar function to create a binary PR76 model for testing
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2
      ar_model = PengRobinson76(tc, pc, w, kij, lij)
   end function binary_PR76

   type(CubicEoS) function ternary_PR76() result(ar_model)
      use yaeos, only: PengRobinson76
      ! Auxiliar function to create a ternary PR76 model for testing
      integer, parameter :: nc=3

      real(pr) :: tc(nc), pc(nc), w(nc)
      real(pr) :: kij(nc, nc), lij(nc, nc)

      tc = [507.6_pr, 514.0_pr, 591.95_pr]
      pc = [30.25_pr, 61.37_pr, 57.86_pr]
      w = [0.301261_pr, 0.643558_pr, 0.466521_pr]

      kij(1, :) = [0.0_pr, 0.1_pr, 0.1_pr]
      kij(2, :) = [0.1_pr, 0.0_pr, 0.1_pr]
      kij(3, :) = [0.1_pr, 0.1_pr, 0.0_pr]

      lij(1, :) = [0.0_pr, 0.05_pr, 0.05_pr]
      lij(2, :) = [0.05_pr, 0.0_pr, 0.05_pr]
      lij(3, :) = [0.05_pr, 0.05_pr, 0.0_pr]

      ar_model = PengRobinson76(tc, pc, w, kij, lij)
   end function ternary_PR76

   ! ===========================================================================
   ! Fugacity VT
   ! ---------------------------------------------------------------------------
   subroutine test_fugacity_VT(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: n = 2

      real(pr) :: lnfug(n), dlnphidp(n), dlnphidt(n), dlnphidn(n, n), z(n)
      real(pr) :: lnfug_val(n), dlnphidp_val(n), dlnphidt_val(n)
      real(pr) ::  V, T, P

      eos = binary_PR76()

      lnfug_val = [2.0785927055052529, -2.2825114783106386]
      dlnphidp_val = [-0.99328668293856137, -0.9965756859512391]
      dlnphidt_val = [3.0263825169536504E-002, 7.6204959316774373E-002]

      z = [0.3_pr, 0.7_pr]
      V = 8.8780451065729321E-002_pr
      T = 150

      call eos%lnphi_vt(z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn)

      check = .true.

      if (maxval(abs(lnfug - lnfug_val)) > 1e-5_pr) check = .false.
      if (maxval(abs(dlnphidp - dlnphidp_val)) > 1e-5_pr) check = .false.
      if (maxval(abs(dlnphidt - dlnphidt_val)) > 1e-5_pr) check = .false.
   end subroutine test_fugacity_VT

   ! ==========================================================================
   ! Hr VT
   ! --------------------------------------------------------------------------
   subroutine test_enthalpy_residual_vt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: nc=3

      real(pr) :: n(nc), nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: v, t, delta_v, delta_t, delta_n
      real(pr) :: Hr, HrT, HrV, Hrn(nc) ! yaeos residuals
      real(pr) :: Hr_fromphi, HrT_num, HrV_num, Hrn_num(nc) ! numeric residuals

      real(pr) :: Hr_plus_delta_t, Hr_plus_delta_v
      real(pr) :: Hr_plus_delta_n1, Hr_plus_delta_n2, Hr_plus_delta_n3

      real(pr) :: Hr_minus_delta_t, Hr_minus_delta_v
      real(pr) :: Hr_minus_delta_n1, Hr_minus_delta_n2, Hr_minus_delta_n3

      real(pr) :: lnfug(nc), dlnphidt(nc)

      eos = ternary_PR76()

      v = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 7.0_pr, 4.0_pr]

      delta_t = 0.0001_pr
      delta_v = 0.00001_pr
      delta_n = 0.00001_pr

      ! yaeos residual enthalpies
      call eos%enthalpy_residual_vt(n, v, t, Hr, HrT=HrT, HrV=HrV, Hrn=Hrn)

      ! test against fugacity coefficient derivatives
      ! (Michelsen and Mollerup chapter 2 eq 37)
      call eos%lnphi_vt(n, v, t, lnPhi=lnfug, dlnphidt=dlnphidt)

      Hr_fromphi = -sum(n * dlnphidT) * R * t**2 ! Hr(T,P) = Hr(T,V)

      if (rel_error(Hr, Hr_fromphi) > 1e-14) then
         print *, "Error Hr_fromphi:", Hr_fromphi, " HrT:", HrT
         check = .false.
         return
      end if

      ! numeric derivatives residual enthalpies
      ! HrT_num
      call eos%enthalpy_residual_vt(n, v, t + delta_t, Hr_plus_delta_t)
      call eos%enthalpy_residual_vt(n, v, t - delta_t, Hr_minus_delta_t)

      HrT_num = (Hr_plus_delta_t - Hr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(HrT, HrT_num) > 1e-6) then
         print *, "Error HrT:", HrT, " HrT_num:", HrT_num
         check = .false.
         return
      end if

      ! HrV_num
      call eos%enthalpy_residual_vt(n, v + delta_v, t, Hr_plus_delta_v)
      call eos%enthalpy_residual_vt(n, v - delta_v, t, Hr_minus_delta_v)

      HrV_num = (Hr_plus_delta_v - Hr_minus_delta_v) / (2.0_pr * delta_v)

      if (rel_error(HrV, HrV_num) > 1e-6) then
         print *, "Error HrV:", HrV, " HrV_num:", HrV_num
         check = .false.
         return
      end if

      ! Hrn_num
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%enthalpy_residual_vt(nd1, v, t, Hr_plus_delta_n1)
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%enthalpy_residual_vt(nd1, v, t, Hr_minus_delta_n1)

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%enthalpy_residual_vt(nd2, v, t, Hr_plus_delta_n2)
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%enthalpy_residual_vt(nd2, v, t, Hr_minus_delta_n2)

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%enthalpy_residual_vt(nd3, v, t, Hr_plus_delta_n3)
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%enthalpy_residual_vt(nd3, v, t, Hr_minus_delta_n3)

      Hrn_num = [&
         (Hr_plus_delta_n1 - Hr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Hr_plus_delta_n2 - Hr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Hr_plus_delta_n3 - Hr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Hrn, Hrn_num, rtol=1e-6_pr)) then
         print *, "Error Hrn:", Hrn, " Hrn_num:", Hrn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_enthalpy_residual_vt

   ! ==========================================================================
   ! Gr VT
   ! --------------------------------------------------------------------------
   subroutine test_gibss_residual_vt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: nc=3

      real(pr) :: n(nc), nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: Zcomp, ntot, v, t, p, delta_v, delta_t, delta_n
      real(pr) :: Gr_tp, Gr, GrT, GrV, Grn(nc) ! yaeos residuals
      real(pr) :: Gr_fromphi, GrT_num, GrV_num, Grn_num(nc) ! numeric residuals

      real(pr) :: Gr_plus_delta_t, Gr_plus_delta_v
      real(pr) :: Gr_plus_delta_n1, Gr_plus_delta_n2, Gr_plus_delta_n3
      real(pr) :: Gr_minus_delta_t, Gr_minus_delta_v
      real(pr) :: Gr_minus_delta_n1, Gr_minus_delta_n2, Gr_minus_delta_n3

      real(pr) :: lnfug(nc)

      eos = ternary_PR76()

      n = [3.0_pr, 7.0_pr, 5.0_pr]
      p = 1.0_pr
      t = 303.15_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      delta_t = 0.000001_pr
      delta_v = 0.000001_pr
      delta_n = 0.001_pr

      ntot = sum(n)
      Zcomp = p*v/(ntot*R*t)

      ! yaeos residual gibbs
      call eos%gibbs_residual_vt(n, V, T, Gr, GrT=GrT, GrV=GrV, Grn=Grn)

      ! test against fugacity coefficient
      ! (Michelsen and Mollerup chapter 2 eq 31)
      call eos%lnphi_vt(n, V, T, lnPhi=lnfug)

      Gr_tp = Gr - ntot*R*T*log(Zcomp) ! M and M chapter 1 Table 6
      Gr_fromphi = sum(n * lnfug) * R * T

      if (rel_error(Gr_tp, Gr_fromphi) > 1e-14) then
         print *, "Error Gr_fromphi:", Gr_fromphi, " Gr_tp:", Gr_tp
         check = .false.
         return
      end if

      ! numeric derivatives residual enthalpies
      ! GrT_num
      call eos%gibbs_residual_vt(n, v, t + delta_t, Gr_plus_delta_t)
      call eos%gibbs_residual_vt(n, v, t - delta_t, Gr_minus_delta_t)

      GrT_num = (Gr_plus_delta_t - Gr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(GrT, GrT_num) > 1e-6) then
         print *, "Error GrT:", GrT, " GrT_num:", GrT_num
         check = .false.
         return
      end if

      ! GrV_num
      call eos%gibbs_residual_vt(n, v + delta_v, t, Gr_plus_delta_v)
      call eos%gibbs_residual_vt(n, v - delta_v, t, Gr_minus_delta_v)

      GrV_num = (Gr_plus_delta_v - Gr_minus_delta_v) / (2.0_pr * delta_v)

      if (rel_error(GrV, GrV_num) > 1e-6) then
         print *, "Error GrV:", GrV, " GrV_num:", GrV_num
         check = .false.
         return
      end if

      ! Grn_num
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_residual_vt(nd1, v, t, Gr_plus_delta_n1)
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_residual_vt(nd1, v, t, Gr_minus_delta_n1)

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_residual_vt(nd2, v, t, Gr_plus_delta_n2)
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_residual_vt(nd2, v, t, Gr_minus_delta_n2)

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_residual_vt(nd3, v, t, Gr_plus_delta_n3)
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_residual_vt(nd3, v, t, Gr_minus_delta_n3)

      Grn_num = [&
         (Gr_plus_delta_n1 - Gr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Gr_plus_delta_n2 - Gr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Gr_plus_delta_n3 - Gr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Grn, Grn_num, rtol=1e-6_pr)) then
         print *, "Error Grn:", Grn, " Grn_num:", Grn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_gibss_residual_vt

   ! ==========================================================================
   ! Sr VT
   ! --------------------------------------------------------------------------
   subroutine test_entropy_residual_vt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: nc=3

      real(pr) :: n(nc), nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: Zcomp, ntot, v, t, p, delta_v, delta_t, delta_n
      real(pr) :: Sr_tp, Sr_tp_hg, Sr, SrT, SrV, Srn(nc)
      real(pr) :: Hr
      real(pr) :: Gr_tp, Gr
      real(pr) :: SrT_num, SrV_num, Srn_num(nc) ! numeric residuals

      real(pr) :: Sr_plus_delta_t, Sr_plus_delta_v
      real(pr) :: Sr_plus_delta_n1, Sr_plus_delta_n2, Sr_plus_delta_n3
      real(pr) :: Sr_minus_delta_t, Sr_minus_delta_v
      real(pr) :: Sr_minus_delta_n1, Sr_minus_delta_n2, Sr_minus_delta_n3

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 7.0_pr, 2.0_pr]
      delta_t = 0.000001_pr
      delta_v = 0.000001_pr
      delta_n = 0.00001_pr

      eos = ternary_PR76()

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ntot = sum(n)

      Zcomp = p*v/(ntot*R*t)

      ! yaeos residual gibbs
      call eos%entropy_residual_vt(n, v, t, Sr, SrT=SrT, SrV=SrV, Srn=Srn)
      call eos%enthalpy_residual_vt(n, v, t, Hr)
      call eos%gibbs_residual_vt(n, v, t, Gr)

      ! test against Hr and Gr
      ! (Michelsen and Mollerup chapter 2 eq 22)
      Gr_tp = Gr - ntot*R*t*log(Zcomp)
      Sr_tp_hg = (Hr - Gr_tp)/t ! Hr(T,P) = Hr(T,V)
      Sr_tp = Sr + ntot*R*log(Zcomp)

      if (rel_error(Sr_tp, Sr_tp_hg) > 1e-14) then
         print *, "Error Sr_tp:", Sr_tp, " Sr_tp_hg:", Sr_tp_hg
         check = .false.
         return
      end if

      ! SrT_num
      call eos%entropy_residual_vt(n, v, t + delta_t, Sr_plus_delta_t)
      call eos%entropy_residual_vt(n, v, t - delta_t, Sr_minus_delta_t)

      SrT_num = (Sr_plus_delta_t - Sr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(SrT, SrT_num) > 1e-6) then
         print *, "Error SrT:", SrT, " SrT_num:", SrT_num
         check = .false.
         return
      end if

      ! SrV_num
      call eos%entropy_residual_vt(n, v + delta_v, t, Sr_plus_delta_v)
      call eos%entropy_residual_vt(n, v - delta_v, t, Sr_minus_delta_v)

      SrV_num = (Sr_plus_delta_v - Sr_minus_delta_v) / (2.0_pr * delta_v)

      if (rel_error(SrV, SrV_num) > 1e-6) then
         print *, "Error SrV:", SrV, " SrV_num:", SrV_num
         check = .false.
         return
      end if

      ! Srn_num
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%entropy_residual_vt(nd1, v, t, Sr_plus_delta_n1)
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%entropy_residual_vt(nd1, v, t, Sr_minus_delta_n1)

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%entropy_residual_vt(nd2, v, t, Sr_plus_delta_n2)
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%entropy_residual_vt(nd2, v, t, Sr_minus_delta_n2)

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%entropy_residual_vt(nd3, v, t, Sr_plus_delta_n3)
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%entropy_residual_vt(nd3, v, t, Sr_minus_delta_n3)

      Srn_num = [&
         (Sr_plus_delta_n1 - Sr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Sr_plus_delta_n2 - Sr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Sr_plus_delta_n3 - Sr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Srn, Srn_num, rtol=1e-6_pr)) then
         print *, "Error Srn:", Srn, " Srn_num:", Srn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_entropy_residual_vt

   ! ==========================================================================
   ! Ur VT
   ! --------------------------------------------------------------------------
   subroutine test_internal_energy_vt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: nc=3

      real(pr) :: n(nc), zd1(nc), zd2(nc)
      real(pr) :: Ur, UrT, UrV, Urn(nc) ! yaeos residuals
      real(pr) :: Hr ! yaeos residuals
      real(pr) :: UrT_num, UrV_num, Urn_num(nc) ! numeric residuals

      real(pr) :: ntot, p, t, v
      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: Ur_plus_delta_t, Ur_plus_delta_v
      real(pr) :: Ur_plus_delta_n1, Ur_plus_delta_n2, Ur_plus_delta_n3
      real(pr) :: Ur_minus_delta_t, Ur_minus_delta_v
      real(pr) :: Ur_minus_delta_n1, Ur_minus_delta_n2, Ur_minus_delta_n3

      real(pr) :: delta_t, delta_v, delta_n

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      delta_t = 0.00001_pr
      delta_v = 0.000001_pr
      delta_n = 0.000001_pr

      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ntot = sum(n)

      ! yaeos residual U
      call eos%internal_energy_residual_vt(&
         n, v, t, Ur, UrT=UrT, UrV=UrV, Urn=Urn &
         )

      ! UrT_num
      call eos%internal_energy_residual_vt(n, v, t + delta_t, Ur_plus_delta_t)
      call eos%internal_energy_residual_vt(n, v, t - delta_t, Ur_minus_delta_t)

      UrT_num = (Ur_plus_delta_t - Ur_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(UrT, UrT_num) > 1e-6) then
         print *, "Error UrT:", UrT, " UrT_num:", UrT_num
         check = .false.
         return
      end if

      ! UrV_num
      call eos%internal_energy_residual_vt(n, v + delta_v, t, Ur_plus_delta_v)
      call eos%internal_energy_residual_vt(n, v - delta_v, t, Ur_minus_delta_v)

      UrV_num = (Ur_plus_delta_v - Ur_minus_delta_v) / (2.0_pr * delta_v)

      if (rel_error(UrV, UrV_num) > 1e-6) then
         print *, "Error UrV:", UrV, " UrV_num:", UrV_num
         check = .false.
         return
      end if

      ! Urn_num
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%internal_energy_residual_vt(nd1, v, t, Ur_plus_delta_n1)
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%internal_energy_residual_vt(nd1, v, t, Ur_minus_delta_n1)

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%internal_energy_residual_vt(nd2, v, t, Ur_plus_delta_n2)
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%internal_energy_residual_vt(nd2, v, t, Ur_minus_delta_n2)

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%internal_energy_residual_vt(nd3, v, t, Ur_plus_delta_n3)
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%internal_energy_residual_vt(nd3, v, t, Ur_minus_delta_n3)

      Urn_num = [&
         (Ur_plus_delta_n1 - Ur_minus_delta_n1) / (2.0_pr * delta_n), &
         (Ur_plus_delta_n2 - Ur_minus_delta_n2) / (2.0_pr * delta_n), &
         (Ur_plus_delta_n3 - Ur_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Urn, Urn_num, rtol=1e-6_pr)) then
         print *, "Error Urn:", Urn, " Urn_num:", Urn_num
         check = .false.
         return
      end if

      ! Test against enthalpy
      call eos%enthalpy_residual_vt(n, v, t, Hr)

      if (rel_error(Ur, Hr - p*v + ntot*R*t) > 1e-14) then
         ! MM - Chapter 1 - Table 6
         print *, "Error Ur:", Ur, " Hr - p*v + ntot*R*t:", Hr - p*v + ntot*R*t
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_internal_energy_vt

   ! ==========================================================================
   ! Cpr and Cvr
   ! --------------------------------------------------------------------------
   subroutine test_cp_and_cv_vt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos
      integer, parameter :: nc=3

      real(pr) :: Cpr, Cvr
      real(pr) :: n(nc)
      real(pr) :: ntot, v, t, p, dpdv, dpdt, dvdt

      real(pr) :: left_hs, right_hs

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [4.0_pr, 2.0_pr, 1.5_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p, dpdv=dpdv, dpdt=dpdt)
      dvdt = -dpdt/dpdv

      ntot = sum(n)

      ! Michelsen and Mollerup chapter 2 eq 19
      call eos%Cp_residual_vt(n, v, t, Cpr)
      call eos%Cv_residual_vt(n, v, t, Cvr)

      left_hs = (Cpr - Cvr)/R
      right_hs = t / R * dpdt * dvdt - ntot

      if (rel_error(left_hs, right_hs) > 1e-14) then
         print *, "Error Cp-Cv:", (Cpr - Cvr) / R, &
            "t / R * dpdt * dvdt - ntot:", right_hs

         check = .false.
         return
      end if

      check = .true.
   end subroutine test_cp_and_cv_vt

   ! ==========================================================================
   ! Fugacity PT
   ! --------------------------------------------------------------------------
   subroutine test_fugacity_TP(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      real(pr) :: lnphip(2), lnphi(2), dlnphidp(2), dlnphidt(2), dlnphidn(2, 2)

      real(pr), allocatable :: z(:)
      real(pr) ::  v, t, p

      real(pr) :: lnphi_val(2), dlnphidp_val(2), dlnphidt_val(2)

      character(len=:), allocatable :: root_type


      lnphi_val = [2.0759140949373416, -2.2851989270402058]
      dlnphidp_val = [-0.99059224575177762, -0.99388122357848807]
      dlnphidt_val = [3.0263769083149254E-002, 7.6204871541712640E-002]

      eos = binary_PR76()
      z = [0.3, 0.7]

      p = 1
      t = 150

      root_type = "liquid"

      call eos%lnphi_pt(&
         n=z, P=P, T=T, V=V, root_type=root_type, &
         lnPhi=lnPhi,lnPhiP=lnPhiP, &
         dlnPhidP=dlnPhidP, dlnPhidT=dlnphidT, dlnPhidn=dlnPhidn&
         )

      if (maxval(rel_error(lnphi_val, lnphi)) > 1e-4) then
         print *, "Error lnphi:", lnphi, " lnphi_val:", lnphi_val
         check = .false.
         return
      end if

      if (maxval(rel_error(lnphi_val + log(P), lnPhiP)) > 1e-4) then
         print *, "Error lnPhiP:", lnPhiP, " lnphi_val + log(P):", lnphi_val + log(P)
         check = .false.
         return
      end if

      if (maxval(rel_error(dlnphidp_val, dlnphidp)) > 1e-4) then
         print *, "Error dlnphidp:", dlnphidp, " dlnphidp_val:", dlnphidp_val
         check = .false.
         return
      end if

      if (maxval(rel_error(dlnphidt_val, dlnphidt)) > 1e-4) then
         print *, "Error dlnphidt:", dlnphidt, " dlnphidt_val:", dlnphidt_val
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_fugacity_TP

   ! ==========================================================================
   ! Ar PT
   ! --------------------------------------------------------------------------
   subroutine test_helmholtz_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ar, ArT, ArP, Arn(nc)
      real(pr) :: Ar_v

      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: ntot, v, t, p, n(nc), delta_p, delta_t, delta_n, Z

      real(pr) :: ArT_num, ArP_num, Arn_num(nc) ! numeric residuals
      real(pr) :: Ar_plus_delta_t, Ar_plus_delta_p
      real(pr) :: Ar_plus_delta_n1, Ar_plus_delta_n2, Ar_plus_delta_n3
      real(pr) :: Ar_minus_delta_t, Ar_minus_delta_p
      real(pr) :: Ar_minus_delta_n1, Ar_minus_delta_n2, Ar_minus_delta_n3

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]
      delta_t = 0.00001_pr
      delta_p = 0.00001_pr
      delta_n = 0.00001_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ntot = sum(n)
      Z = p*v/(ntot*R*t)

      ! residual Ar
      call eos%helmholtz_residual_pt(&
         n, p, t, root_type="stable", Ar=Ar, ArT=ArT, ArP=ArP, Arn=Arn&
         )

      call eos%residual_helmholtz(n, v, t, Ar=Ar_v)

      if (rel_error(Ar, Ar_v - ntot * R * T * log(Z)) > 1e-14) then
         print *, "Error Ar:", Ar, &
            " Ar_v - ntot * R * T * log(Z):", Ar_v - ntot * R * T * log(Z)

         check = .false.
         return
      end if

      ! ArT_num
      call eos%helmholtz_residual_pt(&
         n, p, t + delta_t, root_type="stable", Ar=Ar_plus_delta_t &
         )
      call eos%helmholtz_residual_pt(&
         n, p, t - delta_t, root_type="stable", Ar=Ar_minus_delta_t &
         )

      ArT_num = (Ar_plus_delta_t - Ar_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(ArT, ArT_num) > 1e-6) then
         print *, "Error ArT:", ArT, " ArT_num:", ArT_num
         check = .false.
         return
      end if

      ! ArP_num
      call eos%helmholtz_residual_pt(&
         n, p + delta_p, t, root_type="stable", Ar=Ar_plus_delta_p &
         )
      call eos%helmholtz_residual_pt(&
         n, p - delta_p, t, root_type="stable", Ar=Ar_minus_delta_p &
         )

      ArP_num = (Ar_plus_delta_p - Ar_minus_delta_p) / (2.0_pr * delta_p)

      if (rel_error(ArP, ArP_num) > 1e-6) then
         print *, "Error ArP:", ArP, " ArP_num:", ArP_num
         check = .false.
         return
      end if

      ! Arn_num
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%helmholtz_residual_pt(&
         nd1, p, t, root_type="stable", Ar=Ar_plus_delta_n1 &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%helmholtz_residual_pt(&
         nd1, p, t, root_type="stable", Ar=Ar_minus_delta_n1 &
         )

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%helmholtz_residual_pt(&
         nd2, p, t, root_type="stable", Ar=Ar_plus_delta_n2 &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%helmholtz_residual_pt(&
         nd2, p, t, root_type="stable", Ar=Ar_minus_delta_n2 &
         )

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%helmholtz_residual_pt(&
         nd3, p, t, root_type="stable", Ar=Ar_plus_delta_n3 &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%helmholtz_residual_pt(&
         nd3, p, t, root_type="stable", Ar=Ar_minus_delta_n3 &
         )

      Arn_num = [&
         (Ar_plus_delta_n1 - Ar_minus_delta_n1) / (2.0_pr * delta_n), &
         (Ar_plus_delta_n2 - Ar_minus_delta_n2) / (2.0_pr * delta_n), &
         (Ar_plus_delta_n3 - Ar_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Arn, Arn_num, rtol=1e-6_pr)) then
         print *, "Error Arn:", Arn, " Arn_num:", Arn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_helmholtz_pt

   ! ==========================================================================
   ! Sr PT
   ! --------------------------------------------------------------------------
   subroutine test_entropy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Sr, SrT, SrP, Srn(nc)
      real(pr) :: Sr_v

      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: ntot, v, t, p, n(nc), delta_p, delta_t, delta_n, Z

      real(pr) :: SrT_num, SrP_num, Srn_num(nc) ! numeric residuals
      real(pr) :: Sr_plus_delta_t, Sr_plus_delta_p
      real(pr) :: Sr_plus_delta_n1, Sr_plus_delta_n2, Sr_plus_delta_n3
      real(pr) :: Sr_minus_delta_t, Sr_minus_delta_p
      real(pr) :: Sr_minus_delta_n1, Sr_minus_delta_n2, Sr_minus_delta_n3

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]
      delta_t = 0.00001_pr
      delta_p = 0.00001_pr
      delta_n = 0.00001_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ntot = sum(n)
      Z = p*v/(ntot*R*t)

      ! Sr
      call eos%entropy_residual_pt(&
         n, p, t, root_type="stable", Sr=Sr, SrP=SrP, SrT=SrT, Srn=Srn &
         )
      call eos%entropy_residual_vt(n, v, t, Sr_v)

      if (rel_error(Sr, Sr_v + ntot*R*log(Z)) > 1e-14) then
         print *, "Error Sr:", Sr, " Sr_v:", Sr_v + ntot*R*log(Z)
         check = .false.
         return
      end if

      ! SrT
      call eos%entropy_residual_pt(&
         n, p, t + delta_t, root_type="stable", Sr=Sr_plus_delta_t &
         )
      call eos%entropy_residual_pt(&
         n, p, t - delta_t, root_type="stable", Sr=Sr_minus_delta_t &
         )

      SrT_num = (Sr_plus_delta_t - Sr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(SrT, SrT_num) > 1e-6) then
         print *, "Error SrT:", SrT, " SrT_num:", SrT_num
         check = .false.
         return
      end if

      ! SrP
      call eos%entropy_residual_pt(&
         n, p + delta_p, t, root_type="stable", Sr=Sr_plus_delta_p &
         )
      call eos%entropy_residual_pt(&
         n, p - delta_p, t, root_type="stable", Sr=Sr_minus_delta_p &
         )

      SrP_num = (Sr_plus_delta_p - Sr_minus_delta_p) / (2.0_pr * delta_p)

      if (rel_error(SrP, SrP_num) > 1e-6) then
         print *, "Error SrP:", SrP, " SrP_num:", SrP_num
         check = .false.
         return
      end if

      ! Srn
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%entropy_residual_pt(&
         nd1, p, t, root_type="stable", Sr=Sr_plus_delta_n1 &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%entropy_residual_pt(&
         nd1, p, t, root_type="stable", Sr=Sr_minus_delta_n1 &
         )

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%entropy_residual_pt(&
         nd2, p, t, root_type="stable", Sr=Sr_plus_delta_n2 &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%entropy_residual_pt(&
         nd2, p, t, root_type="stable", Sr=Sr_minus_delta_n2 &
         )

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%entropy_residual_pt(&
         nd3, p, t, root_type="stable", Sr=Sr_plus_delta_n3 &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%entropy_residual_pt(&
         nd3, p, t, root_type="stable", Sr=Sr_minus_delta_n3 &
         )

      Srn_num = [&
         (Sr_plus_delta_n1 - Sr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Sr_plus_delta_n2 - Sr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Sr_plus_delta_n3 - Sr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Srn, Srn_num, rtol=1e-6_pr)) then
         print *, "Error Srn:", Srn, " Srn_num:", Srn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_entropy_pt

   ! ==========================================================================
   ! Hr PT
   ! --------------------------------------------------------------------------
   subroutine test_enthalpy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Hr, HrT, HrP, Hrn(nc)
      real(pr) :: Hr_v

      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: v, t, p, n(nc), delta_p, delta_t, delta_n, Z

      real(pr) :: HrT_num, HrP_num, Hrn_num(nc) ! numeric residuals
      real(pr) :: Hr_plus_delta_t, Hr_plus_delta_p
      real(pr) :: Hr_plus_delta_n1, Hr_plus_delta_n2, Hr_plus_delta_n3
      real(pr) :: Hr_minus_delta_t, Hr_minus_delta_p
      real(pr) :: Hr_minus_delta_n1, Hr_minus_delta_n2, Hr_minus_delta_n3

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]
      delta_t = 0.00001_pr
      delta_p = 0.00001_pr
      delta_n = 0.00001_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Hr
      call eos%enthalpy_residual_pt(&
         n, p, t, root_type="stable", Hr=Hr, &
         HrP=HrP, HrT=HrT, Hrn=Hrn &
         )

      call eos%enthalpy_residual_vt(n, v, t, Hr_v)

      if (rel_error(Hr, Hr_v) > 1e-14) then
         print *, "Error Hr:", Hr, " Hr_v:", Hr_v
         check = .false.
         return
      end if

      ! HrT
      call eos%enthalpy_residual_pt(&
         n, p, t + delta_t, root_type="stable", Hr=Hr_plus_delta_t &
         )
      call eos%enthalpy_residual_pt(&
         n, p, t - delta_t, root_type="stable", Hr=Hr_minus_delta_t &
         )

      HrT_num = (Hr_plus_delta_t - Hr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(HrT, HrT_num) > 1e-6) then
         print *, "Error HrT:", HrT, " HrT_num:", HrT_num
         check = .false.
         return
      end if

      ! HrP
      call eos%enthalpy_residual_pt(&
         n, p + delta_p, t, root_type="stable", Hr=Hr_plus_delta_p &
         )
      call eos%enthalpy_residual_pt(&
         n, p - delta_p, t, root_type="stable", Hr=Hr_minus_delta_p &
         )

      HrP_num = (Hr_plus_delta_p - Hr_minus_delta_p) / (2.0_pr * delta_p)

      if (rel_error(HrP, HrP_num) > 1e-6) then
         print *, "Error HrP:", HrP, " HrP_num:", HrP_num
         check = .false.
         return
      end if

      ! Hrn
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%enthalpy_residual_pt(&
         nd1, p, t, root_type="stable", Hr=Hr_plus_delta_n1 &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%enthalpy_residual_pt(&
         nd1, p, t, root_type="stable", Hr=Hr_minus_delta_n1 &
         )

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%enthalpy_residual_pt(&
         nd2, p, t, root_type="stable", Hr=Hr_plus_delta_n2 &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%enthalpy_residual_pt(&
         nd2, p, t, root_type="stable", Hr=Hr_minus_delta_n2 &
         )

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%enthalpy_residual_pt(&
         nd3, p, t, root_type="stable", Hr=Hr_plus_delta_n3 &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%enthalpy_residual_pt(&
         nd3, p, t, root_type="stable", Hr=Hr_minus_delta_n3 &
         )

      Hrn_num = [&
         (Hr_plus_delta_n1 - Hr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Hr_plus_delta_n2 - Hr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Hr_plus_delta_n3 - Hr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Hrn, Hrn_num, rtol=1e-6_pr)) then
         print *, "Error Hrn:", Hrn, " Hrn_num:", Hrn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_enthalpy_pt

   ! ==========================================================================
   ! Gr PT
   ! --------------------------------------------------------------------------
   subroutine test_gibbs_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Gr, GrT, GrP, Grn(nc)
      real(pr) :: Gr_v

      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: ntot, v, t, p, n(nc), delta_p, delta_t, delta_n, Z

      real(pr) :: GrT_num, GrP_num, Grn_num(nc) ! numeric residuals
      real(pr) :: Gr_plus_delta_t, Gr_plus_delta_p
      real(pr) :: Gr_plus_delta_n1, Gr_plus_delta_n2, Gr_plus_delta_n3
      real(pr) :: Gr_minus_delta_t, Gr_minus_delta_p
      real(pr) :: Gr_minus_delta_n1, Gr_minus_delta_n2, Gr_minus_delta_n3

      real(pr) :: Hr, Sr, lnphi(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]
      delta_t = 0.00001_pr
      delta_p = 0.00001_pr
      delta_n = 0.00001_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ntot = sum(n)
      Z = p * v / (ntot * R * t)

      ! Gr
      call eos%gibbs_residual_pt(&
         n, p, t, root_type="stable", Gr=Gr, GrP=GrP, GrT=GrT, Grn=Grn &
         )

      call eos%gibbs_residual_vt(n, v, t, Gr=Gr_v)

      if (rel_error(Gr, Gr_v - ntot * R * t * log(Z)) > 1e-14) then
         print *, "Error Gr:", Gr, " Gr_v:", Gr_v - ntot * R * t * log(Z)
         check = .false.
         return
      end if

      ! Gr from Sr and Hr
      call eos%entropy_residual_pt(n, p, t, root_type="stable", Sr=Sr)
      call eos%enthalpy_residual_pt(n, p, t, root_type="stable", Hr=Hr)

      if (rel_error(Gr, Hr - t * Sr) > 1e-14) then
         print *, "Error Gr:", Gr, " Hr - t * Sr:", Hr - t * Sr
         check = .false.
         return
      end if

      ! Gr from phis
      call eos%lnphi_pt(n, p, t, root_type="stable", lnPhi=lnphi)

      if (rel_error(Gr, R * t * sum(n * lnphi)) > 1e-14) then
         print *, "Error Gr:", Gr, &
            " R * t * sum(n * lnphi):", R * t * sum(n * lnphi)

         check = .false.
         return
      end if

      ! GrT
      call eos%gibbs_residual_pt(&
         n, p, t + delta_t, root_type="stable", Gr=Gr_plus_delta_t &
         )
      call eos%gibbs_residual_pt(&
         n, p, t - delta_t, root_type="stable", Gr=Gr_minus_delta_t &
         )

      GrT_num = (Gr_plus_delta_t - Gr_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(GrT, GrT_num) > 1e-6) then
         print *, "Error GrT:", GrT, " GrT_num:", GrT_num
         check = .false.
         return
      end if

      ! GrP
      call eos%gibbs_residual_pt(&
         n, p + delta_p, t, root_type="stable", Gr=Gr_plus_delta_p &
         )
      call eos%gibbs_residual_pt(&
         n, p - delta_p, t, root_type="stable", Gr=Gr_minus_delta_p &
         )

      GrP_num = (Gr_plus_delta_p - Gr_minus_delta_p) / (2.0_pr * delta_p)

      if (rel_error(GrP, GrP_num) > 1e-6) then
         print *, "Error GrP:", GrP, " GrP_num:", GrP_num
         check = .false.
         return
      end if

      ! Grn
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_residual_pt(&
         nd1, p, t, root_type="stable", Gr=Gr_plus_delta_n1 &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_residual_pt(&
         nd1, p, t, root_type="stable", Gr=Gr_minus_delta_n1 &
         )

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_residual_pt(&
         nd2, p, t, root_type="stable", Gr=Gr_plus_delta_n2 &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_residual_pt(&
         nd2, p, t, root_type="stable", Gr=Gr_minus_delta_n2 &
         )

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_residual_pt(&
         nd3, p, t, root_type="stable", Gr=Gr_plus_delta_n3 &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_residual_pt(&
         nd3, p, t, root_type="stable", Gr=Gr_minus_delta_n3 &
         )

      Grn_num = [&
         (Gr_plus_delta_n1 - Gr_minus_delta_n1) / (2.0_pr * delta_n), &
         (Gr_plus_delta_n2 - Gr_minus_delta_n2) / (2.0_pr * delta_n), &
         (Gr_plus_delta_n3 - Gr_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Grn, Grn_num, rtol=1e-6_pr)) then
         print *, "Error Grn:", Grn, " Grn_num:", Grn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_gibbs_pt

   ! ==========================================================================
   ! Ur PT
   ! --------------------------------------------------------------------------
   subroutine test_internal_energy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ur, UrT, UrP, Urn(nc)
      real(pr) :: Ur_v

      real(pr) :: nd1(nc), nd2(nc), nd3(nc)
      real(pr) :: ntot, v, t, p, n(nc), delta_p, delta_t, delta_n, Z

      real(pr) :: UrT_num, UrP_num, Urn_num(nc) ! numeric residuals
      real(pr) :: Ur_plus_delta_t, Ur_plus_delta_p
      real(pr) :: Ur_plus_delta_n1, Ur_plus_delta_n2, Ur_plus_delta_n3
      real(pr) :: Ur_minus_delta_t, Ur_minus_delta_p
      real(pr) :: Ur_minus_delta_n1, Ur_minus_delta_n2, Ur_minus_delta_n3

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]
      delta_t = 0.00001_pr
      delta_p = 0.00001_pr
      delta_n = 0.00001_pr

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Ur
      call eos%internal_energy_residual_pt(&
         n, p, t, root_type="stable", Ur=Ur, &
         UrP=UrP, UrT=UrT, Urn=Urn &
         )

      call eos%internal_energy_residual_vt(n, v, t, Ur_v)

      if (rel_error(Ur, Ur_v) > 1e-14) then
         print *, "Error Ur:", Ur, " Ur_v:", Ur_v
         check = .false.
         return
      end if

      ! UrT
      call eos%internal_energy_residual_pt(&
         n, p, t + delta_t, root_type="stable", Ur=Ur_plus_delta_t &
         )
      call eos%internal_energy_residual_pt(&
         n, p, t - delta_t, root_type="stable", Ur=Ur_minus_delta_t &
         )

      UrT_num = (Ur_plus_delta_t - Ur_minus_delta_t) / (2.0_pr * delta_t)

      if (rel_error(UrT, UrT_num) > 1e-6) then
         print *, "Error UrT:", UrT, " UrT_num:", UrT_num
         check = .false.
         return
      end if

      ! UrP
      call eos%internal_energy_residual_pt(&
         n, p + delta_p, t, root_type="stable", Ur=Ur_plus_delta_p &
         )
      call eos%internal_energy_residual_pt(&
         n, p - delta_p, t, root_type="stable", Ur=Ur_minus_delta_p &
         )

      UrP_num = (Ur_plus_delta_p - Ur_minus_delta_p) / (2.0_pr * delta_p)

      if (rel_error(UrP, UrP_num) > 1e-6) then
         print *, "Error UrP:", UrP, " UrP_num:", UrP_num
         check = .false.
         return
      end if

      ! Urn
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%internal_energy_residual_pt(&
         nd1, p, t, root_type="stable", Ur=Ur_plus_delta_n1 &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%internal_energy_residual_pt(&
         nd1, p, t, root_type="stable", Ur=Ur_minus_delta_n1 &
         )

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%internal_energy_residual_pt(&
         nd2, p, t, root_type="stable", Ur=Ur_plus_delta_n2 &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%internal_energy_residual_pt(&
         nd2, p, t, root_type="stable", Ur=Ur_minus_delta_n2 &
         )

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%internal_energy_residual_pt(&
         nd3, p, t, root_type="stable", Ur=Ur_plus_delta_n3 &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%internal_energy_residual_pt(&
         nd3, p, t, root_type="stable", Ur=Ur_minus_delta_n3 &
         )

      Urn_num = [&
         (Ur_plus_delta_n1 - Ur_minus_delta_n1) / (2.0_pr * delta_n), &
         (Ur_plus_delta_n2 - Ur_minus_delta_n2) / (2.0_pr * delta_n), &
         (Ur_plus_delta_n3 - Ur_minus_delta_n3) / (2.0_pr * delta_n) &
         ]

      if (.not. allclose(Urn, Urn_num, rtol=1e-6_pr)) then
         print *, "Error Urn:", Urn, " Urn_num:", Urn_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_internal_energy_pt

   ! ==========================================================================
   ! Cv and Cp PT
   ! --------------------------------------------------------------------------
   subroutine test_cv_and_cp_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: cv, cp
      real(pr) :: cv_v, cp_v
      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      call eos%Cv_residual_vt(n, v, t, cv_v)
      call eos%Cp_residual_vt(n, v, t, cp_v)

      call eos%Cv_residual_pt(n, p, t, root_type="stable", Cv=cv)
      call eos%Cp_residual_pt(n, p, t, root_type="stable", Cp=cp)

      if (rel_error(cv, cv_v) > 1e-14) then
         print *, "Error Cv:", cv, " Cv_v:", cv_v
         check = .false.
         return
      end if

      if (rel_error(cp, cp_v) > 1e-14) then
         print *, "Error Cp:", cp, " Cp_v:", cp_v
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_cv_and_cp_pt

   ! ==========================================================================
   ! lngamma PT
   ! --------------------------------------------------------------------------
   subroutine test_lngamma_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: lngamma(nc)
      real(pr) :: dlngammadP(nc), dlngammadT(nc), dlngammadn(nc, nc)
      real(pr) :: lngamma_ind(nc), dlngammadP_ind(nc)
      real(pr) :: dlngammadT_ind(nc), dlngammadn_ind(nc, nc)

      real(pr) :: lnphi(nc), lnphi_pure_i(nc), n_pure(nc)
      real(pr) :: lnphi_pure_temp(nc), lngamma_fromphi(nc)

      real(pr) :: p, t, v, n(nc)
      real(pr) :: delta_p, delta_t, delta_n
      real(pr) :: lng_plus_delta_p(nc), lng_minus_delta_p(nc), lngdP_num(nc)
      real(pr) :: lng_plus_delta_t(nc), lng_minus_delta_t(nc), lngdT_num(nc)
      real(pr) :: n_plus(nc), n_minus(nc), lng_plus(nc), lng_minus(nc)
      real(pr) :: lng_num(nc, nc)

      integer :: i, j

      eos = ternary_PR76()

      delta_p = 0.002_pr
      delta_t = 0.001_pr
      delta_n = 0.001_pr

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", &
         lngamma=lngamma, dlngammadP=dlngammadP, dlngammadT=dlngammadT, &
         dlngammadn=dlngammadn &
         )

      ! ln gamma
      call eos%lnphi_pt(n, p, t, root_type="stable", lnphi=lnphi)

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%lnphi_pt(&
            n_pure, p, t, root_type="stable", lnphi=lnphi_pure_temp &
            )

         lnphi_pure_i(i) = lnphi_pure_temp(i)
      end do

      lngamma_fromphi = lnphi - lnphi_pure_i

      if (.not. allclose(lngamma, lngamma_fromphi, rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      ! Individual calls
      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", lngamma=lngamma_ind &
         )
      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", dlngammadP=dlngammadP_ind &
         )
      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", dlngammadT=dlngammadT_ind &
         )
      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", dlngammadn=dlngammadn_ind &
         )

      if (.not. allclose(lngamma, lngamma_ind, rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadP, dlngammadP_ind, rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadT, dlngammadT_ind, rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadn(1,:), dlngammadn_ind(1,:), rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadn(2,:), dlngammadn_ind(2,:), rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadn(3,:), dlngammadn_ind(3,:), rtol=1e-14_pr)) then
         check = .false.
         return
      end if

      ! dP
      call eos%ln_activity_coefficient(&
         n, p + delta_p, t, root_type="stable", &
         lngamma=lng_plus_delta_p &
         )
      call eos%ln_activity_coefficient(&
         n, p - delta_p, t, root_type="stable", &
         lngamma=lng_minus_delta_p &
         )
      lngdP_num = (lng_plus_delta_p - lng_minus_delta_p) / (2.0_pr * delta_p)

      if (.not. allclose(dlngammadP, lngdP_num, rtol=1e-5_pr)) then
         check = .false.
         return
      end if


      ! dT
      call eos%ln_activity_coefficient(&
         n, p, t + delta_t, root_type="stable", &
         lngamma=lng_plus_delta_t &
         )
      call eos%ln_activity_coefficient(&
         n, p, t - delta_t, root_type="stable", &
         lngamma=lng_minus_delta_t &
         )
      lngdT_num = (lng_plus_delta_t - lng_minus_delta_t) / (2.0_pr * delta_t)

      if (.not. allclose(dlngammadT, lngdT_num, rtol=1e-6_pr)) then
         check = .false.
         return
      end if

      ! dn
      do j = 1, nc
         ! original composition
         n_plus  = n
         n_minus = n

         ! component j perturbation
         n_plus(j)  = n_plus(j)  + delta_n
         n_minus(j) = n_minus(j) - delta_n

         ! +delta evaluation
         call eos%ln_activity_coefficient(&
            n_plus, p, t, root_type="stable", &
            lngamma=lng_plus &
            )

         ! -delta evaluation
         call eos%ln_activity_coefficient(&
            n_minus, p, t, root_type="stable", &
            lngamma=lng_minus &
            )

         ! derivada numerica
         do i = 1, nc
            lng_num(i,j) = (lng_plus(i) - lng_minus(i)) / (2.0_pr * delta_n)
         end do
      end do

      if (.not. allclose(dlngammadn(1,:), lng_num(1,:), rtol=2e-7_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadn(2,:), lng_num(2,:), rtol=2e-7_pr)) then
         check = .false.
         return
      end if

      if (.not. allclose(dlngammadn(3,:), lng_num(3,:), rtol=2e-7_pr)) then
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_lngamma_pt

   ! ==========================================================================
   ! Excess Gibbs PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_gibbs_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ge, GeP, GeT, Gen(nc)
      real(pr) :: Gr_i(nc), n_pure(nc), Gr_temp, Gr

      real(pr) :: Ge_plus_delta, Ge_minus_delta
      real(pr) :: GeP_num, GeT_num, Gen_num(nc)
      real(pr) :: nd1(nc), nd2(nc), nd3(nc)

      real(pr) :: ln_gamma(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)
      real(pr) :: delta_p, delta_t, delta_n

      eos = ternary_PR76()

      delta_p = 0.001_pr
      delta_t = 0.001_pr
      delta_n = 0.001_pr

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Ge
      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", lngamma=ln_gamma &
         )

      call eos%gibbs_excess(n, p, t, root_type="stable", Ge=Ge)

      if (rel_error(Ge, R*t*sum(n*ln_gamma)) > 1e-14_pr) then
         print *, "Ge: ", Ge
         print *, "Expected: ", R*t*sum(n*ln_gamma)
         check = .false.
         return
      end if

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%gibbs_residual_pt(&
            n_pure, p, t, root_type="stable", Gr=Gr_temp &
            )

         Gr_i(i) = Gr_temp
      end do

      call eos%gibbs_residual_pt(n, p, t, root_type="stable", Gr=Gr)

      ! Test of adding pure components residual gibbs to excess gibbs
      ! and compare it to residual gibbs of the mixture
      if (rel_error(Gr, sum(n*Gr_i) + Ge) > 1e-9_pr) then
         print *, "Gr: ", Gr
         print *, "sum(n*Gr_i) + Ge: ", sum(n*Gr_i) + Ge
         check = .false.
         return
      end if

      ! Test if excess property is zero at pure conditions
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%gibbs_excess(n_pure, p, t, root_type="stable", Ge=Ge)

      if (abs(Ge) > 1e-14_pr) then
         print *, "Ge at pure conditions: ", Ge
         check = .false.
         return
      end if

      ! Derivatives
      call eos%gibbs_excess(&
         n, p, t, root_type="stable", Ge=Ge, GeP=GeP, GeT=GeT, Gen=Gen &
         )

      ! dP
      call eos%gibbs_excess(&
         n, p + delta_p, t, root_type="stable", Ge=Ge_plus_delta &
         )
      call eos%gibbs_excess(&
         n, p - delta_p, t, root_type="stable", Ge=Ge_minus_delta &
         )

      GeP_num = (Ge_plus_delta - Ge_minus_delta) / (2.0_pr * delta_p)

      if (rel_error(GeP, GeP_num) > 1e-5_pr) then
         print *, "GeP: ", GeP
         print *, "Expected: ", GeP_num
         check = .false.
         return
      end if

      ! dT
      call eos%gibbs_excess(&
         n, p, t + delta_t, root_type="stable", Ge=Ge_plus_delta &
         )
      call eos%gibbs_excess(&
         n, p, t - delta_t, root_type="stable", Ge=Ge_minus_delta &
         )

      GeT_num = (Ge_plus_delta - Ge_minus_delta) / (2.0_pr * delta_t)

      if (rel_error(GeT, GeT_num) > 1e-5_pr) then
         print *, "GeT: ", GeT
         print *, "Expected: ", GeT_num
         check = .false.
         return
      end if

      ! dn
      nd1 = n + [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_excess(&
         nd1, p, t, root_type="stable", Ge=Ge_plus_delta &
         )
      nd1 = n - [delta_n, 0.0_pr, 0.0_pr]
      call eos%gibbs_excess(&
         nd1, p, t, root_type="stable", Ge=Ge_minus_delta &
         )

      Gen_num(1) = (Ge_plus_delta - Ge_minus_delta) / (2.0_pr * delta_n)

      nd2 = n + [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_excess(&
         nd2, p, t, root_type="stable", Ge=Ge_plus_delta &
         )
      nd2 = n - [0.0_pr, delta_n, 0.0_pr]
      call eos%gibbs_excess(&
         nd2, p, t, root_type="stable", Ge=Ge_minus_delta &
         )

      Gen_num(2) = (Ge_plus_delta - Ge_minus_delta) / (2.0_pr * delta_n)

      nd3 = n + [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_excess(&
         nd3, p, t, root_type="stable", Ge=Ge_plus_delta &
         )
      nd3 = n - [0.0_pr, 0.0_pr, delta_n]
      call eos%gibbs_excess(&
         nd3, p, t, root_type="stable", Ge=Ge_minus_delta &
         )

      Gen_num(3) = (Ge_plus_delta - Ge_minus_delta) / (2.0_pr * delta_n)

      if (.not. allclose(Gen, Gen_num, rtol=1e-6_pr)) then
         print *, "Gen: ", Gen
         print *, "Expected: ", Gen_num
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_gibbs_pt

   ! ==========================================================================
   ! Excess Enthalpy PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_enthalpy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: He
      real(pr) :: Hr_i(nc), n_pure(nc), Hr_temp, Hr

      real(pr) :: lngammadT(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! He
      call eos%enthalpy_excess(n, p, t, root_type="stable", He=He)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", dlngammadT=lngammadT &
         )

      if (rel_error(He, -R*T**2*sum(n * lngammadT)) > 1e-14_pr) then
         print *, "He: ", He
         print *, "Expected: ", -R*T**2*sum(n * lngammadT)
         check = .false.
         return
      end if

      ! Check excess property as mixing
      call eos%enthalpy_residual_pt(n, p, t, root_type="stable", Hr=Hr)

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%enthalpy_residual_pt(&
            n_pure, p, t, root_type="stable", Hr=Hr_temp &
            )

         Hr_i(i) = Hr_temp
      end do

      if (rel_error(Hr, sum(n * Hr_i) + He) > 1e-10_pr) then
         print *, "Hr: ", Hr
         print *, "sum(n * Hr_i) + He: ", sum(n * Hr_i) + He
         check = .false.
         return
      end if

      ! Check excess of pure
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%enthalpy_excess(n_pure, p, t, root_type="stable", He=He)

      if (abs(He) > 1e-10_pr) then
         print *, "He: ", He
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_enthalpy_pt

   ! ==========================================================================
   ! Excess entropy PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_entropy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Se
      real(pr) :: Sr_i(nc), n_pure(nc), Sr_temp, Sr

      real(pr) :: lngamma(nc), lngammadT(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Se
      call eos%entropy_excess(n, p, t, root_type="stable", Se=Se)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", lngamma=lngamma, dlngammadT=lngammadT &
         )

      if (rel_error(Se, -R*T*sum(n * lngammadT) - R * sum(n * lngamma)) > 1e-14_pr) then
         print *, "Se: ", Se
         print *, "Expected: ", &
            -R*T*sum(n * lngammadT) - R * sum(n * lngamma)

         check = .false.
         return
      end if

      ! Check excess property as mixing
      call eos%entropy_residual_pt(n, p, t, root_type="stable", Sr=Sr)

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%entropy_residual_pt(&
            n_pure, p, t, root_type="stable", Sr=Sr_temp &
            )

         Sr_i(i) = Sr_temp
      end do

      if (rel_error(Sr, sum(n * Sr_i) + Se) > 1e-10_pr) then
         print *, "Sr: ", Sr
         print *, "sum(n * Sr_i) + Se: ", sum(n * Sr_i) + Se
         check = .false.
         return
      end if

      ! Check excess of pure
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%entropy_excess(n_pure, p, t, root_type="stable", Se=Se)

      if (abs(Se) > 1e-10_pr) then
         print *, "Se: ", Se
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_entropy_pt

   ! ==========================================================================
   ! Excess internal energy PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_internal_energy_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ue
      real(pr) :: Ur_i(nc), n_pure(nc), Ur_temp, Ur

      real(pr) :: dlngammadP(nc), dlngammadT(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Ue
      call eos%internal_energy_excess(n, p, t, root_type="stable", Ue=Ue)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", &
         dlngammadP=dlngammadP, dlngammadT=dlngammadT &
         )

      if (rel_error(Ue, -P*R*T*sum(n*dlngammadP)-R*T**2*sum(n*dlngammadT)) > 1e-14_pr) then
         print *, "Ue: ", Ue
         print *, "Expected: ", &
            -P*R*T*sum(n*dlngammadP)-R*T**2*sum(n*dlngammadT)

         check = .false.
         return
      end if

      ! Check excess property as mixing
      call eos%internal_energy_residual_pt(n, p, t, root_type="stable", Ur=Ur)

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%internal_energy_residual_pt(&
            n_pure, p, t, root_type="stable", Ur=Ur_temp &
            )

         Ur_i(i) = Ur_temp
      end do

      if (rel_error(Ur, sum(n * Ur_i) + Ue) > 1e-10_pr) then
         print *, "Ur: ", Ur
         print *, "sum(n * Ur_i) + Ue: ", sum(n * Ur_i) + Ue
         check = .false.
         return
      end if

      ! Check excess of pure
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%internal_energy_excess(n_pure, p, t, root_type="stable", Ue=Ue)

      if (abs(Ue) > 1e-10_pr) then
         print *, "Ue: ", Ue
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_internal_energy_pt

   ! ==========================================================================
   ! Excess Helholtz PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_helmholtz_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ae
      real(pr) :: Ar_i(nc), n_pure(nc), Ar_temp, Ar

      real(pr) :: dlngammadP(nc), lngamma(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Ae
      call eos%helmholtz_excess(n, p, t, root_type="stable", Ae=Ae)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", &
         lngamma=lngamma, dlngammadP=dlngammadP &
         )

      if (rel_error(Ae, R*T*sum(n*lngamma)-P*R*T*sum(n*dlngammadP)) > 1e-14_pr) then
         print *, "Ae: ", Ae
         print *, "Expected: ", &
            R*T*sum(n*lngamma)-P*R*T*sum(n*dlngammadP)

         check = .false.
         return
      end if

      ! Check excess property as mixing
      call eos%helmholtz_residual_pt(n, p, t, root_type="stable", Ar=Ar)

      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%helmholtz_residual_pt(&
            n_pure, p, t, root_type="stable", Ar=Ar_temp &
            )

         Ar_i(i) = Ar_temp
      end do

      if (rel_error(Ar, sum(n * Ar_i) + Ae) > 1e-9_pr) then
         print *, "Ar: ", Ar
         print *, "sum(n * Ar_i) + Ae: ", sum(n * Ar_i) + Ae
         check = .false.
         return
      end if

      ! Check excess of pure
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%helmholtz_excess(n_pure, p, t, root_type="stable", Ae=Ae)

      if (abs(Ae) > 1e-10_pr) then
         print *, "Ae: ", Ae
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_helmholtz_pt

   ! ==========================================================================
   ! Excess Volume PT
   ! --------------------------------------------------------------------------
   subroutine test_excess_volume_pt(check)
      logical, intent(out) :: check

      class(ArModel), allocatable :: eos

      integer, parameter :: nc=3

      real(pr) :: Ve
      real(pr) :: Vr_i(nc), n_pure(nc), Vr_temp, Vr

      real(pr) :: dlngammadP(nc)

      integer :: i

      real(pr) :: p, t, v, n(nc)

      eos = ternary_PR76()

      p = 1.0_pr
      t = 303.15_pr
      n = [3.0_pr, 4.0_pr, 2.0_pr]

      call eos%volume(n, p, t, v, root_type="stable")
      call eos%pressure(n, v, t, p)

      ! Ve
      call eos%volume_excess(n, p, t, root_type="stable", Ve=Ve)

      call eos%ln_activity_coefficient(&
         n, p, t, root_type="stable", dlngammadP=dlngammadP &
         )

      if (rel_error(Ve, R*T*sum(n*dlngammadP)) > 1e-14_pr) then
         print *, "Ve: ", Ve
         print *, "Expected: ", R*T*sum(n*dlngammadP)

         check = .false.
         return
      end if

      ! Check excess property as mixing
      do i=1,nc
         n_pure = 0.0_pr
         n_pure(i) = 1.0_pr

         call eos%volume(&
            n_pure, p, t, root_type="stable", V=Vr_temp &
            )

         Vr_i(i) = Vr_temp
      end do

      if (rel_error(v, sum(n * Vr_i) + Ve) > 1e-9_pr) then
         print *, "V: ", v
         print *, "sum(n * Vr_i) + Ve: ", sum(n * Vr_i) + Ve
         check = .false.
         return
      end if

      ! Check excess of pure
      n_pure = 0.0_pr
      n_pure(1) = 1.0_pr

      call eos%volume_excess(n_pure, p, t, root_type="stable", Ve=Ve)

      if (abs(Ve) > 1e-10_pr) then
         print *, "Ve: ", Ve
         check = .false.
         return
      end if

      check = .true.
   end subroutine test_excess_volume_pt
end program test_thermoprops_ar
