module test_thermoprops
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use auxiliar_functions, only: rel_error

   implicit none

contains
   subroutine collect_suite(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("Fugacity VT", test_fugacity_VT), &
         new_unittest("Fugacity tp", test_fugacity_tp), &
         new_unittest("enthalpy_residual_vt", test_enthalpy_residual_vt), &
         new_unittest("test_gibss_residual_vt", test_gibss_residual_vt), &
         new_unittest("test_entropy_residual_vt", test_entropy_residual_vt),&
         new_unittest("test_internal_energy_vt", test_internal_energy_vt), &
         new_unittest("cp_and_cv", cp_and_cv) &
         ]
   end subroutine collect_suite

   subroutine test_fugacity_VT(error)
      use fixtures_models, only: binary_PR76
      use yaeos, only: pr, CubicEoS
      type(error_type), allocatable, intent(out) :: error
      type(CubicEoS) :: eos


      real(pr) :: lnfug(2), dlnphidp(2), dlnphidt(2), dlnphidn(2, 2)
      real(pr), allocatable :: z(:)
      real(pr) ::  v, t, p

      real(pr) :: lnfug_val(2), dlnphidp_val(2), dlnphidt_val(2)


      lnfug_val = [2.0785927055052529, -2.2825114783106386]
      dlnphidp_val = [-0.99328668293856137, -0.9965756859512391]
      dlnphidt_val = [3.0263825169536504E-002, 7.6204959316774373E-002]

      eos = binary_PR76()

      z = [0.3_pr, 0.7_pr]
      v = 8.8780451065729321E-002_pr
      t = 150

      call eos%lnphi_vt(&
         z, V, T, P, lnfug, dlnPhidP, dlnphidT, dlnPhidn &
         )
      call check( &
         error, maxval(abs(lnfug - lnfug_val)) < 1e-5 &
         )
      call check( &
         error, maxval(abs(dlnphidp - dlnphidp_val)) < 1e-5 &
         )
      call check( &
         error, maxval(abs(dlnphidt - dlnphidt_val)) < 1e-5 &
         )
   end subroutine test_fugacity_VT

   subroutine test_fugacity_TP(error)
      use fixtures_models, only: binary_PR76
      use yaeos, only: pr, CubicEoS
      type(error_type), allocatable, intent(out) :: error
      type(CubicEoS) :: eos


      real(pr) :: lnfug(2), dlnphidp(2), dlnphidt(2), dlnphidn(2, 2)

      real(pr), allocatable :: z(:)
      real(pr) ::  v, t, p

      real(pr) :: lnfug_val(2), dlnphidp_val(2), dlnphidt_val(2)

      character(len=:), allocatable :: root_type


      lnfug_val = [2.0759140949373416, -2.2851989270402058]
      dlnphidp_val = [-0.99059224575177762, -0.99388122357848807]
      dlnphidt_val = [3.0263769083149254E-002, 7.6204871541712640E-002]

      eos = binary_PR76()
      z = [0.3, 0.7]

      p = 1
      t = 150

      root_type = "liquid"
      call eos%lnphi_pt(&
         z, P, T, V, root_type, lnfug, dlnPhidP, dlnphidT, dlnPhidn&
         )

      call check(&
         error, maxval(rel_error(lnfug_val, lnfug)) < 1e-4 &
         )
      call check(&
         error, maxval(rel_error(dlnphidp_val, dlnphidp)) < 1e-4 &
         )
      call check(&
         error, maxval(rel_error(dlnphidt_val, dlnphidt)) < 1e-4 &
         )
   end subroutine test_fugacity_TP

   ! ==========================================================================
   ! Hr
   ! --------------------------------------------------------------------------
   subroutine test_enthalpy_residual_vt(error)
      use yaeos, only: pr, R, CubicEoS, PengRobinson76

      type(error_type), allocatable, intent(out) :: error

      type(CubicEoS) :: eos
      integer, parameter :: n=2

      real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)
      real(pr) :: z(n), zd1(n), zd2(n)
      real(pr) :: v, t, delta_v, delta_t, delta_n
      real(pr) :: Hr, HrT, HrV, Hrn(n) ! yaeos residuals
      real(pr) :: Hr_fromphi, HrT_num, HrV_num, Hrn_num(n) ! numeric residuals
      real(pr) :: Hr_delta_t, Hr_delta_v, Hr_delta_n1, Hr_delta_n2
      real(pr) :: lnfug(n), dlnphidt(n)

      v = 1_pr
      t = 150_pr
      delta_t = 0.0001_pr
      delta_v = 0.0001_pr
      delta_n = 0.0001_pr

      z = [0.3_pr, 0.7_pr]
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n,n])
      lij = kij / 2

      eos = PengRobinson76(tc, pc, w, kij, lij)

      ! yaeos residual enthalpies
      call eos%enthalpy_residual_vt(z, v, t, Hr, HrT=HrT, HrV=HrV, Hrn=Hrn)

      ! test against fugacity coefficient derivatives
      ! (Michelsen and Mollerup chapter 2 eq 37)
      call eos%lnphi_vt(z, v, t, lnPhi=lnfug, dlnphidt=dlnphidt)

      Hr_fromphi = -1_pr * sum(z * dlnphidT) * R * t**2 ! Hr(T,P) = Hr(T,V)

      call check(&
         error, rel_error(Hr, Hr_fromphi) < 1e-14 &
         )

      ! numeric derivatives residual enthalpies
      ! HrT_num
      call eos%enthalpy_residual_vt(z, v, t + delta_t, Hr_delta_t)

      HrT_num = (Hr_delta_t - Hr) / delta_t

      call check(&
         error, rel_error(HrT, HrT_num) < 1e-4 &
         )

      ! HrV_num
      call eos%enthalpy_residual_vt(z, v + delta_v, t, Hr_delta_v)

      HrV_num = (Hr_delta_v - Hr) / delta_v

      call check(&
         error, rel_error(HrV, HrV_num) < 1e-4 &
         )

      ! Hrn_num
      zd1 = [0.3_pr + delta_n, 0.7_pr]
      call eos%enthalpy_residual_vt(zd1, v, t, Hr_delta_n1)

      zd2 = [0.3_pr, 0.7_pr + delta_n]
      call eos%enthalpy_residual_vt(zd2, v, t, Hr_delta_n2)

      Hrn_num = [(Hr_delta_n1 - Hr) / delta_n, (Hr_delta_n2 - Hr) / delta_n]

      call check(&
         error, maxval(rel_error(Hrn, Hrn_num)) < 1e-4 &
         )

   end subroutine test_enthalpy_residual_vt

   ! ==========================================================================
   ! Gr
   ! --------------------------------------------------------------------------
   subroutine test_gibss_residual_vt(error)
      use yaeos, only: pr, R, CubicEoS, SoaveRedlichKwong

      type(error_type), allocatable, intent(out) :: error

      type(CubicEoS) :: eos
      integer, parameter :: n=2

      real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)
      real(pr) :: z(n), zd1(n), zd2(n)
      real(pr) :: Zcomp, ntot, v, t, p, delta_v, delta_t, delta_n
      real(pr) :: Gr_tp, Gr, GrT, GrV, Grn(n) ! yaeos residuals
      real(pr) :: Gr_fromphi, GrT_num, GrV_num, Grn_num(n) ! numeric residuals
      real(pr) :: Gr_delta_t, Gr_delta_v, Gr_delta_n1, Gr_delta_n2
      real(pr) :: lnfug(n), lnfugcoeffs(n)

      v = 1_pr
      t = 150_pr
      delta_t = 0.000001_pr
      delta_v = 0.000000001_pr
      delta_n = 0.000000001_pr

      z = [0.3_pr, 0.7_pr]
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n,n])
      lij = kij / 2

      eos = SoaveRedlichKwong(tc, pc, w, kij, lij)

      call eos%pressure(z, v, t, p)
      ntot = sum(z)
      Zcomp = p*v/(ntot*R*t)

      ! yaeos residual gibbs
      call eos%gibbs_residual_vt(z, V, T, Gr, GrT=GrT, GrV=GrV, Grn=Grn)

      ! test against fugacity coefficient
      ! (Michelsen and Mollerup chapter 2 eq 31)
      call eos%lnphi_vt(z, V, T, lnPhi=lnfug)

      lnfugcoeffs = lnfug

      Gr_tp = Gr - ntot*R*T*log(Zcomp) ! M and M chapter 1 Table 6

      Gr_fromphi = sum(z * lnfugcoeffs) * R * T

      call check(&
         error, rel_error(Gr_tp, Gr_fromphi) < 1e-14 &
         )

      ! numeric derivatives residual enthalpies
      ! GrT_num
      call eos%gibbs_residual_vt(z, v, t + delta_t, Gr_delta_t)

      GrT_num = (Gr_delta_t - Gr) / delta_t

      call check(&
         error, rel_error(GrT, GrT_num) < 1e-4 &
         )

      ! GrV_num
      call eos%gibbs_residual_vt(z, v + delta_v, t, Gr_delta_v)

      GrV_num = (Gr_delta_v - Gr) / delta_v

      call check(&
         error, rel_error(GrV, GrV_num) < 1e-4 &
         )

      ! Grn_num
      zd1 = [0.3_pr + delta_n, 0.7_pr]
      call eos%gibbs_residual_vt(zd1, v, t, Gr_delta_n1)

      zd2 = [0.3_pr, 0.7_pr + delta_n]
      call eos%gibbs_residual_vt(zd2, v, t, Gr_delta_n2)

      Grn_num = [(Gr_delta_n1 - Gr) / delta_n, (Gr_delta_n2 - Gr) / delta_n]

      call check(&
         error, maxval(rel_error(Grn, Grn_num)) < 1e-4 &
         )

   end subroutine test_gibss_residual_vt

   ! ==========================================================================
   ! Sr
   ! --------------------------------------------------------------------------
   subroutine test_entropy_residual_vt(error)
      use yaeos, only: pr, R, CubicEoS, SoaveRedlichKwong

      type(error_type), allocatable, intent(out) :: error

      type(CubicEoS) :: eos
      integer, parameter :: n=2

      real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)
      real(pr) :: z(n), zd1(n), zd2(n)
      real(pr) :: Zcomp, ntot, v, t, p, delta_v, delta_t, delta_n
      real(pr) :: Sr_tp, Sr_tp_hg, Sr, SrT, SrV, Srn(n) ! yaeos residuals
      real(pr) :: Hr ! yaeos residuals
      real(pr) :: Gr_tp, Gr ! yaeos residuals
      real(pr) :: SrT_num, SrV_num, Srn_num(n) ! numeric residuals
      real(pr) :: Sr_delta_t, Sr_delta_v, Sr_delta_n1, Sr_delta_n2

      v = 1_pr
      t = 150_pr
      delta_t = 0.000001_pr
      delta_v = 0.000000001_pr
      delta_n = 0.000000001_pr

      z = [0.3_pr, 0.7_pr]
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n,n])
      lij = kij / 2

      ntot = sum(z)

      eos = SoaveRedlichKwong(tc, pc, w, kij, lij)

      call eos%pressure(z, v, t, p)

      Zcomp = p*v/(ntot*R*t)

      ! yaeos residual gibbs
      call eos%entropy_residual_vt(z, v, t, Sr, SrT=SrT, SrV=SrV, Srn=Srn)
      call eos%enthalpy_residual_vt(z, v, t, Hr)
      call eos%gibbs_residual_vt(z, v, t, Gr)

      ! test against Hr and Gr
      ! (Michelsen and Mollerup chapter 2 eq 22)
      Gr_tp = Gr - ntot*R*t*log(Zcomp)
      Sr_tp_hg = (Hr - Gr_tp)/t ! Hr(T,P) = Hr(T,V)
      Sr_tp = Sr + ntot*R*log(Zcomp)

      call check(&
         error, rel_error(Sr_tp, Sr_tp_hg) < 1e-15 &
         )

      ! SrT_num
      call eos%entropy_residual_vt(z, v, t + delta_t, Sr_delta_t)

      SrT_num = (Sr_delta_t - Sr) / delta_t

      call check(&
         error, rel_error(SrT, SrT_num) < 1e-4 &
         )

      ! SrV_num
      call eos%entropy_residual_vt(z, v + delta_v, t, Sr_delta_v)

      SrV_num = (Sr_delta_v - Sr) / delta_v

      call check(&
         error, rel_error(SrV, SrV_num) < 1e-4 &
         )

      ! Srn_num
      zd1 = [0.3_pr + delta_n, 0.7_pr]
      call eos%entropy_residual_vt(zd1, v, t, Sr_delta_n1)

      zd2 = [0.3_pr, 0.7_pr + delta_n]
      call eos%entropy_residual_vt(zd2, v, t, Sr_delta_n2)

      Srn_num = [(Sr_delta_n1 - Sr) / delta_n, (Sr_delta_n2 - Sr) / delta_n]

      call check(&
         error, maxval(rel_error(Srn, Srn_num)) < 1e-4 &
         )
   end subroutine test_entropy_residual_vt

   ! ==========================================================================
   ! Ur
   ! --------------------------------------------------------------------------
   subroutine test_internal_energy_vt(error)
      use yaeos, only: pr, R, CubicEoS, SoaveRedlichKwong

      type(error_type), allocatable, intent(out) :: error

      type(CubicEoS) :: eos
      integer, parameter :: n=2

      real(pr) :: v, t, p, delta_v, delta_t, delta_n, ntot
      real(pr) :: tc(n), pc(n), w(n), kij(n, n), lij(n, n)
      real(pr) :: z(n), zd1(n), zd2(n)
      real(pr) :: Ur, UrT, UrV, Urn(n) ! yaeos residuals
      real(pr) :: Hr ! yaeos residuals
      real(pr) :: UrT_num, UrV_num, Urn_num(n) ! numeric residuals
      real(pr) :: Ur_delta_t, Ur_delta_v, Ur_delta_n1, Ur_delta_n2

      v = 1_pr
      t = 150_pr
      delta_t = 0.000001_pr
      delta_v = 0.000000001_pr
      delta_n = 0.000000001_pr

      z = [0.3_pr, 0.7_pr]
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n,n])
      lij = kij / 2

      ntot = sum(z)

      eos = SoaveRedlichKwong(tc, pc, w, kij, lij)

      ! yaeos residual U
      call eos%internal_energy_residual_vt(&
         z, v, t, Ur, UrT=UrT, UrV=UrV, Urn=Urn &
         )
      
      ! UrT_num
      call eos%internal_energy_residual_vt(z, v, t + delta_t, Ur_delta_t)

      UrT_num = (Ur_delta_t - Ur) / delta_t

      call check(error, rel_error(UrT, UrT_num) < 1e-4)

      ! UrV_num
      call eos%internal_energy_residual_vt(z, v + delta_v, t, Ur_delta_v)

      UrV_num = (Ur_delta_v - Ur) / delta_v

      call check(error, rel_error(UrV, UrV_num) < 1e-4)

      ! Urn_num
      zd1 = [0.3_pr + delta_n, 0.7_pr]
      call eos%internal_energy_residual_vt(zd1, v, t, Ur_delta_n1)

      zd2 = [0.3_pr, 0.7_pr + delta_n]
      call eos%internal_energy_residual_vt(zd2, v, t, Ur_delta_n2)

      Urn_num = [(Ur_delta_n1 - Ur) / delta_n, (Ur_delta_n2 - Ur) / delta_n]

      call check(error, maxval(rel_error(Urn, Urn_num)) < 1e-4)

      ! Test against enthalpy
      call eos%enthalpy_residual_vt(z, v, t, Hr)

      call eos%pressure(z, v, t, p)

      call check(error, rel_error(Ur, Hr - p*v + ntot*R*t) < 1e-14)
   end subroutine test_internal_energy_vt

   ! ==========================================================================
   ! Cpr and Cvr
   ! --------------------------------------------------------------------------
   subroutine cp_and_cv(error)
      ! TODO need derivatives dvdt to complete the test
      use yaeos, only: pr, R, CubicEoS, SoaveRedlichKwong

      type(error_type), allocatable, intent(out) :: error

      type(CubicEoS) :: eos
      integer, parameter :: n=2

      real(pr) :: Cpr, Cvr
      real(pr) :: z(n), tc(n), pc(n), w(n), kij(n, n), lij(n, n)
      real(pr) :: ntot, v, t, p, dpdv, dpdt

      real(pr) :: lefths, righths

      v = 1_pr
      t = 150_pr

      z = [0.3_pr, 0.7_pr]
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n,n])
      lij = kij / 2

      ntot = sum(z)

      eos = SoaveRedlichKwong(tc, pc, w, kij, lij)

      call eos%pressure(z, v, t, p, dpdv=dpdv, dpdt=dpdt)

      ! Dumb test need derivative dvdt to complete the test
      ! Michelsen and Mollerup chapter 2 eq 19
      call eos%Cp_residual_vt(z, v, t, Cpr)
      call eos%Cv_residual_vt(z, v, t, Cvr)

      lefths = (Cpr - Cvr)/R
      righths = -t/R*dpdt**2/dpdv - ntot

      call check(&
         error, rel_error(lefths, righths) < 1e-14 &
         )
   end subroutine cp_and_cv
end module test_thermoprops
