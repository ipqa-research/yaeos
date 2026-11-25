program test_linear_volume_translation
   use yaeos, only: pr, R, ArModel, PengRobinson76, LinearVolumeTranslation, ConstantVolumeTranslation
   use yaeos__consistency, only: numeric_ar_derivatives, ar_consistency, individual_ar_calls
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title
   implicit none

   class(ArModel), allocatable :: base
   class(ArModel), allocatable :: model
   type(LinearVolumeTranslation) :: linear
   type(ConstantVolumeTranslation) :: constant
   real(pr) :: n(4)
   real(pr) :: a(4)
   real(pr) :: b(4)
   real(pr) :: T
   real(pr) :: V
   real(pr) :: W
   real(pr) :: P
   real(pr) :: Pbase
   real(pr) :: C
   real(pr) :: CT
   real(pr) :: CT2
   real(pr) :: Cn(4)
   real(pr) :: CTn(4)
   real(pr) :: Cn2(4,4)
   real(pr) :: Ar
   real(pr) :: ArV
   real(pr) :: ArT
   real(pr) :: ArTV
   real(pr) :: ArV2
   real(pr) :: ArT2
   real(pr) :: Arn(4)
   real(pr) :: ArVn(4)
   real(pr) :: ArTn(4)
   real(pr) :: Arn2(4,4)
   real(pr) :: ArTV_combined
   real(pr) :: ArVn_combined(4)
   real(pr) :: Ar_num
   real(pr) :: ArV_num
   real(pr) :: ArT_num
   real(pr) :: ArTV_num
   real(pr) :: ArV2_num
   real(pr) :: ArT2_num
   real(pr) :: Arn_num(4)
   real(pr) :: ArVn_num(4)
   real(pr) :: ArTn_num(4)
   real(pr) :: Arn2_num(4,4)
   real(pr) :: eq31
   real(pr) :: eq33(4,4)
   real(pr) :: eq34(4)
   real(pr) :: eq36
   real(pr) :: eq37
   real(pr) :: Hr
   real(pr) :: Hrbase
   real(pr) :: Cp
   real(pr) :: Cpbase
   real(pr) :: lnphi(4)
   real(pr) :: lnphibase(4)
   logical :: passed

   write(*, *) test_title("Linear volume translation")
   base = PengRobinson76([190.564_pr, 305.32_pr, 369.83_pr, 425.12_pr], &
      [45.99_pr, 48.72_pr, 42.48_pr, 37.96_pr], [0.0115_pr, 0.0995_pr, 0.1523_pr, 0.2002_pr])
   n = [1.5_pr, 0.2_pr, 0.7_pr, 2.3_pr]
   a = [0.010_pr, 0.020_pr, -0.005_pr, 0.030_pr]
   b = [0.00001_pr, -0.00002_pr, 0.00003_pr, -0.00001_pr]
   T = 600._pr
   V = 0.8_pr
   linear = LinearVolumeTranslation(base, a, b)
   model = linear

   ! The total shift is extensive; optional outputs also work individually.
   call linear%translation(n, T, C, CT, CT2, Cn, CTn, Cn2)
   call assert(abs(C-dot_product(n,a+b*T)) < 1.e-14_pr, "C is the total shift")
   call assert(abs(CT-dot_product(n,b)) < 1.e-14_pr, "CT")
   call assert(CT2 == 0._pr .and. all(Cn2 == 0._pr), "Second derivatives vanish")
   call assert(allclose(Cn,a+b*T,1.e-12_pr), "Cn")
   call assert(allclose(CTn,b,1.e-12_pr), "CTn")
   call linear%translation(2*n, T, C)
   call assert(abs(C-2*dot_product(n,a+b*T)) < 1.e-14_pr, "Extensive C")
   call linear%translation(n, T, C, CT=CT)
   call assert(abs(CT-dot_product(n,b)) < 1.e-14_pr, "Individual CT")
   call linear%translation(n, T, C, CT2=CT2)
   call assert(CT2 == 0._pr, "Individual CT2")
   call linear%translation(n, T, C, Cn=Cn)
   call assert(allclose(Cn,a+b*T,1.e-12_pr), "Individual Cn")
   call linear%translation(n, T, C, CTn=CTn)
   call assert(allclose(CTn,b,1.e-12_pr), "Individual CTn")
   call linear%translation(n, T, C, Cn2=Cn2)
   call assert(all(Cn2 == 0._pr), "Individual Cn2")

   call individual_ar_calls(model, n, V, T, passed)
   call assert(passed, "Individual Helmholtz outputs")
   call model%residual_helmholtz(n,V,T,Ar,ArV,ArT,ArTV,ArV2,ArT2,Arn,ArVn,ArTn,Arn2)
   call model%residual_helmholtz(n,V,T,ArTV=ArTV_combined,ArVn=ArVn_combined)
   call assert(abs(ArTV-ArTV_combined)<1.e-12_pr, "Combined ArTV")
   call assert(allclose(ArVn,ArVn_combined,1.e-12_pr), "Combined ArVn")
   call numeric_ar_derivatives(model,n,V,T,d_n=2.e-4_pr,d_v=1.e-4_pr,d_t=0.05_pr, &
      Ar=Ar_num,ArV=ArV_num,ArT=ArT_num,ArTV=ArTV_num,ArV2=ArV2_num, &
      ArT2=ArT2_num,Arn=Arn_num,ArVn=ArVn_num,ArTn=ArTn_num,Arn2=Arn2_num)
   call assert(rel_error(Ar,Ar_num)<1.e-6_pr, "Ar numerical")
   call assert(rel_error(ArV,ArV_num)<1.e-5_pr, "ArV numerical")
   call assert(rel_error(ArT,ArT_num)<1.e-5_pr, "ArT numerical")
   call assert(rel_error(ArTV,ArTV_num)<1.e-4_pr, "ArTV numerical")
   call assert(rel_error(ArV2,ArV2_num)<1.e-4_pr, "ArV2 numerical")
   call assert(rel_error(ArT2,ArT2_num)<1.e-4_pr, "ArT2 numerical")
   call assert(allclose(Arn,Arn_num,1.e-5_pr), "Arn numerical")
   call assert(allclose(ArVn,ArVn_num,1.e-4_pr), "ArVn numerical")
   call assert(allclose(ArTn,ArTn_num,1.e-4_pr), "ArTn numerical")
   call assert(maxval(rel_error(Arn2,Arn2_num))<1.e-4_pr, "Arn2 numerical")
   call ar_consistency(model,n,V,T,eq31,eq33,eq34,eq36,eq37)
   call assert(abs(eq31)<1.e-12_pr .and. maxval(abs(eq33))<1.e-12_pr, "Gibbs identities")
   call assert(maxval(abs(eq34))<1.e-12_pr, "Gibbs-Duhem")
   call assert(abs(eq36)<1.e-12_pr .and. abs(eq37)<1.e-12_pr, "Pressure and temperature identities")

   ! Independent corresponding-state property relations for a linear shift.
   C = dot_product(n,a+b*T)
   CT = dot_product(n,b)
   W = V+C
   call base%pressure(n,W,T,Pbase)
   call linear%pressure(n,V,T,P)
   call assert(rel_error(P,Pbase)<1.e-12_pr, "Pressure at corresponding volumes")
   call base%lnphi_vt(n,W,T,lnPhi=lnphibase)
   call linear%lnphi_vt(n,V,T,lnPhi=lnphi)
   call assert(allclose(lnphi,lnphibase-P*(a+b*T)/(R*T),1.e-12_pr), "Fugacity correction")
   call base%enthalpy_residual_vt(n,W,T,Hr=Hrbase)
   call linear%enthalpy_residual_vt(n,V,T,Hr=Hr)
   call assert(abs(Hr-Hrbase+P*(C-T*CT))<1.e-10_pr, "Enthalpy correction")
   call base%Cp_residual_vt(n,W,T,Cpbase)
   call linear%Cp_residual_vt(n,V,T,Cp)
   call assert(abs(Cp-Cpbase)<1.e-10_pr, "Linear shift preserves Cp at fixed P")
   call test_property_derivatives(linear, n, V, T)
   call test_fugacity_derivatives(linear, n, V, T)
   call base%volume(n,1._pr,T,W,"vapor")
   call linear%volume(n,1._pr,T,V,"vapor")
   call assert(W>0._pr .and. V>0._pr, "Positive vapor volumes")
   call assert(abs(V-W+C)<1.e-12_pr, "Translated vapor volume")

   ! Zero slopes recover the constant translation, including temperature derivatives.
   linear = LinearVolumeTranslation(base,a,0*b)
   constant = ConstantVolumeTranslation(base,a)
   call linear%residual_helmholtz(n,0.8_pr,T,ArT=ArT,ArT2=ArT2,ArTn=ArTn)
   call constant%residual_helmholtz(n,0.8_pr,T,ArT=ArT_num,ArT2=ArT2_num,ArTn=ArTn_num)
   call assert(abs(ArT-ArT_num)<1.e-12_pr .and. abs(ArT2-ArT2_num)<1.e-12_pr, "Constant limit")
   call assert(maxval(abs(ArTn-ArTn_num))<1.e-12_pr, "Constant limit amount derivatives")

contains

   subroutine test_property_derivatives(eos, amounts, volume, temperature)
      ! Finite differences of residual properties at fixed V,T,n.
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: amounts(:)
      real(pr), intent(in) :: volume
      real(pr), intent(in) :: temperature
      real(pr) :: dv(4)
      real(pr) :: dt(4)
      real(pr) :: dn(size(amounts),4)
      real(pr) :: plus(4)
      real(pr) :: minus(4)
      real(pr) :: numeric(4)
      real(pr) :: n_work(size(amounts))
      real(pr) :: step
      real(pr) :: cv
      real(pr) :: cp
      real(pr) :: pressure
      real(pr) :: volume_plus
      real(pr) :: volume_minus
      real(pr) :: enthalpy_plus
      real(pr) :: enthalpy_minus
      integer :: i

      call eos%enthalpy_residual_vt(amounts,volume,temperature,HrV=dv(1),HrT=dt(1),Hrn=dn(:,1))
      call eos%gibbs_residual_vt(amounts,volume,temperature,GrV=dv(2),GrT=dt(2),Grn=dn(:,2))
      call eos%entropy_residual_vt(amounts,volume,temperature,SrV=dv(3),SrT=dt(3),Srn=dn(:,3))
      call eos%internal_energy_residual_vt(amounts,volume,temperature,UrV=dv(4),UrT=dt(4),Urn=dn(:,4))

      step = 1.e-5_pr
      call residual_values(eos,amounts,volume+step,temperature,plus)
      call residual_values(eos,amounts,volume-step,temperature,minus)
      numeric = (plus-minus)/(2*step)
      call assert(all(abs(dv-numeric)<1.e-5_pr*max(1._pr,abs(dv))), &
         "Residual H,G,S,U derivatives with respect to V")

      step = 0.01_pr
      call residual_values(eos,amounts,volume,temperature+step,plus)
      call residual_values(eos,amounts,volume,temperature-step,minus)
      numeric = (plus-minus)/(2*step)
      call assert(all(abs(dt-numeric)<1.e-6_pr*max(1._pr,abs(dt))), &
         "Residual H,G,S,U derivatives with respect to T")

      step = 1.e-5_pr
      do i = 1, size(amounts)
         n_work = amounts
         n_work(i) = amounts(i)+step
         call residual_values(eos,n_work,volume,temperature,plus)
         n_work(i) = amounts(i)-step
         call residual_values(eos,n_work,volume,temperature,minus)
         numeric = (plus-minus)/(2*step)
         call assert(all(abs(dn(i,:)-numeric)<1.e-5_pr*max(1._pr,abs(dn(i,:)))), &
            "Residual H,G,S,U amount derivatives")
      end do

      ! Cv = (dUr/dT)_V; Cp = (dHr/dT)_P. The latter needs
      ! a new volume at each perturbed temperature.
      call eos%Cv_residual_vt(amounts,volume,temperature,cv)
      call assert(abs(cv-dt(4))<1.e-12_pr, "Cv equals UrT")
      call eos%pressure(amounts,volume,temperature,pressure)
      step = 0.01_pr
      call eos%volume(amounts,pressure,temperature+step,volume_plus,"vapor")
      call eos%volume(amounts,pressure,temperature-step,volume_minus,"vapor")
      call assert(volume_plus>0._pr .and. volume_minus>0._pr, "Cp test volumes")
      call eos%enthalpy_residual_vt(amounts,volume_plus,temperature+step,Hr=enthalpy_plus)
      call eos%enthalpy_residual_vt(amounts,volume_minus,temperature-step,Hr=enthalpy_minus)
      call eos%Cp_residual_vt(amounts,volume,temperature,cp)
      call assert(abs(cp-(enthalpy_plus-enthalpy_minus)/(2*step))<1.e-6_pr, "Cp by fixed-P differences")
   end subroutine test_property_derivatives

   subroutine residual_values(eos, amounts, volume, temperature, values)
      ! Values used for the independent numerical derivatives above.
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: amounts(:)
      real(pr), intent(in) :: volume
      real(pr), intent(in) :: temperature
      real(pr), intent(out) :: values(4)

      call eos%enthalpy_residual_vt(amounts,volume,temperature,Hr=values(1))
      call eos%gibbs_residual_vt(amounts,volume,temperature,Gr=values(2))
      call eos%entropy_residual_vt(amounts,volume,temperature,Sr=values(3))
      call eos%internal_energy_residual_vt(amounts,volume,temperature,Ur=values(4))
   end subroutine residual_values

   subroutine test_fugacity_derivatives(eos, amounts, volume, temperature)
      ! Fugacity derivatives are at fixed P,T,n, not fixed V.
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: amounts(:)
      real(pr), intent(in) :: volume
      real(pr), intent(in) :: temperature
      real(pr) :: pressure
      real(pr) :: dp(size(amounts))
      real(pr) :: dt(size(amounts))
      real(pr) :: dn(size(amounts),size(amounts))
      real(pr) :: plus(size(amounts))
      real(pr) :: minus(size(amounts))
      real(pr) :: numeric(size(amounts))
      real(pr) :: work(size(amounts))
      real(pr) :: step
      integer :: i

      call eos%pressure(amounts,volume,temperature,pressure)
      call eos%lnphi_vt(amounts,volume,temperature,dlnPhidP=dp,dlnPhidT=dt,dlnPhidn=dn)
      step = 0.001_pr
      call fugacity_at_pt(eos,amounts,pressure+step,temperature,plus)
      call fugacity_at_pt(eos,amounts,pressure-step,temperature,minus)
      numeric = (plus-minus)/(2*step)
      call assert(all(abs(dp-numeric)<1.e-7_pr*max(1._pr,abs(dp))), "Fugacity pressure derivatives")

      step = 0.01_pr
      call fugacity_at_pt(eos,amounts,pressure,temperature+step,plus)
      call fugacity_at_pt(eos,amounts,pressure,temperature-step,minus)
      numeric = (plus-minus)/(2*step)
      call assert(all(abs(dt-numeric)<1.e-7_pr*max(1._pr,abs(dt))), "Fugacity temperature derivatives")

      step = 1.e-4_pr
      do i = 1, size(amounts)
         work = amounts
         work(i) = amounts(i)+step
         call fugacity_at_pt(eos,work,pressure,temperature,plus)
         work(i) = amounts(i)-step
         call fugacity_at_pt(eos,work,pressure,temperature,minus)
         numeric = (plus-minus)/(2*step)
         call assert(all(abs(dn(:,i)-numeric)<1.e-6_pr*max(1._pr,abs(dn(:,i)))), &
            "Fugacity composition derivatives")
      end do
   end subroutine test_fugacity_derivatives

   subroutine fugacity_at_pt(eos, amounts, pressure, temperature, values)
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: amounts(:)
      real(pr), intent(in) :: pressure
      real(pr), intent(in) :: temperature
      real(pr), intent(out) :: values(size(amounts))
      real(pr) :: volume

      call eos%volume(amounts,pressure,temperature,volume,"vapor")
      call assert(volume>0._pr, "Fugacity test vapor volume")
      call eos%lnphi_vt(amounts,volume,temperature,lnPhi=values)
   end subroutine fugacity_at_pt
end program test_linear_volume_translation
