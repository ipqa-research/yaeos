program main
   !! Test the calculation of critical lines
   use yaeos
   use yaeos__math, only: solve_system
   use yaeos__equilibria_critical, only: &
      lambda1, F_critical, df_critical, CriticalLine, critical_line, critical_point
   implicit none

   integer, parameter :: nc=12

   type(CubicEoS) :: model
   type(EquilibriumState) :: sat, crit
   type(PTEnvel2) :: env
   type(CriticalLine) :: cl

   real(pr) ::  V, T, P, a

   real(pr) :: z(nc)
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)

   real(pr) :: u(nc)
   integer :: ns
   real(pr) :: S

   real(pr) :: F(3), X(3)
   integer :: i, j

   integer :: unit_cl
   integer :: unit_pt
   integer :: unit_pt_cp, unit_pt_hpl

   model = get_model()

   open(newunit=unit_cl, file="tmp/cl.dat")
   open(newunit=unit_pt, file="tmp/pt.dat")
   open(newunit=unit_pt_cp, file="tmp/pt_cp.dat")
   open(newunit=unit_pt_hpl, file="tmp/pt_hpl.dat")

   print *, "1stCL"
   a = real(1, pr)/100._pr
   cl = critical_line(model, a0=a, z0=z0, zi=zi, ns=spec_CP%a, S=a, dS0=0.1_pr)
   print *, size(cl%a)
   do i=1, size(cl%a)
      write(unit_cl, *) cl%a(i), cl%V(i), cl%T(i), cl%P(i)
   end do
   write (unit_cl, *)
   write (unit_cl, *)

   print *, "2ndCL"
   a = 1-epsilon(1._pr)
   cl = critical_line(model, a0=a, z0=z0, zi=zi, ns=spec_CP%a, S=a, dS0=-0.01_pr)
   print *, size(cl%a)
   do i=1, size(cl%a)
      write(unit_cl, *) cl%a(i), cl%V(i), cl%T(i), cl%P(i)
   end do

   z = a*zi + (1-a)*z0
   T = sum(model%components%Tc * z)
   P = sum(model%components%Pc * z)
   call model%volume(n=z, P=P, T=T, V=V, root_type="stable")
   
   X = [a, log(V), log(T)]
   ns = 1
   S = X(ns)

   a = 0.9
   z = a*zi + (1-a)*z0

   sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
   env = pt_envelope_2ph(model, z, sat)
   write(20, *) env
   write(21, *) env

   !$OMP PARALLEL DO PRIVATE(j, a, z, sat, env, i) shared(model, z0, zi)
   do j=1, 99, 3
      a = real(j, pr)/100
      z = a*zi + (1-a)*z0
      sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
      env = pt_envelope_2ph(model, z, sat)

      do i=1,size(env%points)
         write(unit_pt, *) a, env%points(i)%T, env%points(i)%P
      end do
      write(unit_pt, *)
      write(unit_pt, *)
      
      write(unit_pt_cp, *) a, env%cps

      env = find_hpl(model, z, t0=500._pr, P0=1000._pr)
      do i=1,size(env%points)
         write(unit_pt_hpl, *) a, env%points(i)%T, env%points(i)%P
      end do
      write(unit_pt_hpl, *)
      write(unit_pt_hpl, *)
   end do
   !$OMP END PARALLEL DO

   close(unit_cl)
   close(unit_pt)
   close(unit_pt_cp)
   close(unit_pt_hpl)
contains

   type(CubicEoS) function get_model()
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)
      z0=[0.0656,0.3711,0.0538,0.0373,0.0261,0.0187,0.0218,0.1791,0.091,0.0605,0.0447,0.0302]
      ! zi=[0.1775,0.3878,0.188,0.2196,0.0271,0.,0.,0.,0.,0.,0.,0.]
      ! zi=[0.,-0.1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
      zi=[1.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]

      tc=[304.088888888889,190.6,305.4,369.8,425.2,469.6,507.4,616.2,698.9,770.4,853.1,1001.2]
      pc=[73.7343491450634,45.9196083838941,48.7516547159404,42.3795504688362, &
          37.9291919470491,33.6811224489796,29.6353419746277,28.8261858797573,&
          19.3186017650303,16.5876999448428,15.2728212906784,14.6659542195256]
      w= [0.228,0.008,0.098,0.152,0.193,0.251,0.296,0.454,0.787,1.048,1.276,1.299]
      kij = 0
      kij(1, 2) = 0.12
      kij(1, 3:) = 0.15
      kij(:, 1) = kij(1, :)

      get_model = PengRobinson78(tc, pc, w, kij=kij)
   end function get_model
   ! type(CubicEoS) function get_model()
   !    use yaeos__models, only: SoaveRedlichKwong
   !    real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)
   !    ! Tc=  [190.564, 304.21, 617.7]
   !    ! Pc=  [45.99, 73.83000000000001, 21.1]
   !    ! w=  [0.0115478, 0.223621, 0.492328]
   ! z0 = [0.0, 0.4, 0.3, 0.2, 0.1]
   ! zi = [0.3, 0.7, 0.0, 0.0, 0.0]
   !
   !    Tc=  [304.21, 190.564, 425.12, 617.7, 874.0]
   !    Pc=  [73.83000000000001, 45.99, 37.96, 21.1, 6.800000000000001]
   !    w=  [0.223621, 0.0115478, 0.200164, 0.492328, 1.52596]
   !
   !    kij = 0
   !    kij(1, :) = 0.12
   !    kij(:, 1) = 0.12
   !    get_model = SoaveRedlichKwong(tc, pc, w, kij=kij)
   ! end function get_model
end program main
