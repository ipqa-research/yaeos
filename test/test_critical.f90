program main
   !! Test the calculation of critical lines
   use yaeos
   use yaeos__equilibria_critical, only:  CriticalLine, critical_line, solve_cp
   implicit none

   integer, parameter :: nc=12

   type(CubicEoS) :: model
   type(EquilibriumState) :: sat
   type(PTEnvel2) :: env
   type(CriticalLine) :: cl

   real(pr) ::  V, T, P, a

   real(pr) :: z(nc)
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)

   real(pr) :: u(nc)
   integer :: ns
   real(pr) :: S

   real(pr) :: X(3)
   integer :: i, j

   model = get_model()

   ! Calculate the composition at at specified alpha
   a = real(1, pr)/100._pr
   z = a*zi + (1-a)*z0

   ! Get the full phase envelope of the fluid
   sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
   env = pt_envelope_2ph(model, z, sat)

   ! Calculate the critical point
   T = sum(model%components%Tc * z)
   P = sum(model%components%Pc * z)
   call model%volume(n=z, P=P, T=T, V=V, root_type="stable")

   ! Solve a critical point
   X = [a, log(V), log(T)]
   u = [(1, i=1, size(z0))]
   u = u/sum(u)
   call solve_cp(model, X, 1, X(1), z0, zi, u)
   V = exp(X(2)); T = exp(X(3))
   call model%pressure(z, V=V, T=T, P=P)

   if (sum([T, P] - [env%cps(1)%T, env%cps(1)%P])**2 > 1e-2) then
      print *, "Critical point failed"
      stop 1
   end if

   ! Now test the critical lines
   cl = critical_line(model, a, z0, zi, 0.01_pr)
   do i=1,5
      a = cl%a(i)
      z = a*zi + (1-a)*z0
      sat = saturation_temperature(model, z, P=0.01_pr, kind="dew")
      env = pt_envelope_2ph(model, z, sat)

      if (sum(([cl%T(i), cl%P(i)] - [env%cps(1)%T, env%cps(1)%P]))**2 > 1e-2) then
         print *, "Critical line failed"
         stop 1
      end if
   end do

contains
   type(CubicEoS) function get_model()
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)
      z0=[0.0656,0.3711,0.0538,0.0373,0.0261,0.0187,&
         0.0218,0.1791,0.091,0.0605,0.0447,0.0302]
      zi=[1.0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]

      tc=[304.088888888889,190.6,305.4,369.8,425.2,469.6,507.4,616.2,&
          698.9,770.4,853.1,1001.2]
      pc=[73.7343491450634,45.9196083838941,48.7516547159404,42.3795504688362, &
         37.9291919470491,33.6811224489796,29.6353419746277,28.8261858797573,&
         19.3186017650303,16.5876999448428,15.2728212906784,14.6659542195256]
      w= [0.228,0.008,0.098,0.152,0.193,0.251,0.296,&
          0.454,0.787,1.048,1.276,1.299]
      kij = 0
      kij(1, 2) = 0.12
      kij(1, 3:) = 0.15
      kij(:, 1) = kij(1, :)

      get_model = PengRobinson78(tc, pc, w, kij=kij)
   end function get_model
end program main