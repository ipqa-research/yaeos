program main
   !! Test the calculation of critical lines
   use yaeos
   use yaeos__math, only: interpol
   use fortime, only: Timer
   use testing_aux, only: test_ok, test_title
   implicit none
   
   logical :: WRITE_FILES=.false.

   integer, parameter :: nc=12
   integer :: a_nearest
   integer :: npt=6

   type(CubicEoS) :: model
   type(EquilibriumState) :: sat, crit
   type(PTEnvel2) :: env
   type(CriticalLine) :: cl
   type(Timer) :: tim

   real(pr) ::  V, T, P, a, da

   real(pr) :: z(nc)
   real(pr) :: z0(nc)
   real(pr) :: zi(nc)

   integer :: i

   write(*, *) test_title("CRITICAL POINTS AND LINES")

   model = get_model()

   ! Calculate the composition at at specified alpha
   a = real(1, pr)/100._pr
   z = a*zi + (1-a)*z0

   ! Get the full phase envelope of the fluid
   print *, "Calculating PT envelope"
   call tim%timer_start()
   sat = saturation_temperature(model, z, P=0.00001_pr, kind="dew")
   env = pt_envelope_2ph(model, z, sat, maximum_pressure=1000._pr)
   call tim%timer_stop()

   ! Calculate the critical point
   T = sum(model%components%Tc * z)
   P = sum(model%components%Pc * z)
   call model%volume(n=z, P=P, T=T, V=V, root_type="stable")

   ! Solve a critical point
   crit = critical_point(model, z0, zi, S=a, spec=spec_CP%a, max_iters=300, a0=a)

   if (sum([crit%T, crit%P] - [env%cps(1)%T, env%cps(1)%P])**2 > 1e-2) then
      error stop "Critical point failed"
   end if
   write(*, *) test_ok("Critical point")

   ! Now test the critical lines
   print *, "CL"
   call tim%timer_start()
   cl = critical_line(model, a0=a, z0=z0, zi=zi, ns0=spec_CP%a, S0=a, dS0=0.1_pr, max_points=5000)
   call tim%timer_stop()

   ! if (WRITE_FILES) call write_cl

   da = (cl%a(size(cl%a)) - cl%a(1))/(npt+1)
   do i=1,npt
      a = cl%a(1) + da*i
      z = a*zi + (1-a)*z0
      call tim%timer_start()
      sat = saturation_temperature(model, z, P=0.0001_pr, kind="dew")
      env = pt_envelope_2ph(model, z, sat, maximum_pressure=2000._pr, points=1000, delta_0=1.5_pr)
      print *, "Running PT Envelope", i, size(env%points)
      call tim%timer_stop()
      ! if (WRITE_FILES) call write_env
      ! write(i+10, *) env

      a_nearest = minloc(abs(cl%a - a), dim=1)

      T = interpol(cl%a(a_nearest), cl%a(a_nearest+1), cl%T(a_nearest), cl%T(a_nearest+1), a)
      P = interpol(cl%a(a_nearest), cl%a(a_nearest+1), cl%P(a_nearest), cl%P(a_nearest+1), a)

      if (maxval(([T, P] - [env%cps(1)%T, env%cps(1)%P]) / [T, P]) > 1e-2) then
         write(*, *) [T, P] 
         write(*, *) [env%cps(1)%T, env%cps(1)%P]
         error stop "Critical line failed"
      end if
   end do
   write(*, *) test_ok("Critical line compared to individual PT lines")

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

   ! Write to files when internally testing
   ! subroutine write_cl
   !    do i=1,size(cl%a)
   !       write(1, *) cl%a(i), cl%T(i), cl%P(i)
   !    end do
   !    write(1, *)
   !    write(1, *)
   ! end subroutine

   ! subroutine write_env
   !    integer :: i
   !    do i=1,size(env%points)
   !       write(1, *) a, env%points(i)%T, env%points(i)%P
   !    end do
   !    write(1, *)
   !    write(1, *)
   !    call flush(1)
   ! end subroutine write_env
end program main
