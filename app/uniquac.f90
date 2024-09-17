program cositom
   use yaeos, only: UNIQUAC, setup_uniquac, pr, R

   type(UNIQUAC) :: model

   real(pr) :: qs(2), rs(2), T, z(2), Ge, GeT, GeT2, Gen(2), Gen2(2,2), GeTn(2)
   real(pr) :: Ge1, Ge2, delta, gammas(2)
   real(pr), dimension(2,2) :: A, B, C, D, E

   rs = [0.92_pr, 2.1055_pr]
   qs = [1.4_pr, 1.972_pr]

   A = reshape([[0.0_pr, -75.46_pr, 120.20_pr, 0.0_pr]], [2,2])
   B = reshape([[0.0_pr, -0.10062_pr, 0.44835_pr, 0.0_pr]], [2,2])
   C = reshape([[0.0_pr, -0.0008052_pr, 0.0004704_pr, 0.0_pr]], [2,2])
   D = reshape([[0.0_pr, -0.001_pr, -0.001_pr, 0.0_pr]], [2,2])
   E = reshape([[0.0_pr, -0.00001_pr, -0.00001_pr, 0.0_pr]], [2,2])

   model = setup_uniquac(qs, rs, A, B, C, D, E)

   ! Ge
   T = 298.15_pr
   z = [0.5_pr, 0.5_pr]

   delta = 1e-8_pr

   call model%excess_gibbs(z, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, Gen2=Gen2, GeTn=GeTn)
   call model%excess_gibbs(z + [delta, 0.0_pr], T, Ge=Ge1)
   call model%excess_gibbs(z + [0.0_pr, delta], T, Ge=Ge2)
   call model%ln_activity_coefficient(z, T, gammas)

   print *, 'Ge = ', Ge * 100
   print *, "Ge_thermo = ", -203914.9332908132_pr
   print *, "=================================================================="
   print *, "GeT = ", GeT * 100
   print *, "GeT_thermo = ", -671.8453387157984_pr
   print *, "=================================================================="
   print *, "GeT2 = ", GeT2 * 100
   print *, "GeT2_thermo = ", 0.11574726032863471_pr
   print *, "=================================================================="
   print *, "Gen_yaeos = ", Gen / R / T
   print *, "gammas_yaeos = ", gammas
   print *, "Gen_thermo = ", [-164.27158219423336_pr,   -0.24513258199664273_pr]
   print *, "Gen numeric = ", (Ge1 - Ge) / delta / R / T, (Ge2 - Ge) / delta / R / T
   print *, "=================================================================="
   print *, "Gen2_yaeos = "
   print *, Gen2(1,:) * 100
   print *, Gen2(2,:) * 100
   print *, "Gen2_thermo = "
   print *, [2230.622723883598_pr, -2230.6227238836045_pr]
   print *, [-2230.622723883596_pr, 2230.6227238836022_pr]

end program cositom