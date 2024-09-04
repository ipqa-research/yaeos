program cositom
   use yaeos, only: UNIQUAC, setup_uniquac, pr

   type(UNIQUAC) :: model

   real(pr) :: qs(2), rs(2), T, z(2), Ge, GeT, GeT2
   real(pr), dimension(2,2) :: A, B, C, D, E

   T = 25.0 + 273.15

   rs = [0.92, 2.1055]
   qs = [1.4, 1.972]

   A = reshape([[0.0_pr, -75.46_pr, 120.20_pr, 0.0_pr]], [2,2])
   B = reshape([[0.0_pr, -0.10062_pr, 0.44835_pr, 0.0_pr]], [2,2])
   C = reshape([[0.0_pr, -0.0008052_pr, 0.0004704_pr, 0.0_pr]], [2,2])
   D = reshape([[0.0_pr, -0.001_pr, -0.001_pr, 0.0_pr]], [2,2])
   E = reshape([[0.0_pr, -0.00001_pr, -0.00001_pr, 0.0_pr]], [2,2])

   model = setup_uniquac(qs, rs, A, B, C, D, E)

   ! Ge
   T = 298.15_pr
   z = [0.5_pr, 0.5_pr]

   call model%excess_gibbs(z, T, Ge=Ge, GeT=GeT, GeT2=GeT2)

   print *, 'Ge = ', Ge * 100
   print *, "Ge_thermo = ", -203914.9332908132_pr
   print *, "=================================================================="
   print *, "GeT = ", GeT * 100
   print *, "GeT_thermo = ", -671.8453387157984_pr
   print *, "=================================================================="
   print *, "GeT2 = ", GeT2 * 100
   print *, "GeT2_thermo = ", 0.11574726032863471_pr

end program cositom