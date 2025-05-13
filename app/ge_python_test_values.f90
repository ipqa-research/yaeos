program ge_test_values
   use yaeos, only: pr, GeModel

   use yaeos, only: NRTL
   use yaeos, only: setup_unifac, UNIFAC, Groups
   use yaeos, only: setup_uniquac, UNIQUAC
   use yaeos, only: setup_dortmund
   implicit none

   integer, parameter :: nmodels = 3, nc = 3, file=10

   ! Properties
   real(pr) :: Ge, GeT, GeT2, Gen(nc), GeTn(nc), Gen2(nc, nc)
   real(pr) :: He, HeT, Hen(nc)
   real(pr) :: Se, SeT, Sen(nc)
   real(pr) :: lngamma(nc), dlngamma_dT(nc), dlngamma_dn(nc, nc)

   ! NRTL
   type(NRTL) :: nrtl_model
   real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc)

   ! UNIFAC
   type(UNIFAC) :: unifac_model
   type(Groups) :: molecules(nc)

   ! Dortmund
   type(UNIFAC) :: dortmund_model

   ! UNIQUAC
   type(UNIQUAC) :: uniquac_model
   real(pr) :: aij(nc, nc), bij(nc, nc), cij(nc, nc), dij(nc, nc), eij(nc, nc)
   real(pr) :: rs(nc), qs(nc)

   real(pr) :: n(nc), T


   ! ==========================================================================
   ! Set up the models
   ! --------------------------------------------------------------------------
   ! NRTL
   a(1, :) = [0.0_pr, -0.801_pr, -0.351_pr]
   a(2, :) = [-0.523_pr, 0.0_pr, 0.214_pr]
   a(3, :) = [0.127_pr, 0.211_pr, 0.0_pr]

   b(1, :) = [0.0_pr, -586.1_pr, 246.2_pr]
   b(2, :) = [301.2_pr, 0.0_pr, -104.2_pr]
   b(3, :) = [150.23_pr, -114.78_pr, 0.0_pr]

   c(1, :) = [0.0_pr, 0.3_pr, 0.3_pr]
   c(2, :) = [0.3_pr, 0.0_pr, 0.3_pr]
   c(3, :) = [0.3_pr, 0.3_pr, 0.0_pr]

   nrtl_model = NRTL(a, b, c)

   ! UNIFAC
   ! Hexane [CH3, CH2]
   molecules(1)%groups_ids = [1, 2]
   molecules(1)%number_of_groups = [2, 4]

   ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   ! Toluene [ACH, ACCH3]
   molecules(3)%groups_ids = [9, 11]
   molecules(3)%number_of_groups = [5, 1]

   unifac_model = setup_unifac(molecules)

   ! Dortmund
   dortmund_model = setup_dortmund(molecules)

   ! UNIQUAC
   aij(1,:) = [0.0_pr, -75.46_pr, -60.15_pr]
   aij(2,:) = [120.20_pr, 0.0_pr, 44.22_pr]
   aij(3,:) = [120.20_pr, 33.21_pr, 0.0_pr]

   bij(1,:) = [0.0_pr, -0.10062_pr, 0.2566_pr]
   bij(2,:) = [0.44835_pr, 0.0_pr, -0.01325_pr]
   bij(3,:) = [0.44835_pr, 0.124_pr, 0.0_pr]

   cij(1,:) = [0.0_pr, -0.0008052_pr, 0.00021_pr]
   cij(2,:) = [0.0004704_pr, 0.0_pr, -0.00033_pr]
   cij(3,:) = [0.0004704_pr, -0.000247_pr, 0.0_pr]

   dij(1,:) = [0.0_pr, -0.001_pr, 0.0002_pr]
   dij(2,:) = [-0.001_pr, 0.0_pr, 0.0002_pr]
   dij(3,:) = [-0.001_pr, 0.0002_pr, 0.0_pr]

   eij(1,:) = [0.0_pr, -0.00001_pr, 0.00001_pr]
   eij(2,:) = [-0.00001_pr, 0.0_pr, 0.00001_pr]
   eij(3,:) = [-0.00001_pr, 0.00001_pr, 0.0_pr]

   rs = [0.92_pr, 2.1055_pr, 1.5_pr]
   qs = [1.4_pr, 1.972_pr, 1.4_pr]

   uniquac_model = setup_uniquac(qs, rs, aij, bij, cij, dij, eij)

   ! ==========================================================================
   ! Generate test values
   ! --------------------------------------------------------------------------
   n = [15.9754_pr, 3.125_pr, 24.6721_pr]
   T = 320.0_pr

   ! Open file
   open(&
      unit=file, file="ge_test_vals.txt", status="REPLACE", action="WRITE" &
      )

   ! ==========================================================================
   ! NRTL
   ! --------------------------------------------------------------------------
   call nrtl_model%excess_gibbs(&
      n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
      )

   call nrtl_model%excess_enthalpy(&
      n, T, He=He, HeT=HeT, Hen=Hen &
      )

   call nrtl_model%excess_entropy(&
      n, T, Se=Se, SeT=SeT, Sen=Sen &
      )

   call nrtl_model%ln_activity_coefficient(&
      n, T, lngamma=lngamma, dlngammadT=dlngamma_dT, dlngammadn=dlngamma_dn&
      )

   write(file, *) "NRTL", ",", Ge, ",", GeT, ",", GeT2, ",", Gen(1), ",", &
      Gen(2), ",", Gen(3), ",", GeTn(1), ",", GeTn(2), ",", GeTn(3), ",", &
      Gen2(1, 1), ",", Gen2(1, 2), ",", Gen2(1, 3), ",", Gen2(2, 1), ",", &
      Gen2(2, 2), ",", Gen2(2, 3), ",", Gen2(3, 1), ",", Gen2(3, 2), ",", &
      Gen2(3, 3), ",", He, ",", HeT, ",", Hen(1), ",", Hen(2), ",", Hen(3), &
      ",", Se, ",", SeT, ",", Sen(1), ",", Sen(2), ",", Sen(3), ",", &
      lngamma(1), ",", lngamma(2), ",", lngamma(3), ",", &
      dlngamma_dT(1), ",", dlngamma_dT(2), ",", dlngamma_dT(3), ",", &
      dlngamma_dn(1, 1), ",", dlngamma_dn(1, 2), ",", dlngamma_dn(1, 3), ",", &
      dlngamma_dn(2, 1), ",", dlngamma_dn(2, 2), ",", dlngamma_dn(2, 3), ",", &
      dlngamma_dn(3, 1), ",", dlngamma_dn(3, 2), ",", dlngamma_dn(3, 3)

   ! ==========================================================================
   ! UNIFAC
   ! --------------------------------------------------------------------------
   call unifac_model%excess_gibbs(&
      n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
      )

   call unifac_model%excess_enthalpy(&
      n, T, He=He, HeT=HeT, Hen=Hen &
      )

   call unifac_model%excess_entropy(&
      n, T, Se=Se, SeT=SeT, Sen=Sen &
      )

   call unifac_model%ln_activity_coefficient(&
      n, T, lngamma=lngamma, dlngammadT=dlngamma_dT, dlngammadn=dlngamma_dn&
      )

   write(file, *) "UNIFAC", ",", Ge, ",", GeT, ",", GeT2, ",", Gen(1), ",", &
      Gen(2), ",", Gen(3), ",", GeTn(1), ",", GeTn(2), ",", GeTn(3), ",", &
      Gen2(1, 1), ",", Gen2(1, 2), ",", Gen2(1, 3), ",", Gen2(2, 1), ",", &
      Gen2(2, 2), ",", Gen2(2, 3), ",", Gen2(3, 1), ",", Gen2(3, 2), ",", &
      Gen2(3, 3), ",", He, ",", HeT, ",", Hen(1), ",", Hen(2), ",", Hen(3), &
      ",", Se, ",", SeT, ",", Sen(1), ",", Sen(2), ",", Sen(3), ",", &
      lngamma(1), ",", lngamma(2), ",", lngamma(3), ",", &
      dlngamma_dT(1), ",", dlngamma_dT(2), ",", dlngamma_dT(3), ",", &
      dlngamma_dn(1, 1), ",", dlngamma_dn(1, 2), ",", dlngamma_dn(1, 3), ",", &
      dlngamma_dn(2, 1), ",", dlngamma_dn(2, 2), ",", dlngamma_dn(2, 3), ",", &
      dlngamma_dn(3, 1), ",", dlngamma_dn(3, 2), ",", dlngamma_dn(3, 3)

   ! ==========================================================================
   ! Dortmund
   ! --------------------------------------------------------------------------
   call dortmund_model%excess_gibbs(&
      n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
      )

   call dortmund_model%excess_enthalpy(&
      n, T, He=He, HeT=HeT, Hen=Hen &
      )

   call dortmund_model%excess_entropy(&
      n, T, Se=Se, SeT=SeT, Sen=Sen &
      )

   call dortmund_model%ln_activity_coefficient(&
      n, T, lngamma=lngamma, dlngammadT=dlngamma_dT, dlngammadn=dlngamma_dn&
      )

   write(file, *) "Dortmund", ",", Ge, ",", GeT, ",", GeT2, ",", Gen(1), ",", &
      Gen(2), ",", Gen(3), ",", GeTn(1), ",", GeTn(2), ",", GeTn(3), ",", &
      Gen2(1, 1), ",", Gen2(1, 2), ",", Gen2(1, 3), ",", Gen2(2, 1), ",", &
      Gen2(2, 2), ",", Gen2(2, 3), ",", Gen2(3, 1), ",", Gen2(3, 2), ",", &
      Gen2(3, 3), ",", He, ",", HeT, ",", Hen(1), ",", Hen(2), ",", Hen(3), &
      ",", Se, ",", SeT, ",", Sen(1), ",", Sen(2), ",", Sen(3), ",", &
      lngamma(1), ",", lngamma(2), ",", lngamma(3), ",", &
      dlngamma_dT(1), ",", dlngamma_dT(2), ",", dlngamma_dT(3), ",", &
      dlngamma_dn(1, 1), ",", dlngamma_dn(1, 2), ",", dlngamma_dn(1, 3), ",", &
      dlngamma_dn(2, 1), ",", dlngamma_dn(2, 2), ",", dlngamma_dn(2, 3), ",", &
      dlngamma_dn(3, 1), ",", dlngamma_dn(3, 2), ",", dlngamma_dn(3, 3)

   ! ==========================================================================
   ! UNIQUAC
   ! --------------------------------------------------------------------------
   call uniquac_model%excess_gibbs(&
      n, T, Ge=Ge, GeT=GeT, GeT2=GeT2, Gen=Gen, GeTn=GeTn, Gen2=Gen2 &
      )

   call uniquac_model%excess_enthalpy(&
      n, T, He=He, HeT=HeT, Hen=Hen &
      )

   call uniquac_model%excess_entropy(&
      n, T, Se=Se, SeT=SeT, Sen=Sen &
      )

   call uniquac_model%ln_activity_coefficient(&
      n, T, lngamma=lngamma, dlngammadT=dlngamma_dT, dlngammadn=dlngamma_dn&
      )

   write(file, *) "UNIQUAC", ",", Ge, ",", GeT, ",", GeT2, ",", Gen(1), ",", &
      Gen(2), ",", Gen(3), ",", GeTn(1), ",", GeTn(2), ",", GeTn(3), ",", &
      Gen2(1, 1), ",", Gen2(1, 2), ",", Gen2(1, 3), ",", Gen2(2, 1), ",", &
      Gen2(2, 2), ",", Gen2(2, 3), ",", Gen2(3, 1), ",", Gen2(3, 2), ",", &
      Gen2(3, 3), ",", He, ",", HeT, ",", Hen(1), ",", Hen(2), ",", Hen(3), &
      ",", Se, ",", SeT, ",", Sen(1), ",", Sen(2), ",", Sen(3), ",", &
      lngamma(1), ",", lngamma(2), ",", lngamma(3), ",", &
      dlngamma_dT(1), ",", dlngamma_dT(2), ",", dlngamma_dT(3), ",", &
      dlngamma_dn(1, 1), ",", dlngamma_dn(1, 2), ",", dlngamma_dn(1, 3), ",", &
      dlngamma_dn(2, 1), ",", dlngamma_dn(2, 2), ",", dlngamma_dn(2, 3), ",", &
      dlngamma_dn(3, 1), ",", dlngamma_dn(3, 2), ",", dlngamma_dn(3, 3)

   ! Close file
   close(file)
end program ge_test_values
