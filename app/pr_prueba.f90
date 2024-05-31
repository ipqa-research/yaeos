program pruebita
   use yaeos, only : ArModel, EquilibriaState, saturation_temperature
   use yaeos, only : pr, PengRobinson76

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: sat_t

   real(pr) :: z(2)
   real(pr) :: x_bubble_nokij(18), t_bubble_nokij(18)
   real(pr) :: y_dew_nokij(18), t_dew_nokij(18)

   integer :: i

   model = PengRobinson76(&
      [512.5_pr, 562.05_pr], &
      [80.84_pr, 48.95_pr], &
      [0.565831_pr, 0.2103_pr] &
      ! kij=reshape([0.0_pr, 0.084_pr, 0.084_pr, 0.0_pr], shape=[2, 2]) &
      )

   x_bubble_nokij = [&
      0.010, 0.036, 0.071, 0.118, 0.163, 0.226, &
      0.282, 0.341, 0.410, 0.481, 0.546, 0.626, &
      0.720, 0.796, 0.864, 0.903, 0.948, 0.984  &
   ]

   t_bubble_nokij = [&
      352.584, 351.494, 350.013, 348.455, 347.052, 345.26, &
      343.857, 342.766, 341.675, 340.584, 339.727, 338.87, &
      338.247, 337.857, 337.701, 337.701, 337.701, 337.779 &
   ]

   y_dew_nokij = [&
      0.010, 0.032, 0.058, 0.097, 0.137, 0.182, &
      0.231, 0.270, 0.370, 0.442, 0.492, 0.533, &
      0.583, 0.649, 0.713, 0.809, 0.864, 0.967  &
   ]

   t_dew_nokij = [&
      352.896, 352.429, 351.883, 351.182, 350.325, 349.390, &
      348.377, 347.597, 345.571, 344.013, 343.000, 342.221, &
      341.286, 340.117, 339.104, 338.013, 337.701, 337.623  &
   ]

   ! x_exp = [1e-6, 0.026, 0.050, 0.088, 0.164, 0.333, 0.549, 0.699, 0.782, 0.898, 0.973, 0.9999999]
   ! y_exp = [1e-6, 0.267, 0.371, 0.457, 0.526, 0.559, 0.595, 0.633, 0.665, 0.760, 0.907, 0.9999999]
   ! t_exp = [353.25, 343.82, 339.59, 336.02, 333.35, 331.79, 331.17, 331.25, 331.62, 333.05, 335.85, 337.85]

   print *, "================================================================="
   print *, "Bubble booble"
   print *, "================================================================="
   do i = 1, 18
      z = [x_bubble_nokij(i), 1.0_pr - x_bubble_nokij(i)]
      sat_t = saturation_temperature(model, z, p=1.01325_pr, kind="bubble", t0=t_bubble_nokij(i))

      print *, "T experimental: ", t_bubble_nokij(i), "T yaeos: ", sat_t%T
   end do

   print *, "================================================================="
   print *, "Dew"
   print *, "================================================================="
   do i = 1, 18
      z = [y_dew_nokij(i), 1.0_pr - y_dew_nokij(i)]
      sat_t = saturation_temperature(model, z, p=1.01325_pr, kind="dew", t0=t_dew_nokij(i))

      print *, "T experimental: ", t_dew_nokij(i), "T yaeos: ", sat_t%T
   end do
   
end program pruebita
