program test_dortmund_parameters
   use yaeos, only: pr
   use auxiliar_functions, only: allclose, rel_error
   use testing_aux, only: assert, test_title
   implicit none

   call test_dortmund_subgroups()

contains
   subroutine test_dortmund_subgroups()
      ! OJO EL GRUPO 179 DE LA DDBST NO ESTA EN LA DORTMUND. DA IGUAL, NO SE 
      ! COMO SE USA ESE GRUPO.
      use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

      type(GeGCModelParameters) :: parameters
      integer :: sg_ids(124), mg_ids(124)
      real(pr) :: rs(124), qs(124)

      integer :: i

      parameters = DortmundParameters()

      ! ========================================================================
      ! All the UNIFAC subgroups ids randomly shuffled
      ! ------------------------------------------------------------------------
      sg_ids = [ &
         73, 9, 87, 6, 18, 197, 110, 210, &
         184, 201, 48, 70, 43, 39, 45, 124, &
         195, 209, 47, 123, 62, 50, 69, 54, &
         72, 71, 61, 153, 84, 75, 32, 52, &
         111, 29, 64, 81, 65, 220, 7, 85, &
         60, 112, 66, 38, 94, 109, 114, 15, &
         178, 36, 82, 2, 51, 40, 103, 55, &
         49, 104, 3, 28, 122, 46, 33, 68, &
         26, 17, 30, 21, 80, 16, 86, 101, &
         106, 20, 11, 12, 19, 23, 89, 42, &
         79, 76, 113, 14, 37, 8, 63, 59, &
         211, 108, 119, 31, 13, 5, 196, 74, &
         44, 35, 57, 27, 102, 93, 53, 34, &
         91, 83, 78, 22, 56, 189, 107, 24, &
         10, 92, 77, 4, 25, 105, 1, 100, &
         58, 88, 67, 41 &
         ]

      ! ========================================================================
      ! Check main group search
      ! ------------------------------------------------------------------------
      mg_ids = [ &
         39, 3, 46, 2, 9, 91, 56, 98, &
         84, 93, 22, 2, 44, 18, 21, 61, &
         89, 98, 22, 61, 31, 45, 37, 26, &
         39, 38, 30, 53, 43, 40, 15, 24, &
         56, 14, 33, 5, 34, 90, 2, 14, &
         29, 55, 34, 18, 49, 53, 55, 6, &
         84, 17, 5, 1, 23, 19, 48, 26, &
         22, 52, 1, 14, 61, 21, 15, 36, &
         13, 8, 14, 11, 42, 7, 46, 48, &
         52, 10, 4, 4, 9, 12, 46, 20, &
         42, 40, 55, 5, 18, 2, 32, 29, &
         99, 53, 53, 15, 4, 2, 90, 40, &
         21, 16, 27, 43, 48, 49, 25, 16, &
         47, 43, 42, 11, 26, 87, 53, 13, &
         3, 47, 41, 1, 13, 52, 1, 47, &
         28, 46, 35, 19 &
         ]

      do i=1,size(sg_ids)
         call assert(parameters%get_subgroup_maingroup(sg_ids(i)) == mg_ids(i), "Dortmund Maingroups")
      end do

      ! ========================================================================
      ! Test R value of each subgroup
      ! ------------------------------------------------------------------------
      rs = [ &
         2.381_pr, 0.3763_pr, 3.7543_pr, 1.2832_pr, 1.7048_pr, 3.371_pr, &
         2.687_pr, 1.5654_pr, 1.843_pr, 1.0678_pr, 1.8_pr, 1.2832_pr, &
         0.8_pr, 1.0731_pr, 0.9919_pr, 1.1589_pr, 3.9628_pr, 0.9903_pr, &
         1.8_pr, 1.3863_pr, 2.088_pr, 2.45_pr, 0.5229_pr, 2.644_pr, &
         2.0_pr, 0.8814_pr, 1.299_pr, 0.9104_pr, 1.0413_pr, 1.284_pr, &
         1.368_pr, 2.618_pr, 2.46_pr, 1.6607_pr, 1.209_pr, 1.063_pr, &
         0.9214_pr, 2.4873_pr, 1.2832_pr, 1.6607_pr, 1.535_pr, 2.42_pr, &
         1.303_pr, 1.2393_pr, 2.4617_pr, 0.9104_pr, 2.42_pr, 0.8585_pr, &
         1.3662_pr, 1.1849_pr, 0.6895_pr, 0.6325_pr, 2.65_pr, 1.5575_pr, &
         2.0767_pr, 2.5_pr, 1.8_pr, 1.7943_pr, 0.6325_pr, 1.6607_pr, &
         1.613_pr, 0.9919_pr, 1.368_pr, 1.0_pr, 1.1434_pr, 1.08_pr, &
         1.6607_pr, 1.27_pr, 0.347_pr, 1.7334_pr, 3.981_pr, 2.4748_pr, &
         1.4621_pr, 0.7173_pr, 0.91_pr, 0.91_pr, 1.7048_pr, 1.9_pr, &
         3.2994_pr, 0.8_pr, 0.3479_pr, 0.8215_pr, 2.42_pr, 1.2302_pr, &
         1.4578_pr, 1.2832_pr, 1.076_pr, 1.289_pr, 3.8183_pr, 0.683_pr, &
         1.063_pr, 1.368_pr, 0.91_pr, 1.2832_pr, 2.1094_pr, 1.284_pr, &
         0.9919_pr, 1.0746_pr, 0.4656_pr, 1.7023_pr, 2.2739_pr, 2.4617_pr, &
         0.5365_pr, 1.0746_pr, 1.4515_pr, 1.4046_pr, 0.7136_pr, 1.27_pr, &
         2.887_pr, 2.7867_pr, 1.3601_pr, 1.1434_pr, 0.3763_pr, 1.5_pr, &
         1.6_pr, 0.6325_pr, 1.1434_pr, 1.6282_pr, 0.6325_pr, 1.5_pr, &
         1.24_pr, 3.5268_pr, 3.6_pr, 1.5575_pr &
         ]

      do i=1,size(sg_ids)
         call assert(abs(parameters%get_subgroup_R(sg_ids(i)) - rs(i)) < 1e-10, "Dortmund Rs")
      end do

      ! ========================================================================
      ! Test Q value of each subgroup
      ! ------------------------------------------------------------------------
      qs = [ &
         1.522_pr, 0.4321_pr, 2.892_pr, 1.2489_pr, 1.67_pr, 2.0001_pr, &
         2.12_pr, 3.8076_pr, 1.6997_pr, 2.244_pr, 2.1473_pr, 0.4582_pr, &
         1.2742_pr, 0.353_pr, 1.0127_pr, 0.748_pr, 0.6214_pr, 3.5249_pr, &
         2.5_pr, 1.06_pr, 2.4_pr, 2.8912_pr, 0.7391_pr, 2.5_pr, &
         2.093_pr, 0.7269_pr, 1.289_pr, 0.6538_pr, 1.0116_pr, 1.098_pr, &
         1.0805_pr, 3.1836_pr, 1.808_pr, 1.3377_pr, 1.4_pr, 0.8663_pr, &
         1.3_pr, 2.4457_pr, 1.2489_pr, 0.985_pr, 1.316_pr, 2.4976_pr, &
         1.132_pr, 0.633_pr, 1.842_pr, 0.6538_pr, 2.2497_pr, 0.9938_pr, &
         0.6797_pr, 0.8067_pr, 0.8345_pr, 0.7081_pr, 2.3778_pr, 1.5193_pr, &
         1.1866_pr, 2.304_pr, 1.7946_pr, 1.34_pr, 0.3554_pr, 1.6904_pr, &
         1.368_pr, 0.66_pr, 0.7278_pr, 0.92_pr, 0.8968_pr, 0.975_pr, &
         0.985_pr, 1.6286_pr, 0.0_pr, 2.4561_pr, 3.2_pr, 1.9643_pr, &
         0.78_pr, 0.771_pr, 0.949_pr, 0.7962_pr, 1.5542_pr, 1.8_pr, &
         2.352_pr, 0.9215_pr, 0.1071_pr, 0.5135_pr, 2.0018_pr, 0.8927_pr, &
         0.9022_pr, 0.8962_pr, 0.9169_pr, 1.762_pr, 3.6018_pr, 0.3418_pr, &
         1.123_pr, 1.4332_pr, 0.3769_pr, 1.6016_pr, 2.5106_pr, 1.266_pr, &
         1.3654_pr, 0.824_pr, 0.3589_pr, 1.8784_pr, 1.5754_pr, 2.192_pr, &
         0.3177_pr, 1.176_pr, 1.248_pr, 1.4_pr, 0.8635_pr, 1.4228_pr, &
         2.241_pr, 2.7723_pr, 1.8031_pr, 1.6022_pr, 0.2113_pr, 1.08_pr, &
         0.9_pr, 0.0_pr, 1.2495_pr, 1.06_pr, 1.0608_pr, 1.08_pr, &
         1.068_pr, 2.58_pr, 2.692_pr, 1.1666_pr &
         ]

      do i=1,size(sg_ids)
         call assert(abs(parameters%get_subgroup_Q(sg_ids(i)) - qs(i)) < 1e-10, "Dortmund Qs")
      end do
   end subroutine test_dortmund_subgroups

   subroutine test_dortmund_main_groups()
      use yaeos__models_ge_group_contribution_model_parameters, only: GeGCModelParameters
      use yaeos__models_ge_group_contribution_dortmund_parameters, only: DortmundParameters

      type(GeGCModelParameters) :: parameters

      integer :: ddbst_ij(956, 2)
      real(pr) :: ddbst_aij(956, 2), ddbst_bij(956, 2), ddbst_cij(956, 2)
      real(pr) :: Aij, Aji, Bij, Bji, Cij, Cji

      integer :: i, j, k, l, isg, jsg

      parameters = DortmundParameters()

      ddbst_ij = transpose(reshape([&], [2, 956]))


      ddbst_aij = transpose(reshape([&], [2, 956]))


      ddbst_bij = transpose(reshape([&], [2, 956]))


      ddbst_cij = transpose(reshape([&], [2, 956]))


      do k=1,956
         i = ddbst_ij(k, 1)
         j = ddbst_ij(k, 2)

         Aij = parameters%get_maingroups_aij(i, j)
         Aji = parameters%get_maingroups_aij(j, i)

         Bij = parameters%get_maingroups_bij(i, j)
         Bji = parameters%get_maingroups_bij(j, i)

         Cij = parameters%get_maingroups_cij(i, j)
         Cji = parameters%get_maingroups_cij(j, i)

         call check(error, abs(Aij - ddbst_aij(k, 1)) < 1e-10_pr)
         call check(error, abs(Aji - ddbst_aij(k, 2)) < 1e-10_pr)

         call check(error, abs(Bij - ddbst_bij(k, 1)) < 1e-10_pr)
         call check(error, abs(Bji - ddbst_bij(k, 2)) < 1e-10_pr)

         call check(error, abs(Cij - ddbst_cij(k, 1)) < 1e-10_pr)
         call check(error, abs(Cji - ddbst_cij(k, 2)) < 1e-10_pr)

         ! Check the first occurrence of a subgroup with main groups i, j
         do l=1,size(parameters%subgroups_ids)
            if (parameters%get_subgroup_maingroup(parameters%subgroups_ids(l)) == i) then
               isg = parameters%subgroups_ids(l)
               cycle
            end if
         end do

         do l=1,size(parameters%subgroups_ids)
            if (parameters%get_subgroup_maingroup(parameters%subgroups_ids(l)) == j) then
               jsg = parameters%subgroups_ids(l)
               cycle
            end if
         end do

         Aij = parameters%get_subgroups_aij(isg, jsg)
         Aji = parameters%get_subgroups_aij(jsg, isg)

         Bij = parameters%get_subgroups_bij(isg, jsg)
         Bji = parameters%get_subgroups_bij(jsg, isg)

         Cij = parameters%get_subgroups_cij(isg, jsg)
         Cji = parameters%get_subgroups_cij(jsg, isg)

         call check(error, abs(Aij - ddbst_aij(k, 1)) < 1e-10_pr)
         call check(error, abs(Aji - ddbst_aij(k, 2)) < 1e-10_pr)

         call check(error, abs(Bij - ddbst_bij(k, 1)) < 1e-10_pr)
         call check(error, abs(Bji - ddbst_bij(k, 2)) < 1e-10_pr)

         call check(error, abs(Cij - ddbst_cij(k, 1)) < 1e-10_pr)
         call check(error, abs(Cji - ddbst_cij(k, 2)) < 1e-10_pr)
      end do
   end subroutine test_dortmund_main_groups
end program test_dortmund_parameters
