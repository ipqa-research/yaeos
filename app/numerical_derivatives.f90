program check
   use yaeos, only: pr, PengRobinson76, ArModel
   use yaeos_consistency, only: numeric_ar_derivatives
   implicit none

   class(ArModel), allocatable :: model

   real(pr) :: tc(2), pc(2), w(2)

   real(pr) :: n(2), t, v, k0(2)

   real(pr) :: Ar
   real(pr) :: ArV
   real(pr) :: ArT
   real(pr) :: Arn(size(n))
   real(pr) :: ArV2
   real(pr) :: ArT2
   real(pr) :: ArTV
   real(pr) :: ArVn(size(n))
   real(pr) :: ArTn(size(n))
   real(pr) :: Arn2(size(n), size(n))

   real(pr) :: Ar_num
   real(pr) :: ArV_num
   real(pr) :: ArT_num
   real(pr) :: Arn_num(size(n))
   real(pr) :: ArV2_num
   real(pr) :: ArT2_num
   real(pr) :: ArTV_num
   real(pr) :: ArVn_num(size(n))
   real(pr) :: ArTn_num(size(n))
   real(pr) :: Arn2_num(size(n), size(n))

   n = [0.4, 0.6]
   tc = [190.564, 425.12]
   pc = [45.99, 37.96]
   w = [0.0115478, 0.200164]
   model = PengRobinson76(tc, pc, w)

   t = 600_pr
   v = 0.5_pr

   call model%residual_helmholtz(&
      n, v, t, Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn,&
      ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2)

   call numeric_ar_derivatives(&
      model, n, v, t, d_n = 0.001_pr, d_v = 0.0001_pr, d_t = 0.01_pr, &
      Ar=Ar_num, ArV=ArV_num, ArT=ArT_num, ArTV=ArTV_num, ArV2=ArV2_num, &
      ArT2=ArT2_num, Arn=Arn_num, ArVn=ArVn_num, ArTn=ArTn_num, Arn2=Arn2_num)

    print *, "Ar: ", Ar, "Ar_num: ", Ar_num
    print *, "ArV: ",  ArV, "  ArV_num: ", ArV_num
    print *, "ArT: ",  ArT, "  ArT_num: ", ArT_num
    print *, "Arn: ",  Arn, "  Arn_num: ", Arn_num
    print *, "ArV2: ", ArV2, "  ArV2_num: ", ArV2_num
    print *, "ArT2: ", ArT2, "  ArT2_num: ", ArT2_num
    print *, "ArTV: ", ArTV, "  ArTV_num: ", ArTV_num
    print *, "ArVn: ", ArVn, "  ArVn_num: ", ArVn_num
    print *, "ArTn: ", ArTn, "  ArTn_num: ", ArTn_num
    print *, "Arn2: ", Arn2, "  Arn2_num: ", Arn2_num

end program check
