program main
   use yaeos, only: pr, CubicEoS, SoaveRedlichKwong, CriticalLine, critical_line, EquilibriumState, critical_point
   use YAEOS__MODELS_AR_GERG2008, only: Gerg2008, gerg_2008
   use yaeos__consistency_armodel, only: numeric_ar_derivatives
   use testing_aux, only: test_title, assert
   implicit none
   type(Gerg2008) :: model
   type(CubicEoS) :: cubic

   real(pr) :: n(2), v, t, n0(2)
   real(pr) :: ar, arv, arv2
   real(pr) :: art, art2, artv
   real(pr) :: arvn(2), artn(2), arn(2), arn2(2,2)
   real(pr) :: f1, f2, f3, f4, dx
   integer :: comps(2) = [1, 4]

   real(pr) :: Arnum, ArVnum, ArV2Num, ArTnum, ArT2num, ArTVnum, ArNnum(2), ArN2num(2,2)

   model = gerg_2008(comps)
   cubic = model%srk

   print *, test_title("GERG 2008")

   n = [0.5, 0.5]

   v = 1
   T = 150

   call model%residual_helmholtz(n, v, t, ar=ar, arv=arv, art=art, artv=artv, arv2=arv2, art2=art2, arn=arn, arvn=arvn, artn=artn, arn2=arn2)
   call assert(abs(ar - (-11.1819)) < 1e-4, "Ar Value from literature")

   call numeric_ar_derivatives(model, n, V, T, d_n=0.00001_pr, d_v=0.00001_pr, d_t=0.001_pr, &
      Ar=ArNum, ArV=ArVNum, ArV2=ArV2Num, ArT=ArTNum, ArT2=ArT2Num, ArTV=ArTVNum, ArN=ArNNum, ArN2=ArN2Num)

   call assert(abs(ArNum - Ar) < 1e-5, "Ar Value")
   call assert(abs(ArVnum - ArV) < 1e-5, "ArV Value")
   call assert(abs(ArV2Num - ArV2) < 1e-3, "ArV2 Value")
   call assert(abs(ArTnum - ArT) < 1e-5, "ArT Value")
   call assert(abs(ArT2num - ArT2) < 1e-5, "ArT2 Value")
   call assert(abs(ArTVnum - ArTV) < 1e-5, "ArTV Value")
   call assert(all(abs(ArNnum - ArN) < 1e-5), "ArN Value")
   call assert(all(abs(ArN2num - ArN2) < 1e-3), "ArN2 Value")


   call model%volume(n, 1.0_pr, T, f1, root_type="liquid")
   call cubic%volume(n, 1.0_pr, T, f2, root_type="liquid")
   call assert(abs(f1 - f2) < 1e-2, "Liquid root close to SRK")
   
   call model%volume(n, 1.0_pr, T, f1, root_type="vapor")
   call cubic%volume(n, 1.0_pr, T, f2, root_type="vapor")

   call assert(abs(f1 - f2) < 0.1, "Vapor root close to SRK")
   
   call model%volume(n, 1.0_pr, T, f1, root_type="stable")
   call cubic%volume(n, 1.0_pr, T, f2, root_type="stable")
   call assert(abs(f1 - f2) < 1e-1, "Stable root close to SRK")

end program main
