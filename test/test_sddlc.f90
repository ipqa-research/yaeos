program test_sddlc
   use yaeos, only: pr, sDDLC, CubicEoS, RKPR, EquilibriumState, saturation_pressure, px_envelope_2ph, PXEnvel2, R, numeric_ar_derivatives
   use auxiliar_functions, only: allclose
   use testing_aux, only: assert, test_title
   integer, parameter :: nc=2

   real(pr), dimension(nc) :: Tc, Pc, w, zc, q
   real(pr) :: kij(nc, nc), lij(nc, nc), Tstar(nc, nc)
   type(CubicEoS) :: model
   type(sDDLC) :: mr
   type(EquilibriumState) :: sat
   type(PXEnvel2) :: px
   real(pr) :: z0(nc)=[0, 1], zi(nc)=[1, 0], T=323, P, z(nc), V

   real(pr) :: as(4), Ps(4), a

   integer :: i

   write(*, *) test_title("sDDLC Mixing Rule")

   Ps = [57.143, 169.924, 242.857, 242.105]
   as = [0.263, 0.5033, 0.6004, 0.6004]

   Tc =  [190.564_pr, 562.05_pr]
   Pc =  [45.99_pr, 48.95_pr]
   w =  [0.0115478_pr, 0.2103_pr]
   zc =  [0.286_pr, 0.268_pr]

   q = [1.16_pr, 6.00_pr]


   kij = 0
   lij = 0
   Tstar(1, :) = [0, 190]
   Tstar(2, :) = [190, 0]

   mr = sDDLC(q=q, k=kij, k0=kij, Tref=Tstar,l=lij)
   model = RKPR(Tc, Pc, w, zc)
   call model%set_mixrule(mr)
   
   call test_numdiff

   ! do i=1,4
   !    a = as(i)
   !    z = a * zi + (1-a) * z0
   !    sat = saturation_pressure(model=model, n=z, T=323.15_pr, kind="bubble", p0=Ps(i))
   !    call assert(abs(sat%P - Ps(i)) < 2._pr, "Different Psat compared with Cismondi ")
   ! end do

   ! z = [0.01, 0.99]
   ! sat = saturation_pressure(model=model, n=z, T=323.15_pr, kind="bubble", p0=1._pr)
   ! px = px_envelope_2ph(model, z0, 0.01._pr, zi, sat)


contains

   subroutine test_numdiff
      real(pr) :: V, T, z(nc)
      real(pr) :: Ar_num, ArV_num, ArT_num, ArV2_num, ArT2_num, ArTV_num
      real(pr), dimension(nc) :: Arn_num, ArVn_num, ArTn_num
      real(pr), dimension(nc, nc) :: Arn2_num
      real(pr) :: Ar, ArV, ArT, ArV2, ArT2, ArTV
      real(pr), dimension(nc) :: Arn, ArVn, ArTn
      real(pr), dimension(nc, nc) :: Arn2

      real(pr), parameter :: dn = 1e-5, dv=1e-5, dt=1e-2
      real(pr), parameter :: tol=1e-5

      z = [0.5, 0.5]
      V = 5
      T = 350

      call numeric_ar_derivatives(model, z, V, T, dn, dv, dt, &
         Ar_num, ArV_num, ArT_num, Arn_num, ArV2_num, ArT2_num, ArTV_num, ArVn_num, ArTn_num, Arn2_num)

      call model%residual_helmholtz(z, V, T, &
         Ar=Ar, ArV=ArV, ArT=ArT, Arn=Arn, ArV2=ArV2, ArT2=ArT2, ArTV=ArTV, &
         ArVn=ArVn, ArTn=ArTn, Arn2=Arn2 &
         )


      call assert(allclose([Ar, Ar], tol), "Ar")
      call assert(allclose([ArT_num , ArT], tol), "ArV")
      call assert(allclose([ArT2_num, ArT2], tol), "ArV2")
   end subroutine test_numdiff
end program test_sddlc