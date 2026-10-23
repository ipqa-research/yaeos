program main
   use yaeos
   implicit none
   integer, parameter :: nv = 6
   type(GERG2008) :: eos
   integer :: i, j
   real(pr) :: P, T, V, Gr
   real(pr) :: MM = 14*2

   real(pr) :: Vs(nv)

   eos = get_modelgerg()

   T = 130

   V = 0.1

   vs = [ 0.035001, .040, 0.046666, 0.056, 0.093336, 0.1399972 ]
   print *, "P V T Gr"
   do j=1, nv
      V = vs(j)! /MM * 1e-3_pr
      do i=50, 1000, 50
         T = real(i, pr)
         call eos%pressure([1._pr], V=V, T=T, P=P)
         call eos%gibbs_residual_vt(n=[1._pr], V=V, T=T, Gr=Gr)
         Gr =  Gr - R * T * log(P*V/(R*T))
         print *, P, V, T, Gr
      end do
   end do

contains
   type(GERG2008) function get_modelgerg() result(model)
      integer :: ids(1)
      ! ids = [G2008Components%methane
      ids = [G2008Components%nitrogen]
      model = gerg_2008(ids)
   end function get_modelgerg
end program main
