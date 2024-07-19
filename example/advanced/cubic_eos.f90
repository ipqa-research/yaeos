!> We have seen that some cubic models are available in yaeos. But will show
!> here how it is possible to make your own CubicEoS by choosing which piece
!> of the model you want to use.
!>
!> This is just an example model and it is not intended to be used in real
!> applications. It is just to show how to create a CubicEoS model.
program main
   use yaeos
   integer, parameter :: nc = 2 !! Number of components
   type(CubicEoS) :: eos !! The CubicEoS model
   
   real(pr) :: tc(nc), pc(nc), w(nc) !! Critical constants
   real(pr) :: kij(nc, nc), lij(nc, nc) !! Binary interaction parameters
   character(len=50) :: name(nc) !! Name of the components
   
   type(Substances) :: composition  
      !! All models should have their own composition
   type(AlphaRKPR) :: alpha_function  
      !! The RKPR Alpha Function
   type(QMR) :: mixrule
      !! A Quadratic mixing rule (ClassicVdW)

   integer :: i
   real(pr) :: n(nc), V, P, T

   ! First we define the components critical constants
   Tc = [190._pr,  310._pr]
   Pc = [14._pr,   30._pr ]
   w = [0.001_pr, 0.03_pr]
   name = ["Component 1", "Component 2"]
   composition = Substances(name, Tc, Pc, w)


   ! RKPR Alpha function uses a set of k parameters
   alpha_function%k = [0.48_pr, 0.58_pr]

   ! The mixrule
   kij = 0
   lij = 0
   kij(1,2) = 0.1_pr
   kij(2,1) = kij(1, 2)
   mixrule = QMR(k=kij, l=lij)

   ! We need to set up the model parameters
   eos%ac = [0.042748_pr, 0.052748_pr]
   eos%b  = [0.005, 0.001]
   eos%del1 = [0.1, 0.2]
   eos%del2 = [0.4, 0.5]

   ! We now add the before defined components and mixrule
   eos%alpha = alpha_function
   eos%components = composition
   eos%mixrule = mixrule


   T = 50
   n = [0.5, 0.5]
   do i=1,100
      V = real(i, pr)/100
      call eos%pressure(n, V, T, P)
      print *, V, P

   end do


end program