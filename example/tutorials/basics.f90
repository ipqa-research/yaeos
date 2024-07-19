program basics
   use yaeos ! First the module is imported (used). This will bring all the
             ! types and functions defined in the module to this program.
   implicit none
   ! ===========================================================================
   
   ! In yaeos, all residual Helmholtz models are based on the structure
   !  `ArModel` and are derivatives from it.
   ! Defining the model variable as an allocatable ArModel means that this
   ! variable can be redefined as any other kind of variable that derives from
   ! `ArModel`, for example `CubicEoS`.
   class(ArModel), allocatable :: model

   integer, parameter :: nc=2 !! In this case we define a constant number of
                              !! components.

   ! Used input variables
   real(pr) :: n(nc), Tc(nc), Pc(nc), w(nc)

   ! Output variables
   real(pr) :: P, dPdV
   ! ---------------------------------------------------------------------------

   ! Here the executed code starts

   ! A classic Cubic model can be defined with each component critial constants
   tc = [190._pr,  310._pr]
   pc = [14._pr,   30._pr ]
   w  = [0.001_pr, 0.03_pr]

   ! Here we chose to model with the PengRobinson EoS
   model = PengRobinson78(tc, pc, w)

   ! Number of moles vector
   n = [0.3, 0.7]

   ! Pressure calculation
   call model%pressure(n, V=2.5_pr, T=150._pr, P=P)
   print *, "P: ", P

   ! Derivatives can also be calculated when included as optional arguments!
   call model%pressure(n, V=2.5_pr, T=150._pr, P=P, dPdV=dPdV)
   print *, "dPdV: ", dPdV

   thermoproperties: block
      real(pr) :: Cpr, Cvr
      real(pr) :: Gr, GrT, GrV, Grn(nc)
      real(pr) :: Hr, HrT, HrV, Hrn(nc)
      real(pr) :: lnPhi(nc), dlnPhidT(nc), dlnPhidP(nc), dlnPhidn(nc, nc)

      real(pr) :: n(nc), P, V, T

      n = [0.3_pr, 0.7_pr]
      V = 3.5_pr
      P = 5.0_pr
      call model%Cp_residual_vt(n, V=2.5_pr, T=150._pr, Cp=Cpr)
      call model%Cv_residual_vt(n, V=2.5_pr, T=150._pr, Cv=Cvr)
      call model%gibbs_residual_vt(n, V=2.5_pr, T=150._pr, Gr=Gr, GrT=GrT, GrV=GrV, Grn=Grn)
      call model%enthalpy_residual_vt(n, V=2.5_pr, T=150._pr, Hr=Hr, HrT=HrT, HrV=HrV, Hrn=Hrn)
      call model%lnphi_pt(&
         n, P=P, T=T, root_type="stable", &
         lnPhi=lnPhi, dlnPhidT=dlnphidT, dlnPhidP=dlnPhidP, dlnPhidn=dlnPhidn &
      )
      call model%lnphi_vt(&
         n, V=V, T=T, &
         lnPhi=lnPhi, dlnPhidT=dlnphidT, dlnPhidP=dlnPhidP, dlnPhidn=dlnPhidn &
      )
      call model%volume(n, P, T, root_type="vapor", V=V)
   end block thermoproperties
end program basics