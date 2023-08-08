module yaeos_mixrule_classicvdw
   !-|  Mixing rules module
   !
   !    Cubic EoS Mixing Rule
   !    Contains a single subroutine as an attribute that receives the
   !    set of parameters of pure compounds and returns the mixture's
   !    parameters and their derivatives
   use yaeos_constants, only: pr
   use yaeos_autodiff
   
   implicit none

   real(pr), allocatable :: kij(:, :) !| \(k_{ij}\) BIP
   real(pr), allocatable :: lij(:, :) !| \(l_{ij}\) BIP

contains
   ! ===========================================================================
   ! Constructors
   ! ---------------------------------------------------------------------------
   subroutine setup_ClassicVdW(kij_in, lij_in)
      !-| Set up the module's parameters.
      !   this allocates and assign the values to the parameters to be used.
      real(pr) :: kij_in(:, :)
      real(pr) :: lij_in(:, :)

      kij = kij_in
      lij = lij_in

   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Mixing subroutines
   ! ---------------------------------------------------------------------------
   pure subroutine mix_ClassicVdW(z, v, t, a, b, c, amix, bmix, cmix)
      !-| Mix the parameters.
      !   This subroutine mix all the parameters for the generic Cubic EoS.
      type(hyperdual), intent(in)  :: z(:)       !| Composition
      type(hyperdual), intent(in)  :: v          !| Volume
      type(hyperdual), intent(in)  :: t          !| Temperature
      type(hyperdual), intent(in)  :: a(size(z)) !| Attractive parameter
      type(hyperdual), intent(in)  :: b(size(z)) !| Repulsive parameter
      type(hyperdual), intent(in)  :: c(size(z)) !| Volume traslation
      
      type(hyperdual), intent(out) :: amix       !| Mixture's attractive param
      type(hyperdual), intent(out) :: bmix       !| Mixture's repulsive param
      type(hyperdual), intent(out) :: cmix       !| Mixture's VT
      
      call a_mix_classic(z, v, t, a, amix)
      call b_mix_classic(z, v, t, b, bmix)
      call c_mix_classic(z, v, t, c, cmix)
   end subroutine

   pure subroutine a_mix_classic(z, v, t, a, amix)
      !| Calculate the mixture's attractive parameter.
      !
      ! \[a_{mix} = \sum_i z_i\sum_j z_j \sqrt{(a_{i} a_{j})} (1 - k_{ij})\]
      type(hyperdual), intent(in)  :: z(:) !| Composition
      type(hyperdual), intent(in)  :: v    !| Volume
      type(hyperdual), intent(in)  :: t    !| Temperature
      type(hyperdual), intent(in)  :: a(size(z)) !| Pure's attractive parameter
      type(hyperdual), intent(out) :: amix       !| Mixture's Attractive param

      type(hyperdual) :: ai(size(z)), z2(size(z)), zij

      integer :: i, j

      ai = sqrt(a)
      z2 = z * z

      do i=1,size(z)-1
         do j=i+1,size(z)
            zij = z(i) * z(j)
            amix = amix + zij * (ai(i) * ai(j)) * (1 - kij(i, j))
         end do
      end do
      amix = 2 * amix + sum(z2 * a)
   end subroutine

   pure subroutine b_mix_classic(z, v, t, b, bmix)
      !-| Calculate the mixture's repulsive parameter.
      !
      ! \[ b_{mix} = \sum_i z_i \sum_j z_j (b_i + b_j)/2 (1 - l_{ij}) \]
      !
      type(hyperdual), intent(in)  :: z(:)       !| Composition
      type(hyperdual), intent(in)  :: v          !| Volume
      type(hyperdual), intent(in)  :: t          !| Temperature
      type(hyperdual), intent(in)  :: b(size(z)) !| Pure's repulsive parameter
      type(hyperdual), intent(out) :: bmix       !| Mixture's repulsive param

      integer :: i, j

      bmix = 0.0_pr
      
      do i=1,size(z)
         do j=i+1, size(z)
            bmix = bmix + z(i) * z(j) * (b(i) + b(j)) * (1 - lij(i, j))
         end do
      end do
      bmix = bmix + sum(z * z * b)

      bmix = bmix/sum(z)
   end subroutine

   pure subroutine c_mix_classic(z, v, t, c, cmix)
      !-| Calculate the mixture's volume traslation parameter.
      !
      ! \[ c_{mix} = \sum_i z_i c_i \]
      !
      type(hyperdual), intent(in)  :: z(:)       !| Composition
      type(hyperdual), intent(in)  :: v          !| Volume
      type(hyperdual), intent(in)  :: t          !| Temperature
      type(hyperdual), intent(in)  :: c(size(z)) !| Pure's VT parameter
      type(hyperdual), intent(out) :: cmix       !| Mixture's VT param

      cmix = sum(z*c)
   end subroutine
   
   pure subroutine del_mix_classic(z, v, t, del, delmix)
      !-| For a classic Cubic EoS \(\delta_i\) parameter is assumed as the first
      !   component parameter (since it's the same constant for all components)
      type(hyperdual), intent(in)  :: z(:)         !| Composition
      type(hyperdual), intent(in)  :: v            !| Volume
      type(hyperdual), intent(in)  :: t            !| Temperature
      type(hyperdual), intent(in)  :: del(size(z)) !| Pure's VT parameter
      type(hyperdual), intent(out) :: delmix       !| Mixture's VT param

      delmix = del(1) ! sum(z * del)
   end subroutine
   ! ===========================================================================
end module
