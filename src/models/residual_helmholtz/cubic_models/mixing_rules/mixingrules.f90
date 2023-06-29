module mixrule_classicvdw
   !->  Mixing rules module
   !
   !    Cubic EoS Mixing Rule
   !    Contains a single subroutine as an attribute that receives the
   !    set of parameters of pure compounds and returns the mixture's
   !    parameters and their derivatives
   use constants, only: pr
   use hyperdual_mod
   
   implicit none

   real(pr), allocatable :: kij(:, :) !- \(k_{ij}\) BIP
   real(pr), allocatable :: lij(:, :) !- \(l_{ij}\) BIP

contains
   ! ===========================================================================
   ! Constructors
   ! ---------------------------------------------------------------------------
   subroutine setup_ClassicVdW(kij_in, lij_in)
      !-> Set up the module's parameters.
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
      !-> Mix the parameters.
      !   This subroutine mix all the parameters for the generic Cubic EoS.
      type(hyperdual), intent(in)  :: z(:)       !- Composition
      type(hyperdual), intent(in)  :: v          !- Volume
      type(hyperdual), intent(in)  :: t          !- Temperature
      type(hyperdual), intent(in)  :: a(size(z)) !- Attractive parameter
      type(hyperdual), intent(in)  :: b(size(z)) !- Repulsive parameter
      type(hyperdual), intent(in)  :: c(size(z)) !- Volume traslation
      
      type(hyperdual), intent(out) :: amix       !- Mixture's attractive param
      type(hyperdual), intent(out) :: bmix       !- Mixture's repulsive param
      type(hyperdual), intent(out) :: cmix       !- Mixture's VT
      
      call a_mix_classic(z, v, t, a, amix)
      call b_mix_classic(z, v, t, b, bmix)
      call c_mix_classic(z, v, t, c, cmix)
   end subroutine

   pure subroutine a_mix_classic(z, v, t, a, amix)
      !- Calculate the mixture's attractive parameter.
      !
      ! \[a_{mix} = \sum_i z_i\sum_j z_j \sqrt{(a_{i} a_{j})} (1 - k_{ij})\]
      type(hyperdual), intent(in)  :: z(:) !- Composition
      type(hyperdual), intent(in)  :: v    !- Volume
      type(hyperdual), intent(in)  :: t    !- Temperature
      type(hyperdual), intent(in)  :: a(size(z)) !- Pure's attractive parameter
      type(hyperdual), intent(out) :: amix       !- Mixture's Attractive param

      type(hyperdual) :: aij(size(z), size(z))

      integer :: i, j

      aij = 0.0_pr
      amix = 0.0_pr
      do j = 1, size(z)
         aij(:, j) = sqrt(a(:) * a(j)) * (1.0_pr - kij(:, j))
         amix = amix + sum(z(:) * z(j) * aij(:, j))
      end do
   end subroutine

   pure subroutine b_mix_classic(z, v, t, b, bmix)
      !-> Calculate the mixture's repulsive parameter.
      !
      ! \[ b_{mix} = \sum_i z_i \sum_j z_j (b_i + b_j)/2 (1 - l_{ij}) \]
      !
      type(hyperdual), intent(in)  :: z(:)       !- Composition
      type(hyperdual), intent(in)  :: v          !- Volume
      type(hyperdual), intent(in)  :: t          !- Temperature
      type(hyperdual), intent(in)  :: b(size(z)) !- Pure's repulsive parameter
      type(hyperdual), intent(out) :: bmix       !- Mixture's repulsive param

      type(hyperdual) :: bij(size(z), size(z))

      integer :: i

      bmix = 0.0_pr
      do i = 1, size(z)
         bij(:, i) = (b(:) + b(i))/2.0_pr*(1.0_pr - lij(:, i))
         bmix = bmix + sum(z(:)*z(i)*bij(:, i))
      end do

      bmix = bmix/sum(z)
   end subroutine

   pure subroutine c_mix_classic(z, v, t, c, cmix)
      !-> Calculate the mixture's volume traslation parameter.
      !
      ! \[ c_{mix} = \sum_i z_i c_i \]
      !
      type(hyperdual), intent(in)  :: z(:)       !- Composition
      type(hyperdual), intent(in)  :: v          !- Volume
      type(hyperdual), intent(in)  :: t          !- Temperature
      type(hyperdual), intent(in)  :: c(size(z)) !- Pure's VT parameter
      type(hyperdual), intent(out) :: cmix       !- Mixture's VT param

      cmix = sum(z*c)
   end subroutine
   
   pure subroutine del_mix_classic(z, v, t, del, delmix)
      !-> For a classic Cubic EoS \(\delta_i\) parameter is assumed as the first
      !   component parameter (since it's the same constant for all components)
      type(hyperdual), intent(in)  :: z(:)         !- Composition
      type(hyperdual), intent(in)  :: v            !- Volume
      type(hyperdual), intent(in)  :: t            !- Temperature
      type(hyperdual), intent(in)  :: del(size(z)) !- Pure's VT parameter
      type(hyperdual), intent(out) :: delmix       !- Mixture's VT param

      delmix = del(1)
   end subroutine
   ! ===========================================================================
end module