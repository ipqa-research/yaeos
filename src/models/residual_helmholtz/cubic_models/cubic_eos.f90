module yaeos_generic_cubic
   !-|
   !# Generic Cubic Equation of State.
   !
   !   This module implements a generic cubic based on the definition given
   !   by Michelsen and Møllerup:
   !
   !   \[ 
   !    \alpha^r(T, V, n) = - n \ln(1 - B/V) 
   !                        - \frac{D(T)}{RTB(\delta_1 - \delta_2)} \ln\left(\frac{1+\delta_1 B/V}{1+\delta_2 B/V}\right)
   !   \]
   !
   !   Which considers the parameters:
   !   
   !   - \(D\): Attractive parameter
   !   - \(B\): Repulsive parameter
   !   - \(\delta_i\): EoS parameters
   !
   !   All the parameters are defined as internal procedures via pointers, which
   !   should point at the corresponding function. There is no need to understand
   !   the internal work of Fortran pointers, just definining the important
   !   functions and using the `set_functions` subroutine is enough to make
   !   this method to work. For example, adding the Van der Waals Cubic EoS 
   !   from scratch would be something like:
   ! 
   !
   ! @note First make a module that use the relevant data types, like 
   ! real precision and dual numbers to use automatic differentiation and also
   ! include the setter functions
   ! 
   !```fortran
   ! module cubic_van_der_waals
   !    use yaeos_constants, only: pr ! Import the real precision
   !    use yaeos_autodiff ! Use the automatic differentiation types
   ! 
   !    ! From the generic cubic eos use the `set_functions` subroutine and the
   !    ! Ar function one.
   !    use yaeos_generic_cubic, only: set_functions, a_res 
   ! 
   !    ! We define the model parameters to be used
   !    ! Model parameters: 
   !    type(hyperdual), allocatable :: a(:) ! Attractive
   !    type(hyperdual), allocatable :: b(:) ! Repulsive
   !```
   !
   ! @note Then, inside the module (or not), define all the relevant subroutines
   ! for parameters calculation and their mixing rule. In this case we define
   ! the Cubic Van der Waals EoS, where all parameters are held constant, and
   ! \(\delta_i=0\)
   !
   !```fortran
   !    ! Define how each required parameter is calculated
   !    subroutine a_parameter(z, v, t, a_out)
   !       type(hyperdual), intent(in) :: z(:), v, t
   !       type(hyperdual) :: a_out(size(z))
   !       ! In the VdW EoS the a parameter is held constant
   !       a_out = a ! Use the module's defined a parameter
   !    end subroutine
   !    
   !    subroutine b_parameter(z, v, t, b_out)
   !       type(hyperdual), intent(in) :: z(:), v, t
   !       type(hyperdual) :: b_out(size(z))
   !       ! In the VdW EoS the b parameter is held constant
   !       b_out = b ! Use the module's defined a parameter
   !    end subroutine
   !    
   !    subroutine no_parameter(z, v, t, no)
   !       type(hyperdual), intent(in) :: z(:), v, t
   !       type(hyperdual) :: no(size(z))
   !       ! We'll assign zero to non used parameters, like volume traslation
   !       ! or delta_1 and delta_2
   !       no = 0.0_pr
   !    end subroutine
   ! 
   !    ! And the mixing rule
   !    subroutine mixrule(z, v, t, a_p, b_p, c_p, amix, bmix, cmix)
   !       type(hyperdual), intent(in) :: z(:), v, z, a_p(size(z)), b_p(size(z)), c_p(size(z))
   !       type(hyperdual), intent(out) :: amix, bmix, cmix
   !       ! For simplicity we'll assume linear mixing
   ! 
   !       amix = sum(z * a_p)
   !       bmix = sum(z * b_p)
   !       cmix = 0.0_pr ! No VT
   !    end subroutine
   !    ! ==========================================================================
   !```
   !
   ! @note Now having all the relevant procedures defined, we add a simple
   ! `setup` subroutine to receive from somewhere (an input interpetator routine
   ! for example) the number of components and the values of their parameters.
   ! Inside this subroutine the Generic CEoS are setted up, and finally the 
   ! Generic Cubic EoS Ar subroutine is selected
   !
   !```fortran
   !
   !    ! Finally setup the generic cubic model:
   !    subroutine setup(n, a_in, b_in)
   !       use yaeos_generic_cubic, only: set_functions
   !       use yaeos_ar_models, only: set_ar_function
   !       integer, intent(in) :: n
   !       real(pr), intent(in) :: a_in(n), b_in(n)
   ! 
   !       a = a_in; b = b_in ! Setup the model's parameters
   !       
   !       call set_functions(& ! Setup the generic cubic functions
   !          a_parameter, b_parameter, & ! a and b subroutines
   !          no_parameter, no_parameter, no_parameter, &  ! null delta1, delta2 and c
   !          mixrule) ! Mixing rule
   !       call set_ar_function(a_res) ! Use the generic cubic eos subroutine as the
   !                                   ! main Ar subroutine.
   !    end subroutine
   !end module
   !```
   
   use yaeos_constants, only: pr, R
   use yaeos_autodiff
   use yaeos_interfaces, only: pures_property, abs_cubic_mix
   use yaeos_ar_models, only: set_ar_function
   use yaeos_thermo_properties, only: vinit

   implicit none

   private

   public :: set_functions !| Setter of functions to use parameters
   public :: a_res !| Residual free Helmholtz energy function

   procedure(pures_property), pointer :: attractive_parameter !| Attractive parameter function
   procedure(pures_property), pointer :: repulsive_parameter !| Repulsive parameter functoin
   procedure(pures_property), pointer :: delta1_parameter !| \(\delta_1\) parameter
   procedure(pures_property), pointer :: delta2_parameter !| \(\delta_2\) parameter
   procedure(pures_property), pointer :: volume_traslation !| Volume traslation parameter
   procedure(abs_cubic_mix), pointer :: mix !| Mixing function

contains

   subroutine set_functions(a, b, c, d1, d2, mixfun)
      !-| Setter of Generic Cubic EoS prcoedures.
      !
      !   This subroutine receives the desired subroutines to be
      !   used in the calculation of each parameter, and assings the 
      !   corresponding pointers to them. As well as the mixing rule
      !   subroutine to be used and assigns volume initializer function.
      procedure(pures_property) :: a !| Attractive parameter
      procedure(pures_property) :: b !| Repulsive parameter
      procedure(pures_property) :: c !| Volume traslation parameter
      procedure(pures_property) :: d1 !| \(\delta_1\) parameter
      procedure(pures_property) :: d2 !| \(\delta_2\) parameter
      procedure(abs_cubic_mix) :: mixfun !| Mixing rule

      attractive_parameter => a
      repulsive_parameter => b
      delta1_parameter => d1
      delta2_parameter => d2
      volume_traslation => c
      mix => mixfun

      vinit => v0
   end subroutine
      
   pure subroutine a_res(z, v, t, ar)
      !-| Generic Cubic EoS Ar function
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v !| Volume
      type(hyperdual), intent(in) :: t !| Temperature
      type(hyperdual), intent(out) :: ar !| Residual Helmholtz energy

      type(hyperdual), dimension(size(z)) :: a_pures, b_pures, c_pures
      type(hyperdual), dimension(size(z)) ::  del1_pures, del2_pures
      type(hyperdual) :: a, b, c, del1, del2

      type(hyperdual) :: b_v

      ! Pure components parameters
      call attractive_parameter(z, v, t, a_pures)
      call repulsive_parameter(z, v, t, b_pures)

      call delta1_parameter(z, v, t, del1_pures)
      call delta2_parameter(z, v, t, del2_pures)

      ! Call the model's mixing rule
      ! TODO: Update it to also mix delta 1 and 2
      call mix( &
               z, v, t, &
               a_pures, b_pures, c_pures, &
               a, b, c &
      )

      del1 = del1_pures(1)
      del2 = del2_pures(1)

      b_v = b/v

      ar = (&
            -sum(z) * log(1.0_pr - b_v) &
            - a/(R*t*b)*1.0_pr/(del1 - del2) & 
            * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
      ) ! * R * t
   end subroutine

   pure function v0(z, p, t)
      !-| Volume initializer function.
      !
      ! The classic cubic equation uses the mixture's covolume as
      ! an initializer
      real(pr), intent(in) :: z(:) !| Compositon
      real(pr), intent(in) :: p !| Pressure
      real(pr), intent(in) :: t !| Temperature
      real(pr) :: v0 !| Initial volume

      type(hyperdual) :: z_d(size(z)), p_d, t_d, b(size(z))

      z_d = z
      p_d = p
      t_d = t

      call repulsive_parameter(z_d, p_d, t_d, b)
      v0 = sum(z*b%f0)*1.2_pr
   end function
end module

module yaeos_cubic_eos
   !-| Set of subroutine that are part of the classic Cubic EoS, like
   ! PR and SRK alpha function and constant rest of parameters.
   use yaeos_constants, only: pr, R
   use yaeos_autodiff

   implicit none

   ! Model parameters
   real(pr), allocatable :: ac(:) !| Critical atractive parameter)
   real(pr), allocatable :: b(:)  !| Repulsive parameter)
   real(pr), allocatable :: c(:)  !| Volume traslation)
   real(pr), allocatable :: k(:)  !| alpha k-parameter)

   real(pr), allocatable :: del1(:) !| \(\delta_1\) Parameter)
   real(pr), allocatable :: del2(:) !| \(\delta_2\) Parameter)

   ! Critical constants)
   real(pr), allocatable :: pc(:) !| Critical Pressure)
   real(pr), allocatable :: tc(:) !| Critical Temperature)
   real(pr), allocatable :: w(:)  !| Acentric factor)

contains

   subroutine set_parameters(&
      ac_in, b_in, c_in, k_in, del1_in, del2_in, &
      pc_in, tc_in, w_in & 
   )
      real(pr), intent(in) :: ac_in(:) !| Critical atractive parameter
      real(pr), intent(in) :: b_in(size(ac_in))  !| Repulsive parameter
      real(pr), intent(in) :: c_in(size(ac_in))  !| Volume traslation
      real(pr), intent(in) :: k_in(size(ac_in))  !| alpha k-parameter

      real(pr), intent(in) :: del1_in(size(ac_in)) !| \(\delta_1\) Parameter
      real(pr), intent(in) :: del2_in(size(ac_in)) !| \(\delta_2\) Parameter

      real(pr), intent(in) :: pc_in(size(ac_in)) !| Critical Pressure
      real(pr), intent(in) :: tc_in(size(ac_in)) !| Critical Temperature
      real(pr), intent(in) :: w_in(size(ac_in))  !| Acentric factor

      ac = ac_in
      b = b_in
      c = c_in
      k = k_in

      del1 = del1_in
      del2 = del2_in

      pc = pc_in
      tc = tc_in
      w = w_in
   end subroutine

   ! Allocators
   subroutine alloc(n)
      !-| Allocate the module's parameters to the desired number of components
      integer, intent(in) :: n !| Number of components
      call destroy()
      allocate(ac(n), b(n), c(n), k(n), del1(n), del2(n), pc(n), tc(n), w(n))
   end subroutine
   
   subroutine destroy()
      !-| Deallocate all the module's parameters
      if (allocated(ac)) deallocate(ac)
      if (allocated(b)) deallocate(b)
      if (allocated(c)) deallocate(c)
      if (allocated(k)) deallocate(k)
      if (allocated(del1)) deallocate(del1)
      if (allocated(del2)) deallocate(del2)
      if (allocated(pc)) deallocate(pc)
      if (allocated(tc)) deallocate(tc)
      if (allocated(w)) deallocate(w)
   end subroutine

   pure subroutine a_classic(z, v, t, a_out)
      !-| Atractive parameter.
      !
      !  We assume as classic the PR/SRK equation for the alpha function:
      !  VdW EoS could be a special case where \(k=0\)
      !  \[a = a_c \left(1 + k \left(1 - \sqrt{\frac{T}{T_c}}\right)\right)^2\]
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v    !| Volume
      type(hyperdual), intent(in) :: t    !| Temperature
      type(hyperdual), intent(out) :: a_out(size(z)) !| Attractive parameter

      a_out = ac*(1.0_pr + k*(1.0_pr - sqrt(t/tc)))**2
   end subroutine

   pure subroutine b_classic(z, v, t, b_out)
      !-| Repulsive parameter
      !
      ! The repulsive parameter is held constant.
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v    !| Volume
      type(hyperdual), intent(in) :: t    !| Temperature
      type(hyperdual), intent(out) :: b_out(size(z)) !| Repulsive parameter
      b_out = b
   end subroutine
   
   pure subroutine c_classic(z, v, t, c_out)
      !-| Volume traslation
      !
      !   Most of volume traslation methods keep the parameter constant
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v    !| Volume
      type(hyperdual), intent(in) :: t    !| Temperature
      type(hyperdual), intent(out) :: c_out(size(z)) !| c (VT) parameter
      c_out = c
   end subroutine

   pure subroutine del1_classic(z, v, t, del1_out)
      !-| \(\delta_1\) Parameter.
      !
      !   In most cubic equation of state systems this parameter is
      !   held constant.
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v    !| Volume
      type(hyperdual), intent(in) :: t    !| Temperature
      type(hyperdual), intent(out) :: del1_out(size(z)) !| \(\delta_1\)
      del1_out = del1
   end subroutine

   pure subroutine del2_classic(z, v, t, del2_out)
      !-| \(\delta_2\) Parameter.
      !
      !   In most cubic equation of state systems this parameter is
      !   held constant.
      type(hyperdual), intent(in) :: z(:) !| Composition
      type(hyperdual), intent(in) :: v    !| Volume
      type(hyperdual), intent(in) :: t    !| Temperature
      type(hyperdual), intent(out) :: del2_out(size(z)) !| \(\delta_2\)
      del2_out = del2
   end subroutine

end module
