module yaeos__models_external_apis_feos
   !! # `FeOs` API
   !! API to work with the `FeOs` rust library models.
   !!
   !! # Description
   !! Work in progress
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  A basic code example
   !! ```
   !!
   !! # References
   !! [FeOs](https://github.com/feos-org/feos)
   use yaeos__models_ar, only: ArModel
   use yaeos__constants, only: pr
   implicit none

   type, abstract, extends(ArModel) :: ArModelFeOs
      real(pr), allocatable :: parameters(:)
      !! Relevant model parameteres, should be specific from model
   contains
      procedure :: residual_helmholtz
   end type ArModelFeOs

   type, extends(ArModelFeOs) :: PCSAFT
      real(pr), allocatable :: m(:)
      real(pr), allocatable :: sigma(:)
      real(pr), allocatable :: epsilon_k(:)
      real(pr), allocatable :: kij(:, :)
   contains
      procedure :: get_v0 => pcsaft_v0
   end type PCSAFT

contains

   real(pr) function pcsaft_v0(self, n, P, T) result(v0)
      class(PCSAFT), intent(in) :: self
      real(pr), intent(in) :: n(:), P, T
      error stop 1
   end function pcsaft_v0

   function pc_saft_to_str(self) result(json_str)
      use json_module, only: json_core, json_value, json_array
      use json_kinds, only: CK, IK
      class(PCSAFT), intent(in) :: self
      character(kind=CK, len=:), allocatable :: json_str

      type(json_core) :: json
      type(json_value), pointer :: base, molecule, molecules, params, ident

      integer :: i

      call json%initialize()

      call json%create_object(base,'')
      call json%add(base, "residual_model", "PC-SAFT")
      call json%add(base, "ideal_gas_model", "")
      call json%create_array(molecules,'residual_substance_parameters')
      do i=1,size(self%m)
         call json%create_object(molecule, "")
         call json%create_object(params, 'model_record')
         call json%create_object(ident, 'identifier')

         call json%add(params, "m", self%m(i))
         call json%add(params, "sigma", self%sigma(i))
         call json%add(params, "epsilon_k", self%epsilon_k(i))

         call json%add(ident, "cas", "")
         call json%add(ident, "name", "")
         call json%add(ident, "iupac_name", "")
         call json%add(ident, "smiles", "")
         call json%add(ident, "inchi", "")
         call json%add(ident, "formula", "")

         call json%add(molecule, "molarweight", 0.0)
         call json%add(molecule, params)
         call json%add(molecule, ident)

         call json%add(molecules, molecule)
      end do
      
      call json%add(base, molecules)

      call json%print_to_string(base, json_str)
   end function pc_saft_to_str


!       "{ \
!     \"residual_model\": \"PC-SAFT\", \
!     \"ideal_gas_model\": \"\", \
!     \"residual_substance_parameters\": [ \
!         { \
!             \"identifier\": { \
!                 \"cas\": \"74-82-8\", \
!                 \"name\": \"methane\", \
!                 \"iupac_name\": \"methane\", \
!                 \"smiles\": \"C\", \
!                 \"inchi\": \"InChI=1/CH4/h1H4\", \
!                 \"formula\": \"CH4\" \
!             }, \
!             \"model_record\": { \
!                 \"m\": 1.0, \
!                 \"sigma\": 3.7039, \
!                 \"epsilon_k\": 150.03 \
!             }, \
!             \"molarweight\": 16.043 \
!         } \
!     ] \
! }";

   subroutine model_setter(self, setup_string)
      !! Setup the model from a string
      class(ArModelFeoS), intent(out) :: self
      character(len=*), intent(in) :: setup_string

      ! Setup logic from some input string?
      ! I'd prefer a more specific implementation but I'm not sure how easy
      ! to extend it could be
   end subroutine model_setter

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz model generic interface.
      !!
      !! This interface represents how an Ar model should be implemented.
      !! By our standard, a Resiudal Helmholtz model takes as input:
      !!
      !! - The mixture's number of moles vector.
      !! - Volume, by default in liters.
      !! - Temperature, by default in Kelvin.
      !!
      !! All the output arguments are optional. While this keeps a long
      !! signature for the implementation, this is done this way to take
      !! advantage of any inner optimizations to calculate derivatives
      !! inside the procedure.
      !!
      !! Once the model is implemented, the signature can be short like
      !! `model%residual_helmholtz(n, v, t, ArT2=dArdT2)`
      class(ArModelFeOs), intent(in) :: self !! ArModel
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: v !! Volume [L]
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(out) :: Ar !! Residual Helmoltz energy
      real(pr), optional, intent(out) :: ArV !! \(\frac{dAr}{dV}\)
      real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
      real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
      real(pr), optional, intent(out) :: ArTV !! \(\frac{d^2Ar}{dTV}\)
      real(pr), optional, intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
      real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
      real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVn_i}\)
      real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTn_i}\)
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))!! \(\frac{d^2Ar}{dn_{ij}}\)

      ! Here the function that calculates residual helmholtz energy
      ! from FeOs would be called.
   end subroutine residual_helmholtz
end module yaeos__models_external_apis_feos
