module hyperdual_mod
   !! # Hyperdual arithmetic
   !! Type and overloaded operations for exact first- and second-order directional derivatives.
   !!
   !! # Description
   !! A hyperdual value stores a primal value and two independent tangent seeds:
   !! \[x=x_0+x_1\varepsilon_1+x_2\varepsilon_2+x_{12}\varepsilon_1\varepsilon_2.\]
   !! The nilpotent rules \(\varepsilon_1^2=\varepsilon_2^2=0\) make the
   !! Taylor expansion exact through the mixed second derivative. Overloaded
   !! arithmetic and elementary functions propagate these coefficients, so
   !! derivatives are obtained without a finite-difference step-size error.
   !! For a scalar function the propagation rule is
   !! \[f(x)=f(x_0)+x_1f'(x_0)\varepsilon_1+x_2f'(x_0)\varepsilon_2
   !!   +\left(x_{12}f'(x_0)+x_1x_2f''(x_0)\right)\varepsilon_1\varepsilon_2.\]
   !! Seeding \(x_1=x_2=1\), \(x_{12}=0\) gives \(f'\) in `f1`/`f2` and
   !! \(f''\) in `f12`; seeding different variables in \(\varepsilon_1\) and
   !! \(\varepsilon_2\) gives mixed second derivatives.
   !!
   !! # Examples
   !! ```fortran
   !! type(hyperdual) :: x, y
   !! x = hyperdual(f0=2.0_pr, f1=1.0_pr, f2=1.0_pr)
   !! y = exp(x)
   !! print *, y%f0, y%f1, y%f2, y%f12
   !! ```
   !!
   !! # References
   !! Fike, J. A., and A. C. Alonso (2011). The Development of Hyper-Dual
   !! Numbers for Exact Second-Derivative Calculations. AIAA 49th Aerospace
   !! Sciences Meeting. DOI: 10.2514/6.2011-886.
   !! Rehner, P., and G. Bauer (2021). Application of Generalized (Hyper-) Dual
   !! Numbers in Equation of State Modeling. Frontiers in Chemical Engineering.
   !! DOI: 10.3389/fceng.2021.758090.
   !
   ! Original code provided by Philipp Rehner and Gernot Bauer,
   ! Institute of Thermodynamics and Thermal Process Engineering (ITT),
   ! University of Stuttgart, Stuttgart, Germany
   use yaeos__constants, only: pr
   implicit none

   type, bind(c) :: hyperdual
      !! # hyperdual
      !! Scalar number carrying a primal value and two derivative seeds.
      !!
      !! # Description
      !! The fields store coefficients of 1, epsilon1, epsilon2, and their
      !! product. When `f1` and `f2` are seeded with unit directions, `f12`
      !! stores the corresponding mixed second derivative.
      !!
      !! # Examples
      !! ```fortran
      !! x = hyperdual(f0=x0, f1=1.0_pr, f2=1.0_pr)
      !! y = exp(x)
      !! print *, y%f12
      !! ```
      real(pr) :: f0 = 0 !! Primal value.
      real(pr) :: f1 = 0 !! Coefficient of the first tangent direction.
      real(pr) :: f2 = 0 !! Coefficient of the second tangent direction.
      real(pr) :: f12 = 0 !! Coefficient of the mixed second derivative.
   end type hyperdual

   !---------------------------------------------------------------------
   !--- Operator interfaces ---------------------------------------------
   !---------------------------------------------------------------------

   ! Equal assignment
   !> # assignment (=)
   !> Assigns a scalar or hyperdual value while preserving the four-component representation.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = x
   !> ```
   interface assignment (=)
      procedure EqualHyperDualHyperDual
      procedure EqualHyperDualReal
   end interface

   ! Unary operator +
   !> # operator (+)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x + y
   !> ```
   interface operator (+)
      procedure PlusHyperDualHyperDual
   end interface

   ! Addition operator
   !> # operator (+)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x + y
   !> ```
   interface operator (+)
      procedure AddHyperDualHyperDual
      procedure AddHyperDualReal
      procedure AddRealHyperDual
   end interface

   ! Unary operator -
   !> # operator (-)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x - y
   !> ```
   interface operator (-)
      procedure MinusHyperDualHyperDual
   end interface

   ! Subtraction operator
   !> # operator (-)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x - y
   !> ```
   interface operator (-)
      procedure SubtractHyperDualHyperDual
      procedure SubtractHyperDualReal
      procedure SubtractRealHyperDual
   end interface

   ! Multiplication operator
   !> # operator (*)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   interface operator (*)
      procedure MultiplyHyperDualHyperDual
      procedure MultiplyHyperDualReal
      procedure MultiplyRealHyperDual
      procedure MultiplyHyperDualInt
      procedure MultiplyIntHyperDual
   end interface

   ! Division operator
   !> # operator (/)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x / y
   !> ```
   interface operator (/)
      procedure DivideHyperDualHyperDual
      procedure DivideHyperDualReal
      procedure DivideRealHyperDual
   end interface

   ! Power operator
   !> # operator (**)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x ** y
   !> ```
   interface operator (**)
      procedure PowerHyperDualInt
      procedure PowerHyperDualHyperDual
      procedure PowerHyperDualReal
   end interface



   !---------------------------------------------------------------------
   !--- Summation interface ---------------------------------------------
   !---------------------------------------------------------------------
   !> # sum
   !> Overloads the Fortran `sum` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = sum(x)
   !> ```
   interface sum
      module procedure SumHyperDual
      module procedure SumHyperDual2
   end interface sum



   !---------------------------------------------------------------------
   !--- Logical operator interfaces -------------------------------------
   !---------------------------------------------------------------------

   ! Equal operator.
   !> # operator (.eq.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .eq. y
   !> ```
   interface operator (.eq.) ! or (==)
      procedure eq_dd
      procedure eq_dr
      procedure eq_rd
      procedure eq_di
      procedure eq_id
   end interface

   ! Not equal operator.
   !> # operator (.ne.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .ne. y
   !> ```
   interface operator (.ne.) ! or (/=)
      procedure ne_dd
      procedure ne_dr
      procedure ne_rd
      procedure ne_di
      procedure ne_id
   end interface

   ! Less than operator.
   !> # operator (.lt.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .lt. y
   !> ```
   interface operator (.lt.) ! or (<)
      procedure lt_dd
      procedure lt_dr
      procedure lt_rd
      procedure lt_di
      procedure lt_id
   end interface

   ! Less than or equal operator.
   !> # operator (.le.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .le. y
   !> ```
   interface operator (.le.) ! or (<=)
      procedure le_dd
      procedure le_dr
      procedure le_rd
      procedure le_di
      procedure le_id
   end interface

   ! Greater than operator.
   !> # operator (.gt.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .gt. y
   !> ```
   interface operator (.gt.) ! or (>)
      procedure gt_dd
      procedure gt_dr
      procedure gt_rd
      procedure gt_di
      procedure gt_id
   end interface

   ! Greater than or equal operator.
   !> # operator (.ge.)
   !> Applies this arithmetic operation to hyperdual values and supported real/integer operands, propagating derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = x .ge. y
   !> ```
   interface operator (.ge.) ! or (>=)
      procedure ge_dd
      procedure ge_dr
      procedure ge_rd
      procedure ge_di
      procedure ge_id
   end interface



   !---------------------------------------------------------------------
   !--- Math function interfaces ----------------------------------------
   !---------------------------------------------------------------------

   ! Absolute value function
   !> # abs
   !> Overloads the Fortran `abs` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = abs(x)
   !> ```
   interface abs
      module procedure absHyperDual
   end interface

   ! Integer function
   !> # int
   !> Overloads the Fortran `int` intrinsic; returns the conversion of the primal value `f0` and discards derivative information.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = int(x)
   !> ```
   interface int
      module procedure intHyperDual
   end interface

   ! Nearest integer function
   !> # nint
   !> Overloads the Fortran `nint` intrinsic; returns the conversion of the primal value `f0` and discards derivative information.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = nint(x)
   !> ```
   interface nint
      module procedure nintHyperDual
   end interface

   ! Real function
   !> # real
   !> Overloads the Fortran `real` intrinsic; returns the conversion of the primal value `f0` and discards derivative information.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = real(x)
   !> ```
   interface real
      module procedure realHyperDual
   end interface

   ! Sign function
   !> # sign
   !> Overloads the Fortran `sign` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = sign(x, s)
   !> ```
   interface sign
      module procedure sign_dd
      module procedure sign_dr
      module procedure sign_rd
   end interface

   ! Sine function
   !> # sin
   !> Overloads the Fortran `sin` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = sin(x)
   !> ```
   interface sin
      module procedure sinHyperDual
   end interface

   ! Cosine function
   !> # cos
   !> Overloads the Fortran `cos` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = cos(x)
   !> ```
   interface cos
      module procedure cosHyperDual
   end interface

   ! Tangent function
   !> # tan
   !> Overloads the Fortran `tan` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = tan(x)
   !> ```
   interface tan
      module procedure tanHyperDual
   end interface

   ! Sqrt function
   !> # sqrt
   !> Overloads the Fortran `sqrt` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = sqrt(x)
   !> ```
   interface sqrt
      module procedure sqrtHyperDual
   end interface

   ! Log function
   !> # log
   !> Overloads the Fortran `log` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = log(x)
   !> ```
   interface log
      module procedure logHyperDual
   end interface

   ! Log10 function
   !> # log10
   !> Overloads the Fortran `log10` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = log10(x)
   !> ```
   interface log10
      module procedure log10HyperDual
   end interface

   ! Exp function
   !> # exp
   !> Overloads the Fortran `exp` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = exp(x)
   !> ```
   interface exp
      module procedure expHyperDual
   end interface

   ! Sinh function
   !> # sinh
   !> Overloads the Fortran `sinh` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = sinh(x)
   !> ```
   interface sinh
      module procedure sinhHyperDual
   end interface

   ! Cosh function
   !> # cosh
   !> Overloads the Fortran `cosh` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = cosh(x)
   !> ```
   interface cosh
      module procedure coshHyperDual
   end interface

   ! Tanh function
   !> # tanh
   !> Overloads the Fortran `tanh` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = tanh(x)
   !> ```
   interface tanh
      module procedure tanhHyperDual
   end interface

   ! Acos function
   !> # acos
   !> Overloads the Fortran `acos` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = acos(x)
   !> ```
   interface acos
      module procedure acosHyperDual
   end interface

   ! Asin function
   !> # asin
   !> Overloads the Fortran `asin` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = asin(x)
   !> ```
   interface asin
      module procedure asinHyperDual
   end interface

   ! Atan function
   !> # atan
   !> Overloads the Fortran `atan` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> y = atan(x)
   !> ```
   interface atan
      module procedure atanHyperDual
   end interface

   ! Atan2 function
   !> # atan2
   !> Overloads the Fortran `atan2` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = atan2(y, x)
   !> ```
   interface atan2
      module procedure atan2HyperDual
   end interface

   ! Max function (limited to combinations below, but that
   ! can be extended)
   !> # max
   !> Overloads the Fortran `max` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = max(x, y)
   !> ```
   interface max
      module procedure max_dd
      module procedure max_ddd
      module procedure max_dr
      module procedure max_rd
   end interface

   ! Min function (limited for now to 2 arguments, but that
   ! can be extended)
   !> # min
   !> Overloads the Fortran `min` intrinsic for hyperdual arguments, propagating supported derivative coefficients.
   !>
   !> # Description
   !> The generic interface selects a type-specific implementation based on the argument types.
   !>
   !> # Examples
   !> ```fortran
   !> z = min(x, y)
   !> ```
   interface min
      module procedure min_dd
      module procedure min_dr
      module procedure min_rd
   end interface

   !=====================================================================

contains


   !-------------------------------------------------------------------
   !--- Functions for the equal assignment. ---------------------------
   !-------------------------------------------------------------------

   !> # EqualHyperDualHyperDual
   !> Assign a real or hyperdual value.
   !>
   !> # Description
   !> Copies the source scalar or all four hyperdual coefficients into the destination.
   !>
   !> # Examples
   !> ```fortran
   !> result = value
   !> ```
   elemental subroutine EqualHyperDualHyperDual(res, inp)
      implicit none
      type (hyperdual), intent(out) :: res !! Output operation result.
      type (hyperdual), intent(in)  :: inp !! First input operand.

      res%f0  = inp%f0
      res%f1  = inp%f1
      res%f2  = inp%f2
      res%f12 = inp%f12
   end subroutine EqualHyperDualHyperDual

   !> # EqualHyperDualReal
   !> Assign a real or hyperdual value.
   !>
   !> # Description
   !> Copies the source scalar or all four hyperdual coefficients into the destination.
   !>
   !> # Examples
   !> ```fortran
   !> result = value
   !> ```
   elemental subroutine EqualHyperDualReal(res, inp)
      implicit none
      type (hyperdual), intent(out) :: res !! Output operation result.
      real(pr),         intent(in)  :: inp !! First input operand.

      res%f0  = inp
      res%f1  = 0.0_pr
      res%f2  = 0.0_pr
      res%f12 = 0.0_pr
   end subroutine EqualHyperDualReal



   !-------------------------------------------------------------------
   !--- Function for the unary operator +. ----------------------------
   !-------------------------------------------------------------------

   !> # PlusHyperDualHyperDual
   !> Apply unary plus to a hyperdual value.
   !>
   !> # Description
   !> Returns the input with derivative coefficients unchanged.
   !>
   !> # Examples
   !> ```fortran
   !> y = +x
   !> ```
   elemental function PlusHyperDualHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2

      v2%f0  = v1%f0
      v2%f1  = v1%f1
      v2%f2  = v1%f2
      v2%f12 = v1%f12
   end function PlusHyperDualHyperDual



   !-------------------------------------------------------------------
   !--- Functions for the addition operator. --------------------------
   !-------------------------------------------------------------------

   !> # AddHyperDualHyperDual
   !> Add operands.
   !>
   !> # Description
   !> Adds the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x + y
   !> ```
   elemental function AddHyperDualHyperDual(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0  + v2%f0
      v3%f1  = v1%f1  + v2%f1
      v3%f2  = v1%f2  + v2%f2
      v3%f12 = v1%f12 + v2%f12
   end function AddHyperDualHyperDual

   !> # AddHyperDualReal
   !> Add operands.
   !>
   !> # Description
   !> Adds the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x + y
   !> ```
   elemental function AddHyperDualReal(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0 + v2
      v3%f1  = v1%f1
      v3%f2  = v1%f2
      v3%f12 = v1%f12
   end function AddHyperDualReal

   !> # AddRealHyperDual
   !> Add operands.
   !>
   !> # Description
   !> Adds the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x + y
   !> ```
   elemental function AddRealHyperDual(v1,v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1 + v2%f0
      v3%f1  =      v2%f1
      v3%f2  =      v2%f2
      v3%f12 =      v2%f12
   end function AddRealHyperDual



   !-------------------------------------------------------------------
   !--- Function for the unary operator -. ----------------------------
   !-------------------------------------------------------------------

   !> # MinusHyperDualHyperDual
   !> Apply unary minus to a hyperdual value.
   !>
   !> # Description
   !> Negates the primal, tangent, and mixed-derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = -x
   !> ```
   elemental function MinusHyperDualHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2

      v2%f0  = -v1%f0
      v2%f1  = -v1%f1
      v2%f2  = -v1%f2
      v2%f12 = -v1%f12
   end function MinusHyperDualHyperDual



   !-------------------------------------------------------------------
   !--- Functions for the subtraction operator. -----------------------
   !-------------------------------------------------------------------

   !> # SubtractHyperDualHyperDual
   !> Subtract operands.
   !>
   !> # Description
   !> Subtracts the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x - y
   !> ```
   elemental function SubtractHyperDualHyperDual(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0  - v2%f0
      v3%f1  = v1%f1  - v2%f1
      v3%f2  = v1%f2  - v2%f2
      v3%f12 = v1%f12 - v2%f12
   end function SubtractHyperDualHyperDual

   !> # SubtractHyperDualReal
   !> Subtract operands.
   !>
   !> # Description
   !> Subtracts the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x - y
   !> ```
   elemental function SubtractHyperDualReal(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0 - v2
      v3%f1  = v1%f1
      v3%f2  = v1%f2
      v3%f12 = v1%f12
   end function SubtractHyperDualReal

   !> # SubtractRealHyperDual
   !> Subtract operands.
   !>
   !> # Description
   !> Subtracts the primal and corresponding derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x - y
   !> ```
   elemental function SubtractRealHyperDual(v1,v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1 - v2%f0
      v3%f1  =    - v2%f1
      v3%f2  =    - v2%f2
      v3%f12 =    - v2%f12
   end function SubtractRealHyperDual



   !-------------------------------------------------------------------
   !--- Functions for the multiplication operator. --------------------
   !-------------------------------------------------------------------

   !> # MultiplyHyperDualHyperDual
   !> Multiply operands.
   !>
   !> # Description
   !> Applies the product rule to propagate first derivatives and the mixed second derivative.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   elemental function MultiplyHyperDualHyperDual(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0 * v2%f0
      v3%f1  = v1%f0 * v2%f1 + v1%f1 * v2%f0
      v3%f2  = v1%f0 * v2%f2 + v1%f2 * v2%f0
      v3%f12 = v1%f0*v2%f12 + v1%f1*v2%f2 + v1%f2*v2%f1 + v1%f12*v2%f0
   end function MultiplyHyperDualHyperDual

   !> # MultiplyHyperDualReal
   !> Multiply operands.
   !>
   !> # Description
   !> Applies the product rule to propagate first derivatives and the mixed second derivative.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   elemental function MultiplyHyperDualReal(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0  * v2
      v3%f1  = v1%f1  * v2
      v3%f2  = v1%f2  * v2
      v3%f12 = v1%f12 * v2
   end function MultiplyHyperDualReal

   !> # MultiplyRealHyperDual
   !> Multiply operands.
   !>
   !> # Description
   !> Applies the product rule to propagate first derivatives and the mixed second derivative.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   elemental function MultiplyRealHyperDual(v1,v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1 * v2%f0
      v3%f1  = v1 * v2%f1
      v3%f2  = v1 * v2%f2
      v3%f12 = v1 * v2%f12
   end function MultiplyRealHyperDual

   !> # MultiplyHyperDualInt
   !> Multiply operands.
   !>
   !> # Description
   !> Applies the product rule to propagate first derivatives and the mixed second derivative.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   elemental function MultiplyHyperDualInt(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      integer,          intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1%f0  * v2
      v3%f1  = v1%f1  * v2
      v3%f2  = v1%f2  * v2
      v3%f12 = v1%f12 * v2
   end function MultiplyHyperDualInt

   !> # MultiplyIntHyperDual
   !> Multiply operands.
   !>
   !> # Description
   !> Applies the product rule to propagate first derivatives and the mixed second derivative.
   !>
   !> # Examples
   !> ```fortran
   !> z = x * y
   !> ```
   elemental function MultiplyIntHyperDual(v1,v2) result (v3)
      integer,          intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3%f0  = v1 * v2%f0
      v3%f1  = v1 * v2%f1
      v3%f2  = v1 * v2%f2
      v3%f12 = v1 * v2%f12
   end function MultiplyIntHyperDual



   !-------------------------------------------------------------------
   !--- Functions for the division operator. --------------------------
   !-------------------------------------------------------------------

   !> # DivideHyperDualHyperDual
   !> Divide operands.
   !>
   !> # Description
   !> Applies the quotient rule to propagate primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x / y
   !> ```
   elemental function DivideHyperDualHyperDual(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      v3 = v1 * v2**(-1)
   end function DivideHyperDualHyperDual

   !> # DivideHyperDualReal
   !> Divide operands.
   !>
   !> # Description
   !> Applies the quotient rule to propagate primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x / y
   !> ```
   elemental function DivideHyperDualReal(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3
      real(pr)                     :: invV2

      invV2 = 1.0_pr / v2
      v3 = v1 * invV2
   end function DivideHyperDualReal

   !> # DivideRealHyperDual
   !> Divide operands.
   !>
   !> # Description
   !> Applies the quotient rule to propagate primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x / y
   !> ```
   elemental function DivideRealHyperDual(v1,v2) result (v3)
      real(pr), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual) :: invV2
      type (hyperdual) :: v3

      invV2 = 1.0_pr * v2**(-1.0_pr)
      v3 = v1 * invV2
   end function DivideRealHyperDual



   !-------------------------------------------------------------------
   !--- Functions for the power operator. -----------------------------
   !-------------------------------------------------------------------

   !> # PowerHyperDualInt
   !> Raise a hyperdual value to a power.
   !>
   !> # Description
   !> Applies the power rule to the primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x ** 2
   !> ```
   elemental function PowerHyperDualInt(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      integer,          intent(in) :: v2 !! Second input operand.
      integer :: i !! Loop/index variable.
      integer :: vv2
      type (hyperdual)             :: v3

      v3  = 1.0_pr
      vv2 = abs(v2)
      do i=1,vv2
         v3 = v3*v1
      enddo

      if(v2 < 0) v3 = 1.0_pr / v3
   end function PowerHyperDualInt

   !> # PowerHyperDualHyperDual
   !> Raise a hyperdual value to a power.
   !>
   !> # Description
   !> Applies the power rule to the primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x ** y
   !> ```
   elemental function PowerHyperDualHyperDual(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual) :: v3
      type (hyperdual) :: v4

      v4 = logHyperDual(v1)
      v3 = expHyperDual(v2*v4)
   end function PowerHyperDualHyperDual

   !> # PowerHyperDualReal
   !> Raise a hyperdual value to a power.
   !>
   !> # Description
   !> Applies the power rule to the primal and derivative coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> z = x ** 2.5_pr
   !> ```
   elemental function PowerHyperDualReal(v1,v2) result (v3)
      type (hyperdual),      intent(in) :: v1 !! First input operand.
      real(pr),              intent(in) :: v2 !! Second input operand.
      type (hyperdual)                  :: v3

      real(pr), parameter               :: tol = 1.0e-15_pr !! Numerical convergence tolerance.
      real(pr) :: xval
      real(pr) :: deriv

      xval = v1%f0
      if(abs(xval) < tol) then
         if(xval >= 0.0_pr) then
            xval = tol
         else
            xval = -tol
         endif
      endif

      deriv  = v2*(xval**(v2-1.0_pr))
      v3%f0  = (v1%f0)**v2
      v3%f1  =  v1%f1 * deriv
      v3%f2  =  v1%f2 * deriv
      v3%f12 = v1%f12 * deriv &
      & + v2*(v2 - 1.0_pr)*v1%f1*v1%f2*xval**(v2-2.0_pr)
   end function PowerHyperDualReal



   !-------------------------------------------------------------------
   !--- Sum -----------------------------------------------------------
   !-------------------------------------------------------------------
   !> # SumHyperDual
   !> Sum hyperdual array elements.
   !>
   !> # Description
   !> Reduces the hyperdual coefficients using the requested mask or dimension.
   !>
   !> # Examples
   !> ```fortran
   !> total = sum(values)
   !> ```
   pure type(hyperdual) function SumHyperDual(v1, mask)
      type(hyperdual), intent(in) :: v1(:) !! First input operand.
      logical, intent(in), optional :: mask(:) !! Optional logical mask for the reduction.
      integer :: i !! Loop/index variable.

      SumHyperDual = hyperdual(0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)
      if (present(mask)) then
         do i = 1, size(v1)
            if(mask(i)) SumHyperDual = SumHyperDual + v1(i)
         end do
      else
         do i = 1, size(v1)
            SumHyperDual = SumHyperDual + v1(i)
         end do
      end if
   end function SumHyperDual

   !> # SumHyperDual2
   !> Sum hyperdual array elements.
   !>
   !> # Description
   !> Reduces the hyperdual coefficients using the requested mask or dimension.
   !>
   !> # Examples
   !> ```fortran
   !> total = sum(values)
   !> ```
   pure function SumHyperDual2(v1, dim)
      type(hyperdual), intent(in) :: v1(:,:) !! First input operand.
      integer, intent(in) :: dim !! Array dimension along which to reduce.
      type(hyperdual), allocatable :: SumHyperDual2(:)
      integer                     :: i !! Loop/index variable.

      allocate(SumHyperDual2(size(v1)/size(v1,dim)))

      SumHyperDual2 = hyperdual(0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr)
      do i = 1, size(v1,dim)
         if (dim == 1) then
            SumHyperDual2 = SumHyperDual2 + v1(i,:)
         else
            SumHyperDual2 = SumHyperDual2 + v1(:,i)
         end if
      end do
   end function SumHyperDual2



   !-------------------------------------------------------------------
   !--- Functions for the equal operator. -----------------------------
   !-------------------------------------------------------------------

   !> # eq_dd
   !> Compare operands using equality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function eq_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      eq_dd = lhs%f0 == rhs%f0
   end function eq_dd

   !> # eq_dr
   !> Compare operands using equality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   elemental logical function eq_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      eq_dr = lhs%f0 == rhs
   end function eq_dr

   !> # eq_rd
   !> Compare operands using equality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   elemental logical function eq_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      eq_rd = lhs == rhs%f0
   end function eq_rd

   !> # eq_di
   !> Compare operands using equality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function eq_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      eq_di = lhs%f0 == rhs
   end function eq_di

   !> # eq_id
   !> Compare operands using equality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function eq_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      eq_id = lhs == rhs%f0
   end function eq_id



   !-------------------------------------------------------------------
   !--- Functions for the not equal operator. -------------------------
   !-------------------------------------------------------------------

   !> # ne_dd
   !> Compare operands using inequality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ne_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      ne_dd = lhs%f0 /= rhs%f0
   end function ne_dd

   !> # ne_dr
   !> Compare operands using inequality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ne_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      ne_dr = lhs%f0 /= rhs
   end function ne_dr

   !> # ne_rd
   !> Compare operands using inequality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ne_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      ne_rd = lhs /= rhs%f0
   end function ne_rd

   !> # ne_di
   !> Compare operands using inequality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ne_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      ne_di = lhs%f0 /= rhs
   end function ne_di

   !> # ne_id
   !> Compare operands using inequality.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ne_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      ne_id = lhs /= rhs%f0
   end function ne_id



   !-------------------------------------------------------------------
   !--- Functions for the less than operator. -------------------------
   !-------------------------------------------------------------------

   !> # lt_dd
   !> Compare operands using less-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function lt_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      lt_dd = lhs%f0 < rhs%f0
   end function lt_dd

   !> # lt_dr
   !> Compare operands using less-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function lt_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      lt_dr = lhs%f0 < rhs
   end function lt_dr

   !> # lt_rd
   !> Compare operands using less-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function lt_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      lt_rd = lhs < rhs%f0
   end function lt_rd

   !> # lt_di
   !> Compare operands using less-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function lt_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      lt_di = lhs%f0 < rhs
   end function lt_di

   !> # lt_id
   !> Compare operands using less-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function lt_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      lt_id = lhs < rhs%f0
   end function lt_id



   !-------------------------------------------------------------------
   !--- Functions for the less than or equal operator. ----------------
   !-------------------------------------------------------------------

   !> # le_dd
   !> Compare operands using less-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function le_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      le_dd = lhs%f0 <= rhs%f0
   end function le_dd

   !> # le_dr
   !> Compare operands using less-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function le_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      le_dr = lhs%f0 <= rhs
   end function le_dr

   !> # le_rd
   !> Compare operands using less-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function le_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      le_rd = lhs <= rhs%f0
   end function le_rd

   !> # le_di
   !> Compare operands using less-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function le_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      le_di = lhs%f0 <= rhs
   end function le_di

   !> # le_id
   !> Compare operands using less-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function le_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      le_id = lhs <= rhs%f0
   end function le_id



   !-------------------------------------------------------------------
   !--- Functions for the greater than operator. ----------------------
   !-------------------------------------------------------------------

   !> # gt_dd
   !> Compare operands using greater-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function gt_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      gt_dd = lhs%f0 > rhs%f0
   end function gt_dd

   !> # gt_dr
   !> Compare operands using greater-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function gt_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      gt_dr = lhs%f0 > rhs
   end function gt_dr

   !> # gt_rd
   !> Compare operands using greater-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function gt_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      gt_rd = lhs > rhs%f0
   end function gt_rd

   !> # gt_di
   !> Compare operands using greater-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function gt_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      gt_di = lhs%f0 > rhs
   end function gt_di

   !> # gt_id
   !> Compare operands using greater-than.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function gt_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      gt_id = lhs > rhs%f0
   end function gt_id



   !-------------------------------------------------------------------
   !--- Functions for the greater than or equal operator. -------------
   !-------------------------------------------------------------------

   !> # ge_dd
   !> Compare operands using greater-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ge_dd(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! Left-hand operand.
      type (hyperdual), intent(in) :: rhs !! Right-hand operand.

      ge_dd = lhs%f0 >= rhs%f0
   end function ge_dd

   !> # ge_dr
   !> Compare operands using greater-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ge_dr(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      real(pr),         intent(in) :: rhs !! Second input operand.

      ge_dr = lhs%f0 >= rhs
   end function ge_dr

   !> # ge_rd
   !> Compare operands using greater-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ge_rd(lhs, rhs)
      real(pr),         intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      ge_rd = lhs >= rhs%f0
   end function ge_rd

   !> # ge_di
   !> Compare operands using greater-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ge_di(lhs, rhs)
      type (hyperdual), intent(in) :: lhs !! First input operand.
      integer,          intent(in) :: rhs !! Second input operand.

      ge_di = lhs%f0 >= rhs
   end function ge_di

   !> # ge_id
   !> Compare operands using greater-than-or-equal.
   !>
   !> # Description
   !> Returns a logical result based on primal values; derivative slots do not affect comparisons.
   !>
   !> # Examples
   !> ```fortran
   !> logical_value = x < y
   !> ```
   logical function ge_id(lhs, rhs)
      integer,          intent(in) :: lhs !! First input operand.
      type (hyperdual), intent(in) :: rhs !! Second input operand.

      ge_id = lhs >= rhs%f0
   end function ge_id



   !-------------------------------------------------------------------
   !--- Math functions. -----------------------------------------------
   !-------------------------------------------------------------------

   ! Absolute value function.
   !> # absHyperDual
   !> Evaluate the absolute value of a hyperdual value.
   !>
   !> # Description
   !> Uses the primal sign to propagate derivative coefficients through abs.
   !>
   !> # Examples
   !> ```fortran
   !> y = abs(x)
   !> ```
   elemental function absHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2

      if(v1%f0 >= 0.0) then
         v2%f0  = v1%f0
         v2%f1  = v1%f1
         v2%f2  = v1%f2
         v2%f12 = v1%f12
      else
         v2%f0  = -v1%f0
         v2%f1  = -v1%f1
         v2%f2  = -v1%f2
         v2%f12 = -v1%f12
      endif
   end function absHyperDual

   ! Integer function.
   !> # intHyperDual
   !> Convert a hyperdual value to an integer-valued result.
   !>
   !> # Description
   !> Applies integer conversion to the primal value; this discontinuous operation is not intended for derivative propagation.
   !>
   !> # Examples
   !> ```fortran
   !> y = int(x)
   !> ```
   elemental function intHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      integer                      :: v2

      v2 = int(v1%f0)
   end function intHyperDual

   ! Nearest integer function.
   !> # nintHyperDual
   !> Round a hyperdual value to the nearest integer.
   !>
   !> # Description
   !> Applies nearest-integer conversion to the primal value; this discontinuous operation is not intended for derivative propagation.
   !>
   !> # Examples
   !> ```fortran
   !> y = nint(x)
   !> ```
   elemental function nintHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      integer                      :: v2

      v2 = nint(v1%f0)
   end function nintHyperDual

   ! Real function.
   !> # realHyperDual
   !> Convert a hyperdual value to a real-valued result.
   !>
   !> # Description
   !> Returns the primal real value; derivative information is not retained.
   !>
   !> # Examples
   !> ```fortran
   !> y = real(x)
   !> ```
   elemental function realHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr)                     :: v2

      v2 = v1%f0
   end function realHyperDual

   ! Functions for the sign function.
   !> # sign_dd
   !> Apply the Fortran SIGN intrinsic to hyperdual operands.
   !>
   !> # Description
   !> Selects the sign from primal values and propagates the selected magnitude.
   !>
   !> # Examples
   !> ```fortran
   !> y = sign(x, sign_source)
   !> ```
   elemental function sign_dd(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3
      real(pr)                     :: ssign

      if(v2%f0 < 0.0) then
         ssign = -1.0
      else
         ssign =  1.0
      endif
      v3 = ssign*v1
   end function sign_dd

   !> # sign_dr
   !> Apply the Fortran SIGN intrinsic to hyperdual operands.
   !>
   !> # Description
   !> Selects the sign from primal values and propagates the selected magnitude.
   !>
   !> # Examples
   !> ```fortran
   !> y = sign(x, sign_source)
   !> ```
   elemental function sign_dr(v1,v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3
      real(pr)                     :: ssign

      if(v2 < 0.0) then
         ssign = -1.0
      else
         ssign =  1.0
      endif
      v3 = ssign*v1
   end function sign_dr

   !> # sign_rd
   !> Apply the Fortran SIGN intrinsic to hyperdual operands.
   !>
   !> # Description
   !> Selects the sign from primal values and propagates the selected magnitude.
   !>
   !> # Examples
   !> ```fortran
   !> y = sign(x, sign_source)
   !> ```
   elemental function sign_rd(v1,v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3
      real(pr)                     :: ssign

      if(v2%f0 < 0.0) then
         ssign = -1.0
      else
         ssign =  1.0
      endif
      v3 = ssign*v1
   end function sign_rd

   ! Sine function.
   !> # sinHyperDual
   !> Evaluate the sin intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = sin(x)
   !> ```
   elemental function sinHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: f
      real(pr) :: dx

      f = sin(v1%f0)
      dx = cos(v1%f0)
      v2%f0  = f
      v2%f1  = dx * v1%f1
      v2%f2  = dx * v1%f2
      v2%f12 = dx * v1%f12 - f * v1%f1 * v1%f2
   end function sinHyperDual

   ! Cosine function.
   !> # cosHyperDual
   !> Evaluate the cos intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = cos(x)
   !> ```
   elemental function cosHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: f
      real(pr) :: dx

      f = cos(v1%f0)
      dx = -sin(v1%f0)
      v2%f0  = f
      v2%f1  = dx * v1%f1
      v2%f2  = dx * v1%f2
      v2%f12 = dx * v1%f12 - f * v1%f1 * v1%f2
   end function cosHyperDual

   ! Tangent function.
   !> # tanHyperDual
   !> Evaluate the tan intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = tan(x)
   !> ```
   elemental function tanHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: f
      real(pr) :: dx

      f = tan(v1%f0)
      dx = f * f + 1.0_pr
      v2%f0  = f
      v2%f1  = dx * v1%f1
      v2%f2  = dx * v1%f2
      v2%f12 = dx * v1%f12 + v1%f1 * v1%f2 * 2.0_pr * f * dx
   end function tanHyperDual

   ! Sqrt function
   !> # sqrtHyperDual
   !> Evaluate the sqrt intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = sqrt(x)
   !> ```
   elemental function sqrtHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr), parameter :: expo=3.0_pr/2.0_pr
      real(pr) :: square

      square = sqrt(v1%f0)

      v2%f0 = square
      v2%f1 = 0.5_pr / square * v1%f1
      v2%f2 = 0.5_pr / square * v1%f2

      v2%f12 = 0.5_pr * v1%f12 / square - 0.25_pr * v1%f1 * v1%f2 / (v1%f0**expo)
   end function sqrtHyperDual

   ! Log function
   !> # logHyperDual
   !> Evaluate the log intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = log(x)
   !> ```
   elemental function logHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: dx1
      real(pr) :: dx2

      dx1 = v1%f1 / v1%f0
      dx2 = v1%f2 / v1%f0
      v2%f0  = log(v1%f0)
      v2%f1  = dx1
      v2%f2  = dx2
      v2%f12 = v1%f12 / v1%f0 - (dx1 * dx2)
   end function logHyperDual

   ! Log10 function
   !> # log10HyperDual
   !> Evaluate the log10 intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = log10(x)
   !> ```
   elemental function log10HyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2

      v2 = log(v1)/log(10.0_pr)
   end function log10HyperDual

   ! Exp function
   !> # expHyperDual
   !> Evaluate the exp intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = exp(x)
   !> ```
   elemental function expHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr)                     :: dx

      dx = exp(v1%f0)
      v2%f0  = dx
      v2%f1  = dx * v1%f1
      v2%f2  = dx * v1%f2
      v2%f12 = dx * (v1%f12 + v1%f1 * v1%f2)
   end function expHyperDual

   ! Sinh function
   !> # sinhHyperDual
   !> Evaluate the sinh intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = sinh(x)
   !> ```
   elemental function sinhHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual) :: t1
      type (hyperdual) :: t2
      type (hyperdual) :: v2

      t1 = exp(v1)
      t2 = exp(-v1)
      v2 = 0.5_pr*(t1-t2)
   end function sinhHyperDual

   ! Cosh function
   !> # coshHyperDual
   !> Evaluate the cosh intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = cosh(x)
   !> ```
   elemental function coshHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual) :: t1
      type (hyperdual) :: t2
      type (hyperdual) :: v2

      t1 = exp(v1)
      t2 = exp(-v1)
      v2 = 0.5_pr*(t1+t2)
   end function coshHyperDual

   ! Tanh function
   !> # tanhHyperDual
   !> Evaluate the tanh intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = tanh(x)
   !> ```
   elemental function tanhHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual) :: t1
      type (hyperdual) :: t2
      type (hyperdual) :: v2

      t1 = exp(v1)
      t2 = exp(-v1)
      v2 = (t1-t2)/(t1+t2)
   end function tanhHyperDual

   ! Acos function
   !> # acosHyperDual
   !> Evaluate the acos intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = acos(x)
   !> ```
   elemental function acosHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: deriv
      real(pr) :: deriv1

      deriv1 = 1.0_pr - v1%f0*v1%f0
      deriv  = -1.0_pr / sqrt(deriv1)
      v2%f0  = acos(v1%f0)
      v2%f1  = deriv*v1%f1
      v2%f2  = deriv*v1%f2
      v2%f12 = deriv*v1%f12 &
      & + v1%f1 * v1%f2 * (-v1%f0 * deriv1**(-1.5_pr))
   end function acosHyperDual

   ! Asin function
   !> # asinHyperDual
   !> Evaluate the asin intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = asin(x)
   !> ```
   elemental function asinHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: deriv
      real(pr) :: deriv1

      deriv1 = 1.0_pr - v1%f0*v1%f0
      deriv  = 1.0_pr / sqrt(deriv1)
      v2%f0  = asin(v1%f0)
      v2%f1  = deriv*v1%f1
      v2%f2  = deriv*v1%f2
      v2%f12 = deriv*v1%f12 &
      & + v1%f1 * v1%f2 * (v1%f0 * deriv1**(-1.5_pr))
   end function asinHyperDual

   ! Atan function
   !> # atanHyperDual
   !> Evaluate the atan intrinsic for a hyperdual value.
   !>
   !> # Description
   !> Applies the intrinsic and its derivatives to the hyperdual coefficients.
   !>
   !> # Examples
   !> ```fortran
   !> y = atan(x)
   !> ```
   elemental function atanHyperDual(v1) result (v2)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual)             :: v2
      real(pr) :: deriv
      real(pr) :: deriv1

      deriv1 = 1.0_pr + v1%f0*v1%f0
      deriv  = 1.0_pr / deriv1
      v2%f0  = atan(v1%f0)
      v2%f1  = deriv*v1%f1
      v2%f2  = deriv*v1%f2
      v2%f12 = deriv*v1%f12 &
      & + v1%f1 * v1%f2 * (-2.0_pr * v1%f0 / (deriv1 * deriv1))
   end function atanHyperDual

   ! Atan2 function
   !> # atan2HyperDual
   !> Evaluate `atan2(v1, v2)` for hyperdual arguments.
   !>
   !> # Description
   !> Only the primal value and the \(\varepsilon_1\) coefficient `f1` are
   !> propagated; `f2` and `f12` of the result are left at zero, so this
   !> function must not be used where second derivatives are needed.
   !>
   !> # Examples
   !> ```fortran
   !> z = atan2(y, x)
   !> ```
   elemental function atan2HyperDual(v1, v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3
      real(pr) :: a
      real(pr) :: b
      real(pr) :: c
      real(pr) :: d

      a = v1%f0
      b = v1%f1
      c = v2%f0
      d = v2%f1

      v3%f0 = atan2(a,c)
      v3%f1 = (c*b - a*d)/(a*a + c*c)
   end function atan2HyperDual

   ! Max functions
   !> # max_dd
   !> Select the maximum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = max(x, y)
   !> ```
   elemental function max_dd(v1, v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1%f0 > v2%f0) then
         v3 = v1
      else
         v3 = v2
      endif
   end function max_dd

   !> # max_ddd
   !> Select the maximum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> w = max(x, y, z)
   !> ```
   elemental function max_ddd(v1, v2, v3) result (v4)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual), intent(in) :: v3 !! Result of the operation.
      type (hyperdual)             :: v4

      if(v1%f0 > v2%f0) then
         v4 = v1
      else
         v4 = v2
      endif

      if(v3%f0 > v4%f0) v4 = v3
   end function max_ddd

   !> # max_dr
   !> Select the maximum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = max(x, y)
   !> ```
   elemental function max_dr(v1, v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1%f0 > v2) then
         v3 = v1
      else
         v3 = v2
      endif
   end function max_dr

   !> # max_rd
   !> Select the maximum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = max(x, y)
   !> ```
   elemental function max_rd(v1, v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1 > v2%f0) then
         v3 = v1
      else
         v3 = v2
      endif
   end function max_rd

   ! Min functions
   !> # min_dd
   !> Select the minimum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = min(x, y)
   !> ```
   elemental function min_dd(v1, v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1%f0 < v2%f0) then
         v3 = v1
      else
         v3 = v2
      endif
   end function min_dd

   !> # min_dr
   !> Select the minimum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = min(x, y)
   !> ```
   elemental function min_dr(v1, v2) result (v3)
      type (hyperdual), intent(in) :: v1 !! First input operand.
      real(pr),         intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1%f0 < v2) then
         v3 = v1
      else
         v3 = v2
      endif
   end function min_dr

   !> # min_rd
   !> Select the minimum hyperdual operand.
   !>
   !> # Description
   !> Selects the branch using primal values and returns the selected hyperdual operand.
   !>
   !> # Examples
   !> ```fortran
   !> z = min(x, y)
   !> ```
   elemental function min_rd(v1, v2) result (v3)
      real(pr),         intent(in) :: v1 !! First input operand.
      type (hyperdual), intent(in) :: v2 !! Second input operand.
      type (hyperdual)             :: v3

      if(v1 < v2%f0) then
         v3 = v1
      else
         v3 = v2
      endif
   end function min_rd
end module
