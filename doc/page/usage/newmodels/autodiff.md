---
title: Automatic differentiation
---
# Autodiff
The implementation of new models and all their required derivatives can be
an endeavours and error-prone task. A tool that helps with this, at a small
performance cost, can be automatic differentiation. 

Automatic differentiation can be implemented in two ways:

- Forward Differentiation
- Backward Differentiation

With can be combined to obtain higher order derivatives.

In `yaeos` it is possible to add new models via two different kind of
implementations. Operator overloading with `hyperdual` numbers and
source-to-source automatic differentiation with `tapenade`.

@warn
Remember to use the `R`constant from `yaeos__constants`, and all models
should have a `type(Substances)` attribute!
@endwarn

## Hyperdual autodiff

### ArModel
Automatic differentiation with `hyperdual` numbers can be done with the
[[ArModelAdiff]] derived type. This implementation requires just to extend
that derived type with your own implementation and a volume initializer.

```fortran
module hyperdual_pr76
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   use yaeos__constants, only: pr, R
   use yaeos__substance, only: Substances
   implicit none

   type, extends(ArModelAdiff) :: PR76
      !! PengRobinson 76 EoS
      
      ! Mixing rule Parameters
      real(pr), allocatable :: kij(:, :), lij(:, :)

      ! EoS parameters
      real(pr), allocatable :: ac(:), b(:), k(:)
      real(pr), allocatable :: tc(:), pc(:), w(:)
   contains
      procedure :: Ar => arfun
      procedure :: get_v0 => v0
      procedure :: volume => volume
   end type PR76

   real(pr), parameter :: del1 = 1._pr + sqrt(2._pr)
   real(pr), parameter :: del2 = 1._pr - sqrt(2._pr)

contains

    type(PR76) function setup(tc, pc, w, kij, lij) result(self)
        !! Function to obtain a defined PR76 model with setted up parameters
        !! as function of Tc, Pc, and w
        real(pr) :: tc(:)
        real(pr) :: pc(:)
        real(pr) :: w(:)
        real(pr) :: kij(:, :)
        real(pr) :: lij(:, :)

        self%composition%tc = tc
        self%composition%pc = pc
        self%composition%w = w

        self%ac = 0.45723553_pr * R**2 * tc**2 / pc
        self%b = 0.07779607_pr * R * tc_in/pc_in
        self%k = 0.37464_pr + 1.54226_pr * w - 0.26993_pr * w**2

        self%kij = kij
        self%lij = lij
    end function

    function arfun(self, n, v, t) result(ar)
        !! Residual Helmholtz calculation for a generic cubic with
        !! quadratic mixing rules.
        class(PR76) :: self
        type(hyperdual), intent(in) :: n(:), v, t
        type(hyperdual) :: ar
    
        type(hyperdual) :: amix, a(size(n)), ai(size(n)), n2(size(n))
        type(hyperdual) :: bmix
        type(hyperdual) :: b_v, nij

        integer :: i, j

        ! Associate allows us to keep the later expressions simple.
        associate(&
            pc => self%composition%pc, ac => self%ac, b => self%b, k => self%k,&
            kij => self%kij, lij => self%lij, tc => self%compostion%tc & 
            )

            ! Soave alpha function
            a = 1.0_pr + k * (1.0_pr - sqrt(t/tc))
            a = ac * a ** 2
            ai = sqrt(a)

            ! Quadratic Mixing Rule
            amix = 0.0_pr
            bmix = 0.0_pr

            do i=1,size(n)-1
                do j=i+1,size(n)
                    nij = n(i) * n(j)
                    amix = amix + 2 * nij * (ai(i) * ai(j)) * (1 - kij(i, j))
                    bmix = bmix + nij * (b(i) + b(j)) * (1 - lij(i, j))
                end do
            end do
         end do

         amix = amix + sum(n**2*a)
         bmix = bmix + sum(n**2 * b)

         bmix = bmix/sum(n)

         b_v = bmix/v

         ! Generic Cubic Ar function
         ar = (&
            - sum(n) * log(1.0_pr - b_v) &
            - amix / (R*T*bmix)*1.0_pr / (del1 - del2) &
            * log((1.0_pr + del1 * b_v) / (1.0_pr + del2 * b_v)) &
            ) * (R * T)

      end associate
   end function arfun

   function v0(self, n, p, t)
      !! Initialization of liquid volume solving with covolume. This also
      !! helps the Michelsen volume solver
      class(PR76), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: p
      real(pr), intent(in) :: t
      real(pr) :: v0

      v0 = sum(n * self%b) / sum(n)
   end function v0

   subroutine volume(eos, n, P, T, V, root_type)
      !! In the case of models that have a "covolume" value, using the solver
      !! of Michelsen is a better option that the default.
      use yaeos__models_solvers, only: volume_michelsen
      class(PR76), intent(in) :: eos
      real(pr), intent(in)  :: n(:), P, T
      real(pr), intent(out) :: V
      character(len=*), intent(in) :: root_type

      call volume_michelsen(eos, n, P, T, V, root_type)
   end subroutine
end module hyperdual_pr76
```


## Tapenade Adiff
And alternative to `hyperdual` that takes a bit more work, but can end in a more
performant model, is doing `tapenade` source-to-source differentiation. For
this `tapenade` must be installed and accessible from a terminal 
[donwload link](https://tapenade.gitlabpages.inria.fr/userdoc/build/html/download.html).

{!tools/tapenade_diff/README.md!}