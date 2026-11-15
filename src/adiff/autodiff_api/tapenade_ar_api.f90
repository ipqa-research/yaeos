module yaeos__tapenade_ar_api
   !! # Tapenade residual Helmholtz API
   !! Abstract interface and adapter for source-to-source differentiated EOS models.
   !!
   !! # Description
   !! `ArModelTapenade` requires primal, forward, reverse, mixed, and second
   !! directional derivative routines. Its implementation of
   !! `residual_helmholtz` dispatches to those routines to populate requested
   !! derivative outputs.
   !!
   !! # Examples
   !! See `PCSAFT` in `src/models/residual_helmholtz/saft/pc_saft.f90` for a
   !! production implementation and `tools/tapenade_diff/` for the generator.
   !!
   !! # References
   !! Tapenade documentation: https://tapenade.gitlabpages.inria.fr/
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   implicit none

   private

   public :: ArModelTapenade

   type, abstract, extends(ArModel) :: ArModelTapenade
      !! # ArModelTapenade
      !! Base class for EOS models with Tapenade-generated derivatives.
      !!
      !! # Description
      !! Concrete extensions provide the primal `ar` routine and its
      !! forward/reverse derivative variants. The inherited adapter evaluates
      !! requested `ArModel` derivatives by dispatching to those routines.
      !!
      !! # Examples
      !! See the `PCSAFT` model for an integrated implementation.
   contains
      procedure(tapenade_ar), deferred :: ar
      procedure(tapenade_ar_d), deferred :: ar_d
      procedure(tapenade_ar_b), deferred :: ar_b
      procedure(tapenade_ar_d_b), deferred :: ar_d_b
      procedure(tapenade_ar_d_d), deferred :: ar_d_d
      procedure :: residual_helmholtz => residual_helmholtz
   end type ArModelTapenade

   abstract interface
      subroutine tapenade_ar(model, n, v, t, arval)
         !! # tapenade_ar
         !! Primal residual Helmholtz callback required by `ArModelTapenade`.
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model !! Model parameters.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(out) :: arval !! Residual Helmholtz energy.
      end subroutine tapenade_ar

      subroutine tapenade_ar_d(model, n, nd, v, vd, t, td, arval, arvald)
         !! # tapenade_ar_d
         !! Forward directional derivative callback for the residual Helmholtz energy.
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model !! Model parameters.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(in) :: nd(:) !! Tangent seed for component mole amounts (nd).
         real(pr), intent(in) :: vd !! Tangent seed for volume (vd).
         real(pr), intent(in) :: td !! Tangent seed for temperature (td).
         real(pr), intent(out) :: arval !! Primal residual Helmholtz energy.
         real(pr), intent(out) :: arvald !! Forward directional derivative of residual Helmholtz energy.
      end subroutine tapenade_ar_d

      subroutine tapenade_ar_b(model, n, nb, v, vb, t, tb, arval, arvalb)
         !! # tapenade_ar_b
         !! Reverse-mode adjoint callback for the residual Helmholtz energy.
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model !! Model parameters.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr) :: arvalb !! Incoming adjoint seed on the energy.
         real(pr) :: nb(:) !! Adjoint value for component mole amounts.
         real(pr) :: vb !! Adjoint value for volume.
         real(pr) :: tb !! Adjoint value for temperature.
         real(pr) :: arval !! Primal energy value.
      end subroutine tapenade_ar_b

      subroutine tapenade_ar_d_b(model, &
         n, nb, nd, ndb, v, vb, vd, vdb, t, tb, td, tdb, &
         arval, arvalb, arvald, arvaldb)
         !! # tapenade_ar_d_b
         !! Mixed forward/reverse callback used for Hessian and mixed derivatives.
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model !! Model parameters.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr) :: arval !! Primal energy value.

         real(pr), intent(in) :: nd(:) !! Forward tangent seed for the mole amounts.
         real(pr), intent(in) :: vd !! Forward tangent seed for the volume.
         real(pr), intent(in) :: td !! Forward tangent seed for the temperature.
         real(pr) :: arvald !! Forward directional energy derivative.

         real(pr) :: nb(:) !! Adjoint value for component mole amounts.
         real(pr) :: vb !! Adjoint value for volume.
         real(pr) :: tb !! Adjoint value for temperature.
         real(pr) :: arvalb !! Adjoint seed for the primal energy output.

         real(pr) :: ndb(:) !! Mixed tangent-adjoint value for component mole amounts.
         real(pr) :: vdb !! Mixed tangent-adjoint value for volume.
         real(pr) :: tdb !! Mixed tangent-adjoint value for temperature.
         real(pr) :: arvaldb !! Mixed tangent-adjoint energy output.
      end subroutine tapenade_ar_d_b

      subroutine tapenade_ar_d_d(model, n, nd, v, vd0, vd, t, td0, td, &
         arval, arvald0, arvald, arvaldd)
         !! # tapenade_ar_d_d
         !! Nested forward-mode callback for second directional derivatives.
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model !! Model parameters.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: v !! Total volume [L].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(in) :: vd0 !! Tangent seed for volume (vd0).
         real(pr), intent(in) :: td0 !! Tangent seed for temperature (td0).
         real(pr), intent(in) :: nd(:) !! Tangent seed for component mole amounts (nd).
         real(pr), intent(in) :: vd !! Tangent seed for volume (vd).
         real(pr), intent(in) :: td !! Tangent seed for temperature (td).
         real(pr), intent(out) :: arval !! Primal residual Helmholtz energy.
         real(pr), intent(out) :: arvald0 !! First nested tangent of residual Helmholtz energy.
         real(pr), intent(out) :: arvald !! Forward directional derivative of residual Helmholtz energy.
         real(pr), intent(out) :: arvaldd !! Mixed second tangent of residual Helmholtz energy.
      end subroutine tapenade_ar_d_d
   end interface

contains

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz model generic interface
      class(ArModelTapenade), intent(in) :: self !! Model instance.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: v !! Total volume [L].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), optional, intent(out) :: Ar !! Residual Helmholtz energy [bar L].
      real(pr), optional, intent(out) :: ArV !! \(\partial A^r/\partial V\) [bar].
      real(pr), optional, intent(out) :: ArT !! \(\partial A^r/\partial T\) [bar L/K].
      real(pr), optional, intent(out) :: ArT2 !! \(\partial^2 A^r/\partial T^2\) [bar L/K^2].
      real(pr), optional, intent(out) :: ArTV !! \(\partial^2 A^r/\partial T\partial V\) [bar/K].
      real(pr), optional, intent(out) :: ArV2 !! \(\partial^2 A^r/\partial V^2\) [bar/L].
      real(pr), optional, dimension(size(n)), intent(out) :: Arn !! \(\partial A^r/\partial n_i\) [bar L/mol].
      real(pr), optional, dimension(size(n)), intent(out) :: ArVn !! \(\partial^2 A^r/\partial V\partial n_i\) [bar/mol].
      real(pr), optional, dimension(size(n)), intent(out) :: ArTn !! \(\partial^2 A^r/\partial T\partial n_i\) [bar L/(K mol)].
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\partial^2 A^r/\partial n_i\partial n_j\) [bar L/mol^2].

      real(pr) :: nb(size(n))
      real(pr) :: nd(size(n))
      real(pr) :: ndb(size(n))
      real(pr) :: vb
      real(pr) :: vd
      real(pr) :: vdb
      real(pr) :: vd0
      real(pr) :: tb
      real(pr) :: td
      real(pr) :: tdb
      real(pr) :: td0
      real(pr) :: arval
      real(pr) :: arvalb
      real(pr) :: arvald
      real(pr) :: arvaldb
      real(pr) :: arvald0
      real(pr) :: arvaldd

      integer :: i !! Loop/index variable.
      integer :: nc !! Number of components.
      logical :: need_Arn2
      logical :: need_ArVn
      logical :: need_ArTn
      logical :: need_ArV2
      logical :: need_ArT2
      logical :: need_ArTV
      logical :: need_grad
      logical :: need_val

      nc = size(n)

      ! Initialize tracking flags based on requested arguments
      need_Arn2 = present(Arn2)
      need_ArVn = present(ArVn)
      need_ArTn = present(ArTn)
      need_ArV2 = present(ArV2)
      need_ArT2 = present(ArT2)
      need_ArTV = present(ArTV)
      need_grad = present(Arn) .or. present(ArV) .or. present(ArT)
      need_val  = present(Ar)

      ! 1. Hessian matrix w.r.t composition (Uses backward pass: does NOT yield primal Ar)
      if (need_Arn2) then
         do i = 1, nc
            call reset_vars
            arvaldb = 1
            nd(i) = 1
            call self%ar_d_b( &
               n, nb, nd, ndb, &
               v, vb, vd, vdb, &
               t, tb, td, tdb, &
               arval, arvalb, arvald, arvaldb)

            Arn2(i, :) = nb

            ! Extract free cross-derivatives
            if (need_ArVn) ArVn(i) = vb
            if (need_ArTn) ArTn(i) = tb

            ! Gradients are identical across all iterations of this loop
            if (i == 1) then
               if (present(Arn)) Arn = ndb
               if (present(ArV)) ArV = vdb
               if (present(ArT)) ArT = tdb
               need_grad = .false.
            end if
         end do
         need_ArVn = .false.
         need_ArTn = .false.
      end if

      ! 2. Cross derivative w.r.t Volume and Composition (Uses backward pass)
      if (need_ArVn) then
         call reset_vars
         arvaldb = 1
         vd = 1
         call self%ar_d_b( &
            n, nb, nd, ndb, &
            v, vb, vd, vdb, &
            t, tb, td, tdb, &
            arval, arvalb, arvald, arvaldb)

         ArVn = nb

         ! Free byproducts from vd=1 sweep
         if (need_ArV2) then
            ArV2 = vb
            need_ArV2 = .false.
         end if
         if (need_ArTV) then
            ArTV = tb
            need_ArTV = .false.
         end if
         if (need_grad) then
            if (present(Arn)) Arn = ndb
            if (present(ArV)) ArV = vdb
            if (present(ArT)) ArT = tdb
            need_grad = .false.
         end if
      end if

      ! 3. Cross derivative w.r.t Temperature and Composition (Uses backward pass)
      if (need_ArTn) then
         call reset_vars
         arvaldb = 1
         td = 1
         call self%ar_d_b( &
            n, nb, nd, ndb, &
            v, vb, vd, vdb, &
            t, tb, td, tdb, &
            arval, arvalb, arvald, arvaldb)

         ArTn = nb

         ! Free byproducts from td=1 sweep
         if (need_ArT2) then
            ArT2 = tb
            need_ArT2 = .false.
         end if
         if (need_ArTV) then
            ArTV = vb ! vb is d2Ar/(dt*dv)
            need_ArTV = .false.
         end if
         if (need_grad) then
            if (present(Arn)) Arn = ndb
            if (present(ArV)) ArV = vdb
            if (present(ArT)) ArT = tdb
            need_grad = .false.
         end if
      end if

      ! 4. Pure Gradients (Uses backward pass: does NOT yield primal Ar)
      if (need_grad) then
         call reset_vars
         arvalb = 1
         call self%ar_b(n, nb, v, vb, t, tb, arval, arvalb)
         if (present(Arn)) Arn = nb
         if (present(ArV)) ArV = vb
         if (present(ArT)) ArT = tb
      end if

      ! 5. Scalar 2nd derivatives (Uses forward sweep: DOES yield primal Ar)
      if (need_ArV2) then
         call reset_vars
         vd = 1; vd0 = 1
         call self%ar_d_d(n, nd, v, vd0, vd, t, td0, td, arval, arvald0, arvald, arvaldd)
         ArV2 = arvaldd
      end if

      if (need_ArT2) then
         call reset_vars
         td = 1; td0 = 1
         call self%ar_d_d(n, nd, v, vd0, vd, t, td0, td, arval, arvald0, arvald, arvaldd)
         ArT2 = arvaldd
      end if

      if (need_ArTV) then
         call reset_vars
         vd = 1; td0 = 1
         call self%ar_d_d(n, nd, v, vd0, vd, t, td0, td, arval, arvald0, arvald, arvaldd)
         ArTV = arvaldd
      end if

      ! 6. Pure Evaluation (Triggered if Ar is needed and NO forward sweep caught it)
      if (need_val) then
         call reset_vars
         call self%ar(n, v, t, arval)
         Ar = arval
      end if

   contains
      subroutine reset_vars
         nb=0; nd=0; ndb=0
         vb=0; vd=0; vd0=0; vdb=0
         tb=0; td=0; td0=0; tdb=0
         arval=0; arvalb=0; arvald=0; arvald0=0; arvaldb=0; arvaldd=0
      end subroutine reset_vars
   end subroutine residual_helmholtz
end module yaeos__tapenade_ar_api
