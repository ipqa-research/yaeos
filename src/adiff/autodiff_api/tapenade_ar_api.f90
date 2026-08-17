module yaeos__tapenade_ar_api
   !! Module that wraps tapenade generated routines to calculate !
   !! Ar and derivatives.
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   implicit none

   private

   public :: ArModelTapenade

   type, abstract, extends(ArModel) :: ArModelTapenade
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
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(out) :: arval
      end subroutine tapenade_ar

      subroutine tapenade_ar_d(model, n, nd, v, vd, t, td, arval, arvald)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald
      end subroutine tapenade_ar_d

      subroutine tapenade_ar_b(model, n, nb, v, vb, t, tb, arval, arvalb)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arvalb
         real(pr) :: nb(:), vb, tb
         real(pr) :: arval
      end subroutine tapenade_ar_b

      subroutine tapenade_ar_d_b(model, &
         n, nb, nd, ndb, v, vb, vd, vdb, t, tb, td, tdb, &
         arval, arvalb, arvald, arvaldb)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arval

         real(pr), intent(in) :: nd(:), vd, td
         real(pr) :: arvald

         real(pr) :: nb(:), vb, tb
         real(pr) :: arvalb

         real(pr) :: ndb(:), vdb, tdb
         real(pr) :: arvaldb
      end subroutine tapenade_ar_d_b

      subroutine tapenade_ar_d_d(model, n, nd, v, vd0, vd, t, td0, td, &
         arval, arvald0, arvald, arvaldd)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: vd0, td0
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald0, arvald, arvaldd
      end subroutine tapenade_ar_d_d
   end interface

contains

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz model generic interface
      class(ArModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: v, t
      real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
      real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))

      real(pr) :: nb(size(n)), nd(size(n)), ndb(size(n))
      real(pr) :: vb, vd, vdb, vd0
      real(pr) :: tb, td, tdb, td0
      real(pr) :: arval, arvalb, arvald, arvaldb, arvald0, arvaldd

      integer :: i, nc
      logical :: need_Arn2, need_ArVn, need_ArTn
      logical :: need_ArV2, need_ArT2, need_ArTV
      logical :: need_grad, need_val

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
            ArTV = vb  ! vb is d2Ar/(dt*dv)
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
