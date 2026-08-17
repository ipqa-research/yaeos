module yaeos__tapenade_ar_api
   !! Module that wraps tapenade generated routines to calculate
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
   end type

   abstract interface
      subroutine tapenade_ar(model, n, v, t, arval)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(out) :: arval
      end subroutine

      subroutine tapenade_ar_d(model, n, nd, v, vd, t, td, arval, arvald)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald
      end subroutine

      subroutine tapenade_ar_b(model, n, nb, v, vb, t, tb, arval, arvalb)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arvalb
         real(pr) :: nb(:), vb, tb
         real(pr) :: arval
      end subroutine

      subroutine tapenade_ar_d_b(model, &
         n, nb, nd, ndb, v, vb, vd, vdb, t, tb, td, tdb, &
         arval, arvalb, arvald, arvaldb)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t

         real(pr), intent(in) :: nd(:), vd, td
         real(pr) :: arvald

         real(pr) :: nb(:), vb, tb
         real(pr) :: arvalb

         real(pr) :: ndb(:), vdb, tdb
         real(pr) :: arvaldb

         real(pr) :: arval
      end subroutine

      subroutine tapenade_ar_d_d(model, n, nd, v, vd0, vd, t, td0, td, &
         arval, arvald0, arvald, arvaldd)
         import pr, ArModelTapenade
         class(ArModelTapenade), intent(in) :: model
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: vd0, td0
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald0, arvald, arvaldd
      end subroutine
   end interface

contains

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Residual Helmholtz model generic interface
      !! FIXED: Handles all derivative combinations efficiently with no redundant calls
      class(ArModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: v, t
      real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
      real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))

      ! Working arrays
      real(pr) :: nb(size(n)), nd(size(n)), ndb(size(n))
      real(pr) :: vb, vd, vdb, vd0
      real(pr) :: tb, td, tdb, td0
      real(pr) :: arval, arvalb, arvald, arvaldb, arvald0, arvaldd

      ! Flags to track what's been computed
      logical :: computed_arn, computed_arv, computed_art
      integer :: i, nc

      nc = size(n)

      ! =========================================================================
      ! PATH 1: Arn2 requested - use forward-on-backward loop
      ! This efficiently computes: Arn2, Arn, ArV, ArT, ArVn, ArTn
      ! =========================================================================
      if (present(Arn2)) then
         computed_arn = .false.
         computed_arv = .false.
         computed_art = .false.

         do i = 1, nc
            call reset_vars

            ! Set up forward direction (component i)
            nd(i) = 1

            ! Seed for backward sweep (computes second derivatives)
            arvaldb = 1

            call self%ar_d_b(&
               n, nb, nd, ndb, &
               v, vb, vd, vdb, &
               t, tb, td, tdb, &
               arval, arvalb, arvald, arvaldb &
            )

            ! Extract second derivatives (Hessian)
            ! Arn2(i, :) = ndb
            Arn2(i, :) = nb

            ! Save first derivatives on first iteration (identical for all i)
            if (i == 1) then
               if (present(Arn)) then
                  Arn = nb
                  computed_arn = .true.
               end if
               if (present(ArV)) then
                  ArV = vb
                  computed_arv = .true.
               end if
               if (present(ArT)) then
                  ArT = tb
                  computed_art = .true.
               end if
            end if

            ! Save mixed second derivatives (optimized: avoid redundant calls)
            if (present(ArVn)) ArVn(i) = vdb
            if (present(ArTn)) ArTn(i) = tdb
         end do

      ! =========================================================================
      ! PATH 2: Only Arn requested (no Arn2) - use simple reverse mode
      ! =========================================================================
      else if (present(Arn)) then
         call reset_vars
         arvalb = 1
         call self%ar_b(n, nb, v, vb, t, tb, arval, arvalb)
         Arn = nb
         computed_arn = .true.

         ! Grab first derivatives as byproduct
         if (present(ArV)) then
            ArV = vb
            computed_arv = .true.
         end if
         if (present(ArT)) then
            ArT = tb
            computed_art = .true.
         end if

      ! =========================================================================
      ! PATH 3: No composition derivatives, but V,T first derivatives needed
      ! =========================================================================
      else if ((present(ArV) .or. present(ArT)) .and. &
               .not. present(ArVn) .and. .not. present(ArTn)) then
         call reset_vars
         arvalb = 1
         call self%ar_b(n, nb, v, vb, t, tb, arval, arvalb)

         if (present(ArV)) then
            ArV = vb
            computed_arv = .true.
         end if
         if (present(ArT)) then
            ArT = tb
            computed_art = .true.
         end if
      else
         computed_arn = .false.
         computed_arv = .false.
         computed_art = .false.
      end if

      ! =========================================================================
      ! Mixed second derivatives (ArVn, ArTn) if not already computed above
      ! =========================================================================
      if (present(ArVn) .and. .not. present(Arn2)) then
         ArVn = get_ArnX("V")
      end if

      if (present(ArTn) .and. .not. present(Arn2)) then
         ArTn = get_ArnX("T")
      end if

      ! =========================================================================
      ! Pure second derivatives (ArV2, ArT2, ArTV) using forward-on-forward
      ! =========================================================================
      if (present(ArTV)) ArTV = get_dArdX2("TV")
      if (present(ArT2)) ArT2 = get_dArdX2("T2")
      if (present(ArV2)) ArV2 = get_dArdX2("V2")

      ! =========================================================================
      ! The function value itself (if requested and not yet computed)
      ! =========================================================================
      if (present(Ar)) then
         ! Ar is computed in most paths; if not, compute it now
         if (.not. (present(Arn2) .or. present(Arn) .or. &
                     (present(ArV) .and. computed_arv) .or. &
                     (present(ArT) .and. computed_art))) then
            call reset_vars
            call self%ar(n, v, t, arval)
         end if
         Ar = arval
      end if

   contains

      subroutine reset_vars
         !! Reset all working variables
         nb = 0
         nd = 0
         ndb = 0

         vb = 0
         vd = 0
         vd0 = 0
         vdb = 0

         tb = 0
         td = 0
         td0 = 0
         tdb = 0

         arval = 0
         arvalb = 0
         arvald = 0
         arvald0 = 0
         arvaldb = 0
      end subroutine

      function get_dArdX2(var)
         !! Compute pure second derivatives using forward-on-forward AD
         !! Computes: ArV2, ArT2, ArTV
         character(len=*), intent(in) :: var
         real(pr) :: get_dArdX2

         call reset_vars

         select case (var)
         case ("TV")
            vd = 1
            td0 = 1
         case ("V2")
            vd = 1
            vd0 = 1
         case ("T2")
            td = 1
            td0 = 1
         end select

         call self%ar_d_d(&
            n, nd, v, vd0, vd, t, td0, td, &
            arval, arvald0, arvald, arvaldd &
         )
         get_dArdX2 = arvaldd
      end function

      function get_ArnX(var)
         !! Compute mixed second derivatives: ArVn or ArTn
         !! Uses forward-on-backward AD
         character(len=*), intent(in) :: var
         real(pr) :: get_ArnX(size(n))

         call reset_vars

         arvaldb = 1
         select case (var)
         case ("V")
            vd = 1
         case ("T")
            td = 1
         end select

         call self%ar_d_b(&
            n, nb, nd, ndb, &
            v, vb, vd, vdb, &
            t, tb, td, tdb, &
            arval, arvalb, arvald, arvaldb &
         )
         get_ArnX = nb
      end function

   end subroutine

end module