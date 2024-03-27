module yaeos_tapenade_ge_api
   !! Module that wraps tapenade generated routines to calculate !
   !! Ge and derivatives.
   use yaeos_constants, only: pr, R
   use yaeos_models_ge, only: GeModel
   implicit none

   private

   public :: GeModelTapenade
   type, abstract, extends(GeModel) :: GeModelTapenade
   contains
      procedure(tapenade_ge), deferred  :: ge
      procedure(tapenade_ge_d), deferred :: ge_d
      procedure(tapenade_ge_b), deferred :: ge_b
      procedure(tapenade_ge_d_b), deferred :: ge_d_b
      procedure(tapenade_ge_d_d), deferred :: ge_d_d
      procedure :: excess_gibbs => excess_gibbs
   end type

   abstract interface
      subroutine tapenade_ge(model, n, t, ge)
         import GeModelTapenade, pr
         class(GeModelTapenade) :: model
         real(pr), intent(in) :: n(:), t
         real(pr), intent(out) :: ge
      end subroutine

      subroutine tapenade_ge_d(model, n, nd, t, td, ge, ged)
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model
         real(pr), intent(in) :: n(:), t
         real(pr), intent(in) :: nd(:), td
         real(pr), intent(out) :: ge, ged
      end subroutine

      subroutine tapenade_ge_b(model, n, nb, t, tb, ge, geb)
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model
         real(pr), intent(in) :: n(:), t
         real(pr) :: geb
         real(pr) :: nb(:), tb
         real(pr) :: ge
      end subroutine

      subroutine tapenade_ge_d_b(model, &
         n, nb, nd, ndb, t, tb, td, tdb, &
         ge, geb, ged, gedb)
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model
         real(pr), intent(in) :: n(:), t
         real(pr) :: ge

         real(pr), intent(in) :: nd(:), td
         real(pr) :: ged

         real(pr) :: nb(:), tb
         real(pr) :: geb

         real(pr) :: ndb(:), tdb
         real(pr) :: gedb
      end subroutine

      subroutine tapenade_ge_d_d(model, n, nd, t, td0, td, ge, ged0, ged, gedd)
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model
         real(pr), intent(in) :: n(:), t
         real(pr), intent(in) :: td0
         real(pr), intent(in) :: nd(:), td
         real(pr), intent(out) :: ge, ged0, ged, gedd
      end subroutine
   end interface

contains

   subroutine excess_gibbs(&
      self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2 &
      )
      !! Excess Gibbs model generic interface
      class(GeModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: t
      real(pr), optional, intent(out) :: Ge, GeT, GeT2
      real(pr), optional, dimension(size(n)), intent(out) :: Gen, GeTn
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      real(pr) :: nb(size(n)), nd(size(n)), ndb(size(n))
      real(pr) :: tb, td, tdb, td0
      real(pr) :: geb, ged, gedb, ged0, gedd

      integer :: i, nc

      nc = size(n)

      if (present(Gen2)) then
         do i=1, nc
            call reset_vars

            gedb = 1
            if (i <= nc) then
               nd(i) = 1
            end if

            call self%ge_d_b(&
               n, nb, nd, ndb, &
               t, tb, td, tdb, &
               ge, geb, ged, gedb &
            )

            Gen2(i, :) = nb
         end do

         if(present(Gen)) Gen = ndb
         if(present(GeT)) GeT = tdb
      else
         if (present(Gen)) then
            call reset_vars
            geb = 1
            call self%ge_b(n, nb, t, tb, ge, geb)
            Gen = nb
            if (present(GeT)) GeT = tb
         end if
      end if

      if (present(GeTn)) GeTn = get_GenT()
      if (present(GeT2)) GeT2 = get_dGedT2()

      if (present(Ge)) call self%ge(n, t, ge)

   contains
      subroutine reset_vars
         nb=0
         nd=0
         ndb=0

         tb=0
         td=0
         td0 = 0
         tdb=0

         ge = 0
         geb = 0
         ged = 0
         ged0 = 0
         gedb = 0
      end subroutine

      function get_dGedT2()
         real(pr) :: get_dGedT2

         call reset_vars

         td = 1
         td0 = 1

         call self%ge_d_d(&
            n, nd, t, td0, td, &
            ge, ged0, ged, gedd &
         )
         get_dGedT2 = gedd
      end function

      function get_GenT()
         real(pr) :: get_GenT(size(n))
         call reset_vars
         
         gedb = 1
         td = 1

         call self%ge_d_b(&
            n, nb, nd, ndb, &
            t, tb, td, tdb, &
            ge, geb, ged, gedb &
            )
         get_GenT = nb
      end function
   end subroutine
end module
