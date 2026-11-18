module yaeos__tapenade_ge_api
   !! # Tapenade excess-Gibbs API
   !! Abstract interface and adapter for source-to-source differentiated GE models.
   !!
   !! # Description
   !! `GeModelTapenade` adapts generated primal and derivative routines to the
   !! `GeModel` excess-Gibbs and activity-coefficient methods. Concrete models
   !! implement the deferred procedures declared in this module.
   !!
   !! # Examples
   !! Extend `GeModelTapenade` and call the inherited GE property methods on
   !! the resulting model.
   !!
   !! # References
   !! Tapenade documentation: https://tapenade.gitlabpages.inria.fr/
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   implicit none

   private

   public :: GeModelTapenade
   type, abstract, extends(GeModel) :: GeModelTapenade
      !! # GeModelTapenade
      !! Base class for excess-Gibbs models with Tapenade-generated derivatives.
      !!
      !! # Description
      !! A concrete extension supplies primal, tangent, adjoint, and mixed
      !! derivative routines. The type adapts them to the standard `GeModel`
      !! excess-Gibbs interface.
      !!
      !! # Examples
      !! See the implementation template in `tools/tapenade_diff/ge_models/`.
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
         !! # tapenade_ge
         !! Primal excess-Gibbs callback required by `GeModelTapenade`.
         import GeModelTapenade, pr
         class(GeModelTapenade) :: model !! Concrete activity model.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(out) :: ge !! Excess Gibbs energy.
      end subroutine

      subroutine tapenade_ge_d(model, n, nd, t, td, ge, ged)
         !! # tapenade_ge_d
         !! Forward directional derivative callback for excess Gibbs energy.
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model !! Concrete activity model.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(in) :: nd(:) !! Tangent seed for component mole amounts (nd).
         real(pr), intent(in) :: td !! Tangent seed for temperature (td).
         real(pr), intent(out) :: ge !! Primal excess Gibbs energy.
         real(pr), intent(out) :: ged !! Forward directional derivative of excess Gibbs energy.
      end subroutine

      subroutine tapenade_ge_b(model, n, nb, t, tb, ge, geb)
         !! # tapenade_ge_b
         !! Reverse-mode adjoint callback for excess Gibbs energy.
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model !! Concrete activity model.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr) :: geb !! Incoming adjoint seed on the excess Gibbs output.
         real(pr) :: nb(:) !! Adjoint value for component mole amounts.
         real(pr) :: tb !! Adjoint value for temperature.
         real(pr) :: ge !! Primal excess Gibbs energy.
      end subroutine

      subroutine tapenade_ge_d_b(model, &
         n, nb, nd, ndb, t, tb, td, tdb, &
         ge, geb, ged, gedb)
         !! # tapenade_ge_d_b
         !! Mixed forward/reverse derivative callback for GE Hessian components.
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model !! Concrete activity model.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr) :: ge

         real(pr), intent(in) :: nd(:) !! Tangent seed for component mole amounts (nd).
         real(pr), intent(in) :: td !! Tangent seed for temperature (td).
         real(pr) :: ged !! Forward directional energy derivative.

         real(pr) :: nb(:) !! Adjoint value for component mole amounts.
         real(pr) :: tb !! Adjoint value for temperature.
         real(pr) :: geb !! Adjoint seed for the primal output.

         real(pr) :: ndb(:) !! Mixed tangent-adjoint value for component mole amounts.
         real(pr) :: tdb !! Mixed tangent-adjoint value for temperature.
         real(pr) :: gedb !! Mixed tangent-adjoint output for GE.
      end subroutine

      subroutine tapenade_ge_d_d(model, n, nd, t, td0, td, ge, ged0, ged, gedd)
         !! # tapenade_ge_d_d
         !! Nested forward-mode callback for second directional GE derivatives.
         import pr, GeModelTapenade
         class(GeModelTapenade) :: model !! Concrete activity model.
         real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
         real(pr), intent(in) :: t !! Temperature [K].
         real(pr), intent(in) :: td0 !! First temperature tangent seed.
         real(pr), intent(in) :: nd(:) !! Tangent seed for component mole amounts (nd).
         real(pr), intent(in) :: td !! Tangent seed for temperature (td).
         real(pr), intent(out) :: ge !! Primal excess Gibbs energy.
         real(pr), intent(out) :: ged0 !! First nested tangent of excess Gibbs energy.
         real(pr), intent(out) :: ged !! Forward directional derivative of excess Gibbs energy.
         real(pr), intent(out) :: gedd !! Mixed second tangent of excess Gibbs energy.
      end subroutine
   end interface

contains

   subroutine excess_gibbs(&
      self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2 &
      )
      !! Excess Gibbs model generic interface
      class(GeModelTapenade), intent(in) :: self !! Model instance.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs energy [bar L].
      real(pr), optional, intent(out) :: GeT !! \(\partial G^E/\partial T\) [bar L/K].
      real(pr), optional, intent(out) :: GeT2 !! \(\partial^2 G^E/\partial T^2\) [bar L/K^2].
      real(pr), optional, dimension(size(n)), intent(out) :: Gen !! \(\partial G^E/\partial n_i\) [bar L/mol].
      real(pr), optional, dimension(size(n)), intent(out) :: GeTn !! \(\partial^2 G^E/\partial T\partial n_i\) [bar L/(K mol)].
      real(pr), optional, intent(out) :: Gen2(size(n), size(n)) !! \(\partial^2 G^E/\partial n_i\partial n_j\) [bar L/mol^2].

      real(pr) :: nb(size(n))
      real(pr) :: nd(size(n))
      real(pr) :: ndb(size(n))
      real(pr) :: tb
      real(pr) :: td
      real(pr) :: tdb
      real(pr) :: td0
      real(pr) :: geb
      real(pr) :: ged
      real(pr) :: gedb
      real(pr) :: ged0
      real(pr) :: gedd

      integer :: i !! Loop/index variable.
      integer :: nc !! Number of components.

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
      else if (present(Gen)) then
            call reset_vars
            geb = 1
            call self%ge_b(n, nb, t, tb, ge, geb)
            Gen = nb
            if (present(GeT)) GeT = tb
      else
         if (present(GeT)) GeT = get_dGedT()
      end if

      if (present(GeTn)) GeTn = get_GenT()
      if (present(GeT2)) GeT2 = get_dGedT2()

      call reset_vars
      if (present(Ge)) call self%ge(n, t, ge)

   contains
      subroutine reset_vars
         !! Reset Tapenade adjoint/tangent work values before a new seed pass.
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
         !! # get_dGedT2
         !!
         !! # Description
         !! Extract the second temperature derivative from the nested tangent routine.
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

      function get_dGedT()
         !! # get_dGedT
         !!
         !! # Description
         !! Extract the first temperature derivative from the tangent routine.
         real(pr) :: get_dGedT
         call reset_vars

         td = 1

         call self%ge_d(&
            n, nd, t, td, ge, ged &
            )
         get_dGedT = ged
      end function

      function get_GenT()
         !! # get_GenT
         !!
         !! # Description
         !! Extract mixed mole-number/temperature derivatives from the mixed mode.
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
