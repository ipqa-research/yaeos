module yaeos_tapenade_ar_api
   !! Module that wraps tapenade generated routines to calculate !
   !! Ar and derivatives.
   use yaeos_constants, only: pr
   use yaeos_models_ar, only: ArModel
   implicit none

   private

   public :: ArModelTapenade

   type, extends(ArModel) :: ArModelTapenade
      procedure(tapenade_ar), pointer, nopass :: ar
      procedure(tapenade_ar_d), pointer, nopass :: ar_d
      procedure(tapenade_ar_b), pointer, nopass :: ar_b
      procedure(tapenade_ar_d_b), pointer, nopass :: ar_d_b
      procedure(tapenade_ar_d_d), pointer, nopass :: ar_d_d
      procedure(tapenade_v0), pointer, nopass :: v0
   contains
      procedure :: residual_helmholtz => residual_helmholtz
      procedure :: get_v0 => get_v0
   end type

   abstract interface
      subroutine tapenade_ar(n, v, t, arval)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(out) :: arval
      end subroutine

      subroutine tapenade_ar_d(n, nd, v, vd, t, td, arval, arvald)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald
      end subroutine

      subroutine tapenade_ar_b(n, nb, v, vb, t, tb, arval, arvalb)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arvalb
         real(pr) :: nb(:), vb, tb
         real(pr) :: arval
      end subroutine

      subroutine tapenade_ar_d_b(&
         n, nb, nd, ndb, v, vb, vd, vdb, t, tb, td, tdb, &
         arval, arvalb, arvald, arvaldb)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr) :: arval

         real(pr), intent(in) :: nd(:), vd, td
         real(pr) :: arvald

         real(pr) :: nb(:), vb, tb
         real(pr) :: arvalb

         real(pr) :: ndb(:), vdb, tdb
         real(pr) :: arvaldb
      end subroutine

      subroutine tapenade_ar_d_d(n, nd, v, vd0, vd, t, td0, td, &
         arval, arvald0, arvald, arvaldd)
         import pr
         real(pr), intent(in) :: n(:), v, t
         real(pr), intent(in) :: vd0, td0
         real(pr), intent(in) :: nd(:), vd, td
         real(pr), intent(out) :: arval, arvald0, arvald, arvaldd
      end subroutine

      function tapenade_v0(n, p, t)
         import pr
         real(pr), intent(in) :: n(:), p, t
         real(pr) :: tapenade_v0
      end function
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

      real(pr) :: df(size(n) + 2), df2(size(n) + 2, size(n) + 2)

      real(pr) :: nb(size(n)), nd(size(n)), ndb(size(n))
      real(pr) :: vb, vd, vdb, vd0
      real(pr) :: tb, td, tdb, td0
      real(pr) :: arval, arvalb, arvald, arvaldb, arvald0, arvaldd

      integer :: i, nc

      nc = size(n)

      if (present(Arn2)) then
         do i=1, nc
            call reset_vars

            arvaldb = 1
            if (i <= nc) then
               nd(i) = 1
            end if

            call self%ar_d_b(&
               n, nb, nd, ndb, &
               v, vb, vd, vdb, &
               t, tb, td, tdb, &
               arval, arvalb, arvald, arvaldb &
            )

            Arn2(i, :) = nb
         end do

         if(present(Arn)) Arn = ndb
         if(present(ArV)) ArV = vdb
         if(present(ArT)) ArT = tdb
      else
         if (present(Arn)) then
            call reset_vars
            arvalb = 1
            call self%ar_b(n, nb, v, vb, t, tb, arval, arvalb)
            Arn = nb
            if (present(ArT)) ArT = tb
            if (present(ArV)) ArV = vb
         end if
      end if

      if (present(ArTn)) ArTn = get_ArnX("T")
      if (present(ArVn)) ArVn = get_ArnX("V")
      if (present(ArTV)) ArTV = get_dArdX2("TV")
      if (present(ArT2)) ArT2 = get_dArdX2("T2")
      if (present(ArV2)) ArV2 = get_dArdX2("V2")

      if (present(Ar)) Ar = arval

   contains
      subroutine reset_vars
         nb=0
         nd=0
         ndb=0

         vb=0
         vd=0
         vd0 = 0
         vdb=0

         tb=0
         td=0
         td0 = 0
         tdb=0

         arval = 0
         arvalb = 0
         arvald = 0
         arvald0 = 0
         arvaldb = 0
      end subroutine

      function get_dArdX2(var)
         character(len=*), intent(in) :: var
         real(pr) :: get_dArdX2

         call reset_vars

         select case(var)
          case("TV")
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
         character(len=*), intent(in) :: var
         real(pr) :: get_ArnX(size(n))
         call reset_vars
         
         arvaldb = 1
         select case(var)
          case("V")
            vd = 1
          case("T")
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

   function get_v0(self, n, p, t)
      class(ArModelTapenade), intent(in) :: self
      real(pr), intent(in) :: n(:), p, t
      real(pr) :: get_v0

      get_v0 = self%v0(n, p, t)
   end function
end module
