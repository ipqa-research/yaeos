module yaeos__models_ge_NRTL
   use yaeos__tapenade_ge_api, only: gemodeltapenade
   use yaeos__tapenade_interfaces
   use yaeos__constants, only: pr, R
   implicit none

   type, extends(GeModelTapenade) ::  NRTL
      !! Non-Random-Two-Liquid model
      !!
      !! \[
      !!    G^E = nRT \cdot \sum_i x_i \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}}
      !! \]
      !!
      !! with:
      !!
      !! \[\tau_{ij} = A_{ij} + \frac{B_{ij}}{T}\]
      !!
      !! \[G_{ij} = exp(-\frac{C_{ij}}{\tau_{ij}})\]
      real(pr), allocatable :: a(:, :) !! A_{ij} matrix
      real(pr), allocatable :: b(:, :) !! B_{ij} matrix
      real(pr), allocatable :: c(:, :) !! C_{ij} matrix
   contains
      procedure :: ge => excess_gibbs
      procedure :: ge_b => excess_gibbs_b
      procedure :: ge_d => excess_gibbs_d
      procedure :: ge_d_b => excess_gibbs_d_b
      procedure :: ge_d_d => excess_gibbs_d_d
   end type NRTL

   interface NRTL
      module procedure :: init
   end interface

   
contains

   type(NRTL) function init(a, b, c)
      real(pr), intent(in) :: a(:, :)
      real(pr), intent(in) :: b(:, :)
      real(pr), intent(in) :: c(:, :)

      init%a = a
      init%b = b
      init%c = c
   end function

   subroutine EXCESS_GIBBS_D_D_D(model, n, nd, t, td1, td0, td, ge, ged1&
   &   , ged0, ged0d, ged, gedd0, gedd, geddd)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr), intent(IN) :: nd(:)
      real(pr), intent(IN) :: t
      real(pr), intent(IN) :: td1
      real(pr), intent(IN) :: td0
      real(pr), intent(IN) :: td
      real(pr), intent(OUT) :: ge
      real(pr), intent(OUT) :: ged1
      real(pr), intent(OUT) :: ged0
      real(pr), intent(OUT) :: ged0d
      real(pr), intent(OUT) :: ged
      real(pr), intent(OUT) :: gedd0
      real(pr), intent(OUT) :: gedd
      real(pr), intent(OUT) :: geddd
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: gd1(size(n), size(n)), taud1(size(n), size(n))
      real(pr) :: gd0(size(n), size(n)), taud0(size(n), size(n))
      real(pr) :: gd0d(size(n), size(n)), taud0d(size(n), size(n))
      real(pr) :: xd(size(n)), gd(size(n), size(n)), taud(size(n), size(n))
      real(pr) :: gdd0(size(n), size(n)), taudd0(size(n), size(n))
      real(pr) :: gdd(size(n), size(n)), taudd(size(n), size(n))
      real(pr) :: gddd(size(n), size(n)), tauddd(size(n), size(n))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg1d1
      real(pr), dimension(size(n)) :: arg1d0
      real(pr), dimension(size(n)) :: arg1d0d
      real(pr), dimension(size(n)) :: arg1d
      real(pr), dimension(size(n)) :: arg1dd0
      real(pr), dimension(size(n)) :: arg1dd
      real(pr), dimension(size(n)) :: arg1ddd
      real(pr), dimension(size(n)) :: arg2
      real(pr), dimension(size(n)) :: arg2d1
      real(pr), dimension(size(n)) :: arg2d0
      real(pr), dimension(size(n)) :: arg2d0d
      real(pr), dimension(size(n)) :: arg2d
      real(pr), dimension(size(n)) :: arg2dd0
      real(pr), dimension(size(n)) :: arg2dd
      real(pr), dimension(size(n)) :: arg2ddd
      real(pr) :: temp
      real(pr) :: tempd0
      real(pr) :: tempd
      real(pr) :: tempdd
      real(pr) :: temp0
      real(pr) :: temp0d0
      real(pr) :: temp0d
      real(pr) :: temp0dd
      real(pr) :: temp1
      real(pr) :: temp1d0
      real(pr) :: temp1d
      real(pr) :: temp1dd
      real(pr), dimension(size(b, 1), size(b, 2)) :: temp2
      real(pr), dimension(size(b, 1), size(b, 2)) :: temp2d
      real(pr), dimension(size(n), size(n)) :: temp3
      real(pr), dimension(size(n), size(n)) :: temp3d
      real(pr), dimension(size(n)) :: temp4
      real(pr), dimension(size(n)) :: temp4d
      real(pr) :: temp5
      real(pr) :: temp5d
      real(pr) :: temp6
      real(pr) :: temp6d
      real(pr), dimension(size(b, 1), size(b, 2)) :: temp7
      real(pr), dimension(size(n), size(n)) :: temp8
      real(pr), dimension(size(n)) :: temp9
      real(pr) :: temp10
      real(pr) :: temp11
      temp = sum(n)
      xd = (nd - n*sum(nd)/temp)/temp
      x = n/temp
      temp7 = model%b(:, :)*td/(t*t)
      temp2d = -(temp7*2*td1/t)
      temp2 = temp7
      tauddd = td0*2*(temp2d - temp2*td1/t)/t
      taudd = temp2*2*td0/t
      taudd0 = -temp2d
      taud = -temp2
      temp7 = model%b(:, :)*td0/(t*t)
      taud0d = temp7*2*td1/t
      taud0 = -temp7
      taud1 = -(model%b(:, :)*td1/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      temp3d = -(exp(-(model%c*tau))*model%c*taud1)
      temp3 = exp(-(model%c*tau))
      temp8 = exp(-(model%c*tau))
      gddd = -(model%c*(taudd*temp3d + temp3*tauddd - model%c*(temp8*(taud0*&
      &     taudd0 + taud*taud0d) - taud*taud0*exp(-(model%c*tau))*model%c*taud1))&
      &     )
      gdd = -(model%c*(temp3*taudd - model%c*(temp8*(taud*taud0))))
      gdd0 = -(model%c*(taud*temp3d + temp3*taudd0))
      gd = -(model%c*(temp3*taud))
      temp8 = exp(-(model%c*tau))
      gd0d = -(model%c*(temp8*taud0d - taud0*exp(-(model%c*tau))*model%c*&
      &     taud1))
      gd0 = -(model%c*(temp8*taud0))
      gd1 = -(exp(-(model%c*tau))*model%c*taud1)
      g = exp(-(model%c*tau))
      ge = 0
      ged = 0.0_pr
      gedd = 0.0_pr
      ged0 = 0.0_pr
      gedd0 = 0.0_pr
      geddd = 0.0_pr
      ged0d = 0.0_pr
      ged1 = 0.0_pr
      do i = 1, size(n)
         temp4d = xd(:)*taud1(:, i) + x(:)*taudd0(:, i)
         temp4 = xd(:)*tau(:, i) + x(:)*taud(:, i)
         temp9 = xd(:)*taud0(:, i) + x(:)*taudd(:, i)
         arg1ddd(:) = gd0(:, i)*temp4d + temp4*gd0d(:, i) + temp9*gd1(:, i)&
         &       + g(:, i)*(xd(:)*taud0d(:, i) + x(:)*tauddd(:, i)) + x(:)*(taud0(:&
         &       , i)*gdd0(:, i) + gd(:, i)*taud0d(:, i) + gdd(:, i)*taud1(:, i) + tau(&
         &       :, i)*gddd(:, i))
         arg1dd(:) = temp4*gd0(:, i) + g(:, i)*temp9 + x(:)*(gd(:, i)*taud0&
         &       (:, i) + tau(:, i)*gdd(:, i))
         arg1dd0(:) = temp4*gd1(:, i) + g(:, i)*temp4d + x(:)*(gd(:, i)*&
         &       taud1(:, i) + tau(:, i)*gdd0(:, i))
         arg1d(:) = g(:, i)*temp4 + x(:)*(tau(:, i)*gd(:, i))
         arg1d0d(:) = x(:)*(taud0(:, i)*gd1(:, i) + g(:, i)*taud0d(:, i) + gd0(&
         &       :, i)*taud1(:, i) + tau(:, i)*gd0d(:, i))
         arg1d0(:) = x(:)*(g(:, i)*taud0(:, i) + tau(:, i)*gd0(:, i))
         arg1d1(:) = x(:)*(g(:, i)*taud1(:, i) + tau(:, i)*gd1(:, i))
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2ddd(:) = xd(:)*gd0d(:, i) + x(:)*gddd(:, i)
         arg2dd(:) = xd(:)*gd0(:, i) + x(:)*gdd(:, i)
         arg2dd0(:) = xd(:)*gd1(:, i) + x(:)*gdd0(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2d0d(:) = x(:)*gd0d(:, i)
         arg2d0(:) = x(:)*gd0(:, i)
         arg2d1(:) = x(:)*gd1(:, i)
         arg2(:) = x(:)*g(:, i)
         tempdd = sum(arg2d0d(:))
         tempd = sum(arg2d0(:))
         tempd0 = sum(arg2d1(:))
         temp = sum(arg2(:))
         temp0dd = sum(arg1d0d(:))
         temp0d = sum(arg1d0(:))
         temp0d0 = sum(arg1d1(:))
         temp0 = sum(arg1(:))
         temp10 = temp0*tempd/temp
         temp11 = (temp0d - temp10)/temp
         temp1dd = x(i)*(temp0dd - (tempd*temp0d0 + temp0*tempdd - temp10*tempd0)&
         &       /temp - temp11*tempd0)/temp
         temp1d = x(i)*temp11
         temp1d0 = x(i)*(temp0d0 - temp0*tempd0/temp)/temp
         temp1 = x(i)*temp0/temp
         temp5d = sum(arg2dd0(:))
         temp5 = sum(arg2d(:))
         temp11 = (xd(i)*temp0 + x(i)*sum(arg1d(:)) - temp1*temp5)/temp
         temp6d = (xd(i)*temp0d0 + x(i)*sum(arg1dd0(:)) - temp5*temp1d0 - temp1*&
         &       temp5d - temp11*tempd0)/temp
         temp6 = temp11
         temp11 = sum(arg2dd(:))
         temp10 = (xd(i)*temp0d + x(i)*sum(arg1dd(:)) - temp5*temp1d - temp1*&
         &       temp11 - temp6*tempd)/temp
         geddd = geddd + (xd(i)*temp0dd + x(i)*sum(arg1ddd(:)) - temp1d*temp5d -&
         &       temp5*temp1dd - temp11*temp1d0 - temp1*sum(arg2ddd(:)) - tempd*temp6d -&
         &       temp6*tempdd - temp10*tempd0)/temp
         gedd = gedd + temp10
         gedd0 = gedd0 + temp6d
         ged = ged + temp6
         ged0d = ged0d + temp1dd
         ged0 = ged0 + temp1d
         ged1 = ged1 + temp1d0
         ge = ge + temp1
      end do
      temp1 = sum(n)
      temp6 = sum(nd)
      geddd = r*(temp6*(td0*ged1 + ged0*td1 + t*ged0d) + temp1*(td*ged0d + td0*&
      &     gedd0 + gedd*td1 + t*geddd))
      gedd = r*(temp6*(ge*td0 + t*ged0) + temp1*(td*ged0 + ged*td0 + t*gedd))
      gedd0 = r*(temp6*(ge*td1 + t*ged1) + temp1*(td*ged1 + ged*td1 + t*gedd0))
      ged = r*(temp6*(t*ge) + temp1*(td*ge + t*ged))
      ged0d = r*temp1*(td0*ged1 + ged0*td1 + t*ged0d)
      ged0 = r*temp1*(ge*td0 + t*ged0)
      ged1 = r*temp1*(ge*td1 + t*ged1)
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D_D_D

   subroutine EXCESS_GIBBS_D_D(model, n, nd, t, td0, td, ge, ged0, ged, &
   &   gedd)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr), intent(IN) :: nd(:)
      real(pr), intent(IN) :: t
      real(pr), intent(IN) :: td0
      real(pr), intent(IN) :: td
      real(pr), intent(OUT) :: ge
      real(pr), intent(OUT) :: ged0
      real(pr), intent(OUT) :: ged
      real(pr), intent(OUT) :: gedd
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: gd0(size(n), size(n)), taud0(size(n), size(n))
      real(pr) :: xd(size(n)), gd(size(n), size(n)), taud(size(n), size(n))
      real(pr) :: gdd(size(n), size(n)), taudd(size(n), size(n))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg1d0
      real(pr), dimension(size(n)) :: arg1d
      real(pr), dimension(size(n)) :: arg1dd
      real(pr), dimension(size(n)) :: arg2
      real(pr), dimension(size(n)) :: arg2d0
      real(pr), dimension(size(n)) :: arg2d
      real(pr), dimension(size(n)) :: arg2dd
      real(pr) :: temp
      real(pr) :: tempd
      real(pr) :: temp0
      real(pr) :: temp0d
      real(pr) :: temp1
      real(pr) :: temp1d
      real(pr), dimension(size(b, 1), size(b, 2)) :: temp2
      real(pr), dimension(size(n), size(n)) :: temp3
      real(pr), dimension(size(n)) :: temp4
      real(pr) :: temp5
      real(pr) :: temp6
      temp = sum(n)
      xd = (nd - n*sum(nd)/temp)/temp
      x = n/temp
      temp2 = model%b(:, :)*td/(t*t)
      taudd = temp2*2*td0/t
      taud = -temp2
      taud0 = -(model%b(:, :)*td0/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      temp3 = exp(-(model%c*tau))
      gdd = -(model%c*(temp3*taudd - taud*exp(-(model%c*tau))*model%c*taud0)&
      &     )
      gd = -(model%c*(temp3*taud))
      gd0 = -(exp(-(model%c*tau))*model%c*taud0)
      g = exp(-(model%c*tau))
      ge = 0
      ged = 0.0_pr
      gedd = 0.0_pr
      ged0 = 0.0_pr
      do i = 1, size(n)
         temp4 = xd(:)*tau(:, i) + x(:)*taud(:, i)
         arg1dd(:) = temp4*gd0(:, i) + g(:, i)*(xd(:)*taud0(:, i) + x(:)*&
         &       taudd(:, i)) + x(:)*(gd(:, i)*taud0(:, i) + tau(:, i)*gdd(:, i))
         arg1d(:) = g(:, i)*temp4 + x(:)*(tau(:, i)*gd(:, i))
         arg1d0(:) = x(:)*(g(:, i)*taud0(:, i) + tau(:, i)*gd0(:, i))
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2dd(:) = xd(:)*gd0(:, i) + x(:)*gdd(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2d0(:) = x(:)*gd0(:, i)
         arg2(:) = x(:)*g(:, i)
         tempd = sum(arg2d0(:))
         temp = sum(arg2(:))
         temp0d = sum(arg1d0(:))
         temp0 = sum(arg1(:))
         temp1d = x(i)*(temp0d - temp0*tempd/temp)/temp
         temp1 = x(i)*temp0/temp
         temp5 = sum(arg2d(:))
         temp6 = (xd(i)*temp0 + x(i)*sum(arg1d(:)) - temp1*temp5)/temp
         gedd = gedd + (xd(i)*temp0d + x(i)*sum(arg1dd(:)) - temp5*temp1d - temp1&
                       &       *sum(arg2dd(:)) - temp6*tempd)/temp
         ged = ged + temp6
         ged0 = ged0 + temp1d
         ge = ge + temp1
      end do
      temp1 = sum(n)
      temp6 = sum(nd)
      gedd = r*(temp6*(ge*td0 + t*ged0) + temp1*(td*ged0 + ged*td0 + t*gedd))
      ged = r*(temp6*(t*ge) + temp1*(td*ge + t*ged))
      ged0 = r*temp1*(ge*td0 + t*ged0)
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D_D

   subroutine EXCESS_GIBBS_D_B(model, n, nb, nd, ndb, t, tb, td, tdb, ge&
   &   , geb, ged, gedb)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr) :: nb(:)
      real(pr), intent(IN) :: nd(:)
      real(pr) :: ndb(:)
      real(pr), intent(IN) :: t
      real(pr) :: tb
      real(pr), intent(IN) :: td
      real(pr) :: tdb
      real(pr) :: ge
      real(pr) :: geb
      real(pr) :: ged
      real(pr) :: gedb
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: xb(size(n)), gb(size(n), size(n)), taub(size(n), size(n))
      real(pr) :: xd(size(n)), gd(size(n), size(n)), taud(size(n), size(n))
      real(pr) :: xdb(size(n)), gdb(size(n), size(n)), taudb(size(n), size(n&
      &   ))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg1b
      real(pr), dimension(size(n)) :: arg1d
      real(pr), dimension(size(n)) :: arg1db
      real(pr), dimension(size(n)) :: arg2
      real(pr), dimension(size(n)) :: arg2b
      real(pr), dimension(size(n)) :: arg2d
      real(pr), dimension(size(n)) :: arg2db
      real(pr) :: temp
      real(pr) :: tempb
      real(pr) :: temp0
      real(pr) :: temp0b
      real(pr) :: temp1
      real(pr) :: temp1b
      real(pr), dimension(size(n, 1)) :: tempb0
      real(pr), dimension(size(n, 1)) :: temp2
      real(pr) :: temp3
      real(pr), dimension(size(n, 1)) :: tempb1
      real(pr) :: tempb2
      real(pr), dimension(size(n)) :: tempb3
      real(pr) :: temp4
      real(pr) :: temp5
      real(pr) :: tempb4
      real(pr) :: tempb5
      integer :: ad_to
      integer :: arg10
      real(pr) :: result1
      temp = sum(n)
      xd = (nd - n*sum(nd)/temp)/temp
      x = n/temp
      taud = -(model%b(:, :)*td/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      gd = -(exp(-(model%c*tau))*model%c*taud)
      g = exp(-(model%c*tau))
      ge = 0
      ged = 0.0_pr
      do i = 1, size(n)
         arg10 = size(n)
         call PUSHREAL8ARRAY(arg1d, arg10)
         arg1d(:) = g(:, i)*(tau(:, i)*xd(:) + x(:)*taud(:, i)) + x(:)*tau(:&
         &       , i)*gd(:, i)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2(:) = x(:)*g(:, i)
         call PUSHREAL8(temp)
         temp = sum(arg2(:))
         call PUSHREAL8(temp0)
         temp0 = sum(arg1(:))
         temp1 = x(i)*temp0/temp
         ged = ged + (temp0*xd(i) + x(i)*sum(arg1d(:)) - temp1*sum(arg2d(:)))/&
         &       temp
         ge = ge + temp1
      end do
      call PUSHINTEGER4(i - 1)
      temp1 = sum(n)
      tempb4 = r*geb
      geb = temp1*t*tempb4
      temp1b = t*ge*tempb4
      tb = tb + temp1*ge*tempb4
      tempb4 = r*gedb
      tempb5 = sum(nd)*tempb4
      ndb = ndb + t*ge*tempb4
      temp1b = temp1b + (ge*td + t*ged)*tempb4
      tempb2 = temp1*tempb4
      gedb = t*tempb2
      geb = geb + td*tempb2 + t*tempb5
      tdb = tdb + ge*tempb2
      tb = tb + ged*tempb2 + ge*tempb5
      nb = nb + temp1b
      taudb = 0.0_pr
      taub = 0.0_pr
      gb = 0.0_pr
      xdb = 0.0_pr
      xb = 0.0_pr
      gdb = 0.0_pr
      call POPINTEGER4(ad_to)
      do i = ad_to, 1, -1
         tempb2 = gedb/temp
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         temp5 = sum(arg2d(:))
         temp1 = x(i)*temp0/temp
         temp1b = geb - temp5*tempb2
         arg1db = 0.0_pr
         arg2db = 0.0_pr
         temp4 = sum(arg1d(:))
         temp0b = xd(i)*tempb2
         xdb(i) = xdb(i) + temp0*tempb2
         xb(i) = xb(i) + temp4*tempb2 + temp0*temp1b/temp
         arg1db = x(i)*tempb2
         arg2db = -(temp1*tempb2)
         tempb = -((temp0*xd(i) + x(i)*temp4 - temp1*temp5)*tempb2/temp)
         tempb2 = x(i)*temp1b/temp
         temp0b = temp0b + tempb2
         tempb = tempb - temp0*tempb2/temp
         arg1b = 0.0_pr
         call POPREAL8(temp0)
         arg1b = temp0b
         arg2b = 0.0_pr
         call POPREAL8(temp)
         arg2b = tempb
         gb(:, i) = gb(:, i) + x*arg2b + xd*arg2db + x*tau(:, i)*arg1b + (&
         &       tau(:, i)*xd + x*taud(:, i))*arg1db
         gdb(:, i) = gdb(:, i) + x*arg2db + x*tau(:, i)*arg1db
         arg10 = size(n)
         call POPREAL8ARRAY(arg1d, arg10)
         tempb3 = g(:, i)*arg1db
         xb = xb + g(:, i)*arg2b + gd(:, i)*arg2db + tau(:, i)*g(:, i)*&
         &       arg1b + tau(:, i)*gd(:, i)*arg1db + taud(:, i)*tempb3
         xdb = xdb + g(:, i)*arg2db + tau(:, i)*tempb3
         taub(:, i) = taub(:, i) + x*g(:, i)*arg1b + x*gd(:, i)*arg1db + xd&
                     &       *tempb3
         taudb(:, i) = taudb(:, i) + x*tempb3
      end do
      temp3 = sum(nd)
      tempb0 = xdb/temp
      tempb1 = -(temp3*tempb0/temp)
      temp2 = n/temp
      result1 = sum((nd - temp3*temp2)*tempb0)
      tempb = -(sum(n*xb)/temp**2) - result1/temp - sum(temp2*tempb1)
      taub = taub + model%c**2*exp(-(model%c*tau))*taud*gdb - model%c*exp(&
      &     -(model%c*tau))*gb
      taudb = taudb - exp(-(model%c*tau))*model%c*gdb
      tempb2 = -(sum(model%b*taudb)/t**2)
      tb = tb - sum(model%b*taub)/t**2 - 2*td*tempb2/t
      tdb = tdb + tempb2
      nb = nb + xb/temp + tempb1 + tempb
      ndb = ndb + tempb0 - sum(temp2*tempb0)
      gedb = 0.0_pr
      geb = 0.0_pr
   end subroutine EXCESS_GIBBS_D_B

   subroutine EXCESS_GIBBS_D(model, n, nd, t, td, ge, ged)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr), intent(IN) :: nd(:)
      real(pr), intent(IN) :: t
      real(pr), intent(IN) :: td
      real(pr), intent(OUT) :: ge
      real(pr), intent(OUT) :: ged
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: xd(size(n)), gd(size(n), size(n)), taud(size(n), size(n))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg1d
      real(pr), dimension(size(n)) :: arg2
      real(pr), dimension(size(n)) :: arg2d
      real(pr) :: temp
      real(pr) :: temp0
      real(pr) :: temp1
      temp = sum(n)
      xd = (nd - n*sum(nd)/temp)/temp
      x = n/temp
      taud = -(model%b(:, :)*td/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      gd = -(exp(-(model%c*tau))*model%c*taud)
      g = exp(-(model%c*tau))
      ge = 0
      ged = 0.0_pr
      do i = 1, size(n)
         arg1d(:) = g(:, i)*(tau(:, i)*xd(:) + x(:)*taud(:, i)) + x(:)*tau(:&
         &       , i)*gd(:, i)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2(:) = x(:)*g(:, i)
         temp = sum(arg2(:))
         temp0 = sum(arg1(:))
         temp1 = x(i)*temp0/temp
         ged = ged + (temp0*xd(i) + x(i)*sum(arg1d(:)) - temp1*sum(arg2d(:)))/&
         &       temp
         ge = ge + temp1
      end do
      temp1 = sum(n)
      ged = r*(t*ge*sum(nd) + temp1*(ge*td + t*ged))
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D

   subroutine EXCESS_GIBBS_B(model, n, nb, t, tb, ge, geb)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr) :: nb(:)
      real(pr), intent(IN) :: t
      real(pr) :: tb
      real(pr) :: ge
      real(pr) :: geb
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: xb(size(n)), gb(size(n), size(n)), taub(size(n), size(n))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg1b
      real(pr), dimension(size(n)) :: arg2
      real(pr), dimension(size(n)) :: arg2b
      real(pr) :: temp
      real(pr) :: tempb
      real(pr) :: temp0
      real(pr) :: tempb0
      integer :: ad_to
      x = n/sum(n)
      tau = model%a(:, :) + model%b(:, :)/t
      g = exp(-(model%c*tau))
      ge = 0
      do i = 1, size(n)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         ge = ge + x(i)*sum(arg1(:))/sum(arg2(:))
      end do
      call PUSHINTEGER4(i - 1)
      nb = 0.0_pr
      nb = t*ge*r*geb
      tempb0 = sum(n)*r*geb
      geb = t*tempb0
      tb = ge*tempb0
      taub = 0.0_pr
      gb = 0.0_pr
      xb = 0.0_pr
      call POPINTEGER4(ad_to)
      do i = ad_to, 1, -1
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         arg1b = 0.0_pr
         arg2b = 0.0_pr
         temp = sum(arg2(:))
         temp0 = sum(arg1(:))
         tempb = geb/temp
         xb(i) = xb(i) + temp0*tempb
         arg1b = x(i)*tempb
         arg2b = -(x(i)*temp0*tempb/temp)
         xb = xb + g(:, i)*arg2b + tau(:, i)*g(:, i)*arg1b
         gb(:, i) = gb(:, i) + x*arg2b + x*tau(:, i)*arg1b
         taub(:, i) = taub(:, i) + x*g(:, i)*arg1b
      end do
      taub = taub - model%c*exp(-(model%c*tau))*gb
      tb = tb - sum(model%b*taub)/t**2
      temp = sum(n)
      nb = nb + xb/temp - sum(n*xb)/temp**2
      geb = 0.0_pr
   end subroutine EXCESS_GIBBS_B

   subroutine EXCESS_GIBBS(model, n, t, ge)
      implicit none
      class(NRTL) :: model
      real(pr), intent(IN) :: n(:)
      real(pr), intent(IN) :: t
      real(pr), intent(OUT) :: ge
      real(pr) :: x(size(n)), g(size(n), size(n)), tau(size(n), size(n))
      real(pr) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(&
      &   n))
      real(pr) :: down
      integer :: i, j
      intrinsic SUM
      intrinsic EXP
      intrinsic SIZE
      real(pr), dimension(size(n)) :: arg1
      real(pr), dimension(size(n)) :: arg2
      x = n/sum(n)
      tau = model%a(:, :) + model%b(:, :)/t
      g = exp(-(model%c*tau))
      ge = 0
      do i = 1, size(n)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         ge = ge + x(i)*sum(arg1(:))/sum(arg2(:))
      end do
      ge = sum(n)*r*t*ge
   end subroutine EXCESS_GIBBS
end module yaeos__models_ge_NRTL

