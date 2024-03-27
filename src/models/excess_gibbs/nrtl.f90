MODULE yaeos_models_ge_NRTL
   USE YAEOS_TAPENADE_GE_API, ONLY : gemodeltapenade
   use yaeos_constants, only: pr, R
   IMPLICIT NONE

   TYPE, extends(GeModelTapenade) ::  NRTL
      REAL(pr), ALLOCATABLE :: a(:, :)
      REAL(pr), ALLOCATABLE :: b(:, :)
      REAL(pr), ALLOCATABLE :: c(:, :)
   contains
      procedure :: ge => excess_gibbs
      procedure :: ge_b => excess_gibbs_b
      procedure :: ge_d => excess_gibbs_d
      procedure :: ge_d_b => excess_gibbs_d_b
      procedure :: ge_d_d => excess_gibbs_d_d
   end type NRTL

CONTAINS
   FUNCTION SETUP(a_mat, b_mat, c_mat) RESULT (model)
      IMPLICIT NONE
      REAL(pr), INTENT(IN) :: a_mat(:, :)
      REAL(pr), INTENT(IN) :: b_mat(:, :)
      REAL(pr), INTENT(IN) :: c_mat(:, :)
      type(NRTL) :: model

      model%a = a_mat
      model%b = b_mat
      model%c = c_mat
   end function SETUP

   SUBROUTINE EXCESS_GIBBS_D_D_D(model, n, nd, t, td1, td0, td, ge, ged1&
   &   , ged0, ged0d, ged, gedd0, gedd, geddd)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr), INTENT(IN) :: nd(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr), INTENT(IN) :: td1
      REAL(pr), INTENT(IN) :: td0
      REAL(pr), INTENT(IN) :: td
      REAL(pr), INTENT(OUT) :: ge
      REAL(pr), INTENT(OUT) :: ged1
      REAL(pr), INTENT(OUT) :: ged0
      REAL(pr), INTENT(OUT) :: ged0d
      REAL(pr), INTENT(OUT) :: ged
      REAL(pr), INTENT(OUT) :: gedd0
      REAL(pr), INTENT(OUT) :: gedd
      REAL(pr), INTENT(OUT) :: geddd
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: gd1(SIZE(n), SIZE(n)), taud1(SIZE(n), SIZE(n))
      REAL(pr) :: gd0(SIZE(n), SIZE(n)), taud0(SIZE(n), SIZE(n))
      REAL(pr) :: gd0d(SIZE(n), SIZE(n)), taud0d(SIZE(n), SIZE(n))
      REAL(pr) :: xd(SIZE(n)), gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
      REAL(pr) :: gdd0(SIZE(n), SIZE(n)), taudd0(SIZE(n), SIZE(n))
      REAL(pr) :: gdd(SIZE(n), SIZE(n)), taudd(SIZE(n), SIZE(n))
      REAL(pr) :: gddd(SIZE(n), SIZE(n)), tauddd(SIZE(n), SIZE(n))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d0
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d0d
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d
      REAL(pr), DIMENSION(SIZE(n)) :: arg1dd0
      REAL(pr), DIMENSION(SIZE(n)) :: arg1dd
      REAL(pr), DIMENSION(SIZE(n)) :: arg1ddd
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d1
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d0
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d0d
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d
      REAL(pr), DIMENSION(SIZE(n)) :: arg2dd0
      REAL(pr), DIMENSION(SIZE(n)) :: arg2dd
      REAL(pr), DIMENSION(SIZE(n)) :: arg2ddd
      REAL(pr) :: temp
      REAL(pr) :: tempd0
      REAL(pr) :: tempd
      REAL(pr) :: tempdd
      REAL(pr) :: temp0
      REAL(pr) :: temp0d0
      REAL(pr) :: temp0d
      REAL(pr) :: temp0dd
      REAL(pr) :: temp1
      REAL(pr) :: temp1d0
      REAL(pr) :: temp1d
      REAL(pr) :: temp1dd
      REAL(pr), DIMENSION(SIZE(b, 1), SIZE(b, 2)) :: temp2
      REAL(pr), DIMENSION(SIZE(b, 1), SIZE(b, 2)) :: temp2d
      REAL(pr), DIMENSION(SIZE(n), SIZE(n)) :: temp3
      REAL(pr), DIMENSION(SIZE(n), SIZE(n)) :: temp3d
      REAL(pr), DIMENSION(SIZE(n)) :: temp4
      REAL(pr), DIMENSION(SIZE(n)) :: temp4d
      REAL(pr) :: temp5
      REAL(pr) :: temp5d
      REAL(pr) :: temp6
      REAL(pr) :: temp6d
      REAL(pr), DIMENSION(SIZE(b, 1), SIZE(b, 2)) :: temp7
      REAL(pr), DIMENSION(SIZE(n), SIZE(n)) :: temp8
      REAL(pr), DIMENSION(SIZE(n)) :: temp9
      REAL(pr) :: temp10
      REAL(pr) :: temp11
      temp = SUM(n)
      xd = (nd-n*SUM(nd)/temp)/temp
      x = n/temp
      temp7 = model%b(:, :)*td/(t*t)
      temp2d = -(temp7*2*td1/t)
      temp2 = temp7
      tauddd = td0*2*(temp2d-temp2*td1/t)/t
      taudd = temp2*2*td0/t
      taudd0 = -temp2d
      taud = -temp2
      temp7 = model%b(:, :)*td0/(t*t)
      taud0d = temp7*2*td1/t
      taud0 = -temp7
      taud1 = -(model%b(:, :)*td1/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      temp3d = -(EXP(-(model%c*tau))*model%c*taud1)
      temp3 = EXP(-(model%c*tau))
      temp8 = EXP(-(model%c*tau))
      gddd = -(model%c*(taudd*temp3d+temp3*tauddd-model%c*(temp8*(taud0*&
      &     taudd0+taud*taud0d)-taud*taud0*EXP(-(model%c*tau))*model%c*taud1))&
      &     )
      gdd = -(model%c*(temp3*taudd-model%c*(temp8*(taud*taud0))))
      gdd0 = -(model%c*(taud*temp3d+temp3*taudd0))
      gd = -(model%c*(temp3*taud))
      temp8 = EXP(-(model%c*tau))
      gd0d = -(model%c*(temp8*taud0d-taud0*EXP(-(model%c*tau))*model%c*&
      &     taud1))
      gd0 = -(model%c*(temp8*taud0))
      gd1 = -(EXP(-(model%c*tau))*model%c*taud1)
      g = EXP(-(model%c*tau))
      ge = 0
      ged = 0.0_8
      gedd = 0.0_8
      ged0 = 0.0_8
      gedd0 = 0.0_8
      geddd = 0.0_8
      ged0d = 0.0_8
      ged1 = 0.0_8
      DO i=1,SIZE(n)
         temp4d = xd(:)*taud1(:, i) + x(:)*taudd0(:, i)
         temp4 = xd(:)*tau(:, i) + x(:)*taud(:, i)
         temp9 = xd(:)*taud0(:, i) + x(:)*taudd(:, i)
         arg1ddd(:) = gd0(:, i)*temp4d + temp4*gd0d(:, i) + temp9*gd1(:, i)&
         &       + g(:, i)*(xd(:)*taud0d(:, i)+x(:)*tauddd(:, i)) + x(:)*(taud0(:&
         &       , i)*gdd0(:, i)+gd(:, i)*taud0d(:, i)+gdd(:, i)*taud1(:, i)+tau(&
         &       :, i)*gddd(:, i))
         arg1dd(:) = temp4*gd0(:, i) + g(:, i)*temp9 + x(:)*(gd(:, i)*taud0&
         &       (:, i)+tau(:, i)*gdd(:, i))
         arg1dd0(:) = temp4*gd1(:, i) + g(:, i)*temp4d + x(:)*(gd(:, i)*&
         &       taud1(:, i)+tau(:, i)*gdd0(:, i))
         arg1d(:) = g(:, i)*temp4 + x(:)*(tau(:, i)*gd(:, i))
         arg1d0d(:) = x(:)*(taud0(:, i)*gd1(:, i)+g(:, i)*taud0d(:, i)+gd0(&
         &       :, i)*taud1(:, i)+tau(:, i)*gd0d(:, i))
         arg1d0(:) = x(:)*(g(:, i)*taud0(:, i)+tau(:, i)*gd0(:, i))
         arg1d1(:) = x(:)*(g(:, i)*taud1(:, i)+tau(:, i)*gd1(:, i))
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2ddd(:) = xd(:)*gd0d(:, i) + x(:)*gddd(:, i)
         arg2dd(:) = xd(:)*gd0(:, i) + x(:)*gdd(:, i)
         arg2dd0(:) = xd(:)*gd1(:, i) + x(:)*gdd0(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2d0d(:) = x(:)*gd0d(:, i)
         arg2d0(:) = x(:)*gd0(:, i)
         arg2d1(:) = x(:)*gd1(:, i)
         arg2(:) = x(:)*g(:, i)
         tempdd = SUM(arg2d0d(:))
         tempd = SUM(arg2d0(:))
         tempd0 = SUM(arg2d1(:))
         temp = SUM(arg2(:))
         temp0dd = SUM(arg1d0d(:))
         temp0d = SUM(arg1d0(:))
         temp0d0 = SUM(arg1d1(:))
         temp0 = SUM(arg1(:))
         temp10 = temp0*tempd/temp
         temp11 = (temp0d-temp10)/temp
         temp1dd = x(i)*(temp0dd-(tempd*temp0d0+temp0*tempdd-temp10*tempd0)&
         &       /temp-temp11*tempd0)/temp
         temp1d = x(i)*temp11
         temp1d0 = x(i)*(temp0d0-temp0*tempd0/temp)/temp
         temp1 = x(i)*temp0/temp
         temp5d = SUM(arg2dd0(:))
         temp5 = SUM(arg2d(:))
         temp11 = (xd(i)*temp0+x(i)*SUM(arg1d(:))-temp1*temp5)/temp
         temp6d = (xd(i)*temp0d0+x(i)*SUM(arg1dd0(:))-temp5*temp1d0-temp1*&
         &       temp5d-temp11*tempd0)/temp
         temp6 = temp11
         temp11 = SUM(arg2dd(:))
         temp10 = (xd(i)*temp0d+x(i)*SUM(arg1dd(:))-temp5*temp1d-temp1*&
         &       temp11-temp6*tempd)/temp
         geddd = geddd + (xd(i)*temp0dd+x(i)*SUM(arg1ddd(:))-temp1d*temp5d-&
         &       temp5*temp1dd-temp11*temp1d0-temp1*SUM(arg2ddd(:))-tempd*temp6d-&
         &       temp6*tempdd-temp10*tempd0)/temp
         gedd = gedd + temp10
         gedd0 = gedd0 + temp6d
         ged = ged + temp6
         ged0d = ged0d + temp1dd
         ged0 = ged0 + temp1d
         ged1 = ged1 + temp1d0
         ge = ge + temp1
      END DO
      temp1 = SUM(n)
      temp6 = SUM(nd)
      geddd = r*(temp6*(td0*ged1+ged0*td1+t*ged0d)+temp1*(td*ged0d+td0*&
      &     gedd0+gedd*td1+t*geddd))
      gedd = r*(temp6*(ge*td0+t*ged0)+temp1*(td*ged0+ged*td0+t*gedd))
      gedd0 = r*(temp6*(ge*td1+t*ged1)+temp1*(td*ged1+ged*td1+t*gedd0))
      ged = r*(temp6*(t*ge)+temp1*(td*ge+t*ged))
      ged0d = r*temp1*(td0*ged1+ged0*td1+t*ged0d)
      ged0 = r*temp1*(ge*td0+t*ged0)
      ged1 = r*temp1*(ge*td1+t*ged1)
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D_D_D

   SUBROUTINE EXCESS_GIBBS_D_D(model, n, nd, t, td0, td, ge, ged0, ged, &
   &   gedd)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr), INTENT(IN) :: nd(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr), INTENT(IN) :: td0
      REAL(pr), INTENT(IN) :: td
      REAL(pr), INTENT(OUT) :: ge
      REAL(pr), INTENT(OUT) :: ged0
      REAL(pr), INTENT(OUT) :: ged
      REAL(pr), INTENT(OUT) :: gedd
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: gd0(SIZE(n), SIZE(n)), taud0(SIZE(n), SIZE(n))
      REAL(pr) :: xd(SIZE(n)), gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
      REAL(pr) :: gdd(SIZE(n), SIZE(n)), taudd(SIZE(n), SIZE(n))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d0
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d
      REAL(pr), DIMENSION(SIZE(n)) :: arg1dd
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d0
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d
      REAL(pr), DIMENSION(SIZE(n)) :: arg2dd
      REAL(pr) :: temp
      REAL(pr) :: tempd
      REAL(pr) :: temp0
      REAL(pr) :: temp0d
      REAL(pr) :: temp1
      REAL(pr) :: temp1d
      REAL(pr), DIMENSION(SIZE(b, 1), SIZE(b, 2)) :: temp2
      REAL(pr), DIMENSION(SIZE(n), SIZE(n)) :: temp3
      REAL(pr), DIMENSION(SIZE(n)) :: temp4
      REAL(pr) :: temp5
      REAL(pr) :: temp6
      temp = SUM(n)
      xd = (nd-n*SUM(nd)/temp)/temp
      x = n/temp
      temp2 = model%b(:, :)*td/(t*t)
      taudd = temp2*2*td0/t
      taud = -temp2
      taud0 = -(model%b(:, :)*td0/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      temp3 = EXP(-(model%c*tau))
      gdd = -(model%c*(temp3*taudd-taud*EXP(-(model%c*tau))*model%c*taud0)&
      &     )
      gd = -(model%c*(temp3*taud))
      gd0 = -(EXP(-(model%c*tau))*model%c*taud0)
      g = EXP(-(model%c*tau))
      ge = 0
      ged = 0.0_8
      gedd = 0.0_8
      ged0 = 0.0_8
      DO i=1,SIZE(n)
         temp4 = xd(:)*tau(:, i) + x(:)*taud(:, i)
         arg1dd(:) = temp4*gd0(:, i) + g(:, i)*(xd(:)*taud0(:, i)+x(:)*&
         &       taudd(:, i)) + x(:)*(gd(:, i)*taud0(:, i)+tau(:, i)*gdd(:, i))
         arg1d(:) = g(:, i)*temp4 + x(:)*(tau(:, i)*gd(:, i))
         arg1d0(:) = x(:)*(g(:, i)*taud0(:, i)+tau(:, i)*gd0(:, i))
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2dd(:) = xd(:)*gd0(:, i) + x(:)*gdd(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2d0(:) = x(:)*gd0(:, i)
         arg2(:) = x(:)*g(:, i)
         tempd = SUM(arg2d0(:))
         temp = SUM(arg2(:))
         temp0d = SUM(arg1d0(:))
         temp0 = SUM(arg1(:))
         temp1d = x(i)*(temp0d-temp0*tempd/temp)/temp
         temp1 = x(i)*temp0/temp
         temp5 = SUM(arg2d(:))
         temp6 = (xd(i)*temp0+x(i)*SUM(arg1d(:))-temp1*temp5)/temp
         gedd = gedd + (xd(i)*temp0d+x(i)*SUM(arg1dd(:))-temp5*temp1d-temp1&
         &       *SUM(arg2dd(:))-temp6*tempd)/temp
         ged = ged + temp6
         ged0 = ged0 + temp1d
         ge = ge + temp1
      END DO
      temp1 = SUM(n)
      temp6 = SUM(nd)
      gedd = r*(temp6*(ge*td0+t*ged0)+temp1*(td*ged0+ged*td0+t*gedd))
      ged = r*(temp6*(t*ge)+temp1*(td*ge+t*ged))
      ged0 = r*temp1*(ge*td0+t*ged0)
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D_D

   SUBROUTINE EXCESS_GIBBS_D_B(model, n, nb, nd, ndb, t, tb, td, tdb, ge&
   &   , geb, ged, gedb)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr) :: nb(:)
      REAL(pr), INTENT(IN) :: nd(:)
      REAL(pr) :: ndb(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr) :: tb
      REAL(pr), INTENT(IN) :: td
      REAL(pr) :: tdb
      REAL(pr) :: ge
      REAL(pr) :: geb
      REAL(pr) :: ged
      REAL(pr) :: gedb
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: xb(SIZE(n)), gb(SIZE(n), SIZE(n)), taub(SIZE(n), SIZE(n))
      REAL(pr) :: xd(SIZE(n)), gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
      REAL(pr) :: xdb(SIZE(n)), gdb(SIZE(n), SIZE(n)), taudb(SIZE(n), SIZE(n&
      &   ))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1b
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d
      REAL(pr), DIMENSION(SIZE(n)) :: arg1db
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      REAL(pr), DIMENSION(SIZE(n)) :: arg2b
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d
      REAL(pr), DIMENSION(SIZE(n)) :: arg2db
      REAL(pr) :: temp
      REAL(pr) :: tempb
      REAL(pr) :: temp0
      REAL(pr) :: temp0b
      REAL(pr) :: temp1
      REAL(pr) :: temp1b
      REAL(pr), DIMENSION(SIZE(n, 1)) :: tempb0
      REAL(pr), DIMENSION(SIZE(n, 1)) :: temp2
      REAL(pr) :: temp3
      REAL(pr), DIMENSION(SIZE(n, 1)) :: tempb1
      REAL(pr) :: tempb2
      REAL(pr), DIMENSION(SIZE(n)) :: tempb3
      REAL(pr) :: temp4
      REAL(pr) :: temp5
      REAL(pr) :: tempb4
      REAL(pr) :: tempb5
      INTEGER :: ad_to
      EXTERNAL PUSHREAL8ARRAY
      EXTERNAL PUSHREAL8
      EXTERNAL PUSHINTEGER4
      EXTERNAL POPINTEGER4
      EXTERNAL POPREAL8
      EXTERNAL POPREAL8ARRAY
      INTEGER :: arg10
      REAL(pr) :: result1
      temp = SUM(n)
      xd = (nd-n*SUM(nd)/temp)/temp
      x = n/temp
      taud = -(model%b(:, :)*td/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      gd = -(EXP(-(model%c*tau))*model%c*taud)
      g = EXP(-(model%c*tau))
      ge = 0
      ged = 0.0_8
      DO i=1,SIZE(n)
         arg10 = SIZE(n)
         CALL PUSHREAL8ARRAY(arg1d, arg10)
         arg1d(:) = g(:, i)*(tau(:, i)*xd(:)+x(:)*taud(:, i)) + x(:)*tau(:&
         &       , i)*gd(:, i)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2(:) = x(:)*g(:, i)
         CALL PUSHREAL8(temp)
         temp = SUM(arg2(:))
         CALL PUSHREAL8(temp0)
         temp0 = SUM(arg1(:))
         temp1 = x(i)*temp0/temp
         ged = ged + (temp0*xd(i)+x(i)*SUM(arg1d(:))-temp1*SUM(arg2d(:)))/&
         &       temp
         ge = ge + temp1
      END DO
      CALL PUSHINTEGER4(i - 1)
      temp1 = SUM(n)
      tempb4 = r*geb
      geb = temp1*t*tempb4
      temp1b = t*ge*tempb4
      tb = tb + temp1*ge*tempb4
      tempb4 = r*gedb
      tempb5 = SUM(nd)*tempb4
      ndb = ndb + t*ge*tempb4
      temp1b = temp1b + (ge*td+t*ged)*tempb4
      tempb2 = temp1*tempb4
      gedb = t*tempb2
      geb = geb + td*tempb2 + t*tempb5
      tdb = tdb + ge*tempb2
      tb = tb + ged*tempb2 + ge*tempb5
      nb = nb + temp1b
      taudb = 0.0_8
      taub = 0.0_8
      gb = 0.0_8
      xdb = 0.0_8
      xb = 0.0_8
      gdb = 0.0_8
      CALL POPINTEGER4(ad_to)
      DO i=ad_to,1,-1
         tempb2 = gedb/temp
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         temp5 = SUM(arg2d(:))
         temp1 = x(i)*temp0/temp
         temp1b = geb - temp5*tempb2
         arg1db = 0.0_8
         arg2db = 0.0_8
         temp4 = SUM(arg1d(:))
         temp0b = xd(i)*tempb2
         xdb(i) = xdb(i) + temp0*tempb2
         xb(i) = xb(i) + temp4*tempb2 + temp0*temp1b/temp
         arg1db = x(i)*tempb2
         arg2db = -(temp1*tempb2)
         tempb = -((temp0*xd(i)+x(i)*temp4-temp1*temp5)*tempb2/temp)
         tempb2 = x(i)*temp1b/temp
         temp0b = temp0b + tempb2
         tempb = tempb - temp0*tempb2/temp
         arg1b = 0.0_8
         CALL POPREAL8(temp0)
         arg1b = temp0b
         arg2b = 0.0_8
         CALL POPREAL8(temp)
         arg2b = tempb
         gb(:, i) = gb(:, i) + x*arg2b + xd*arg2db + x*tau(:, i)*arg1b + (&
         &       tau(:, i)*xd+x*taud(:, i))*arg1db
         gdb(:, i) = gdb(:, i) + x*arg2db + x*tau(:, i)*arg1db
         arg10 = SIZE(n)
         CALL POPREAL8ARRAY(arg1d, arg10)
         tempb3 = g(:, i)*arg1db
         xb = xb + g(:, i)*arg2b + gd(:, i)*arg2db + tau(:, i)*g(:, i)*&
         &       arg1b + tau(:, i)*gd(:, i)*arg1db + taud(:, i)*tempb3
         xdb = xdb + g(:, i)*arg2db + tau(:, i)*tempb3
         taub(:, i) = taub(:, i) + x*g(:, i)*arg1b + x*gd(:, i)*arg1db + xd&
         &       *tempb3
         taudb(:, i) = taudb(:, i) + x*tempb3
      END DO
      temp3 = SUM(nd)
      tempb0 = xdb/temp
      tempb1 = -(temp3*tempb0/temp)
      temp2 = n/temp
      result1 = SUM((nd-temp3*temp2)*tempb0)
      tempb = -(SUM(n*xb)/temp**2) - result1/temp - SUM(temp2*tempb1)
      taub = taub + model%c**2*EXP(-(model%c*tau))*taud*gdb - model%c*EXP(&
      &     -(model%c*tau))*gb
      taudb = taudb - EXP(-(model%c*tau))*model%c*gdb
      tempb2 = -(SUM(model%b*taudb)/t**2)
      tb = tb - SUM(model%b*taub)/t**2 - 2*td*tempb2/t
      tdb = tdb + tempb2
      nb = nb + xb/temp + tempb1 + tempb
      ndb = ndb + tempb0 - SUM(temp2*tempb0)
      gedb = 0.0_8
      geb = 0.0_8
   end subroutine EXCESS_GIBBS_D_B

   SUBROUTINE EXCESS_GIBBS_D(model, n, nd, t, td, ge, ged)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr), INTENT(IN) :: nd(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr), INTENT(IN) :: td
      REAL(pr), INTENT(OUT) :: ge
      REAL(pr), INTENT(OUT) :: ged
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: xd(SIZE(n)), gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1d
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      REAL(pr), DIMENSION(SIZE(n)) :: arg2d
      REAL(pr) :: temp
      REAL(pr) :: temp0
      REAL(pr) :: temp1
      temp = SUM(n)
      xd = (nd-n*SUM(nd)/temp)/temp
      x = n/temp
      taud = -(model%b(:, :)*td/t**2)
      tau = model%a(:, :) + model%b(:, :)/t
      gd = -(EXP(-(model%c*tau))*model%c*taud)
      g = EXP(-(model%c*tau))
      ge = 0
      ged = 0.0_8
      DO i=1,SIZE(n)
         arg1d(:) = g(:, i)*(tau(:, i)*xd(:)+x(:)*taud(:, i)) + x(:)*tau(:&
         &       , i)*gd(:, i)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2d(:) = g(:, i)*xd(:) + x(:)*gd(:, i)
         arg2(:) = x(:)*g(:, i)
         temp = SUM(arg2(:))
         temp0 = SUM(arg1(:))
         temp1 = x(i)*temp0/temp
         ged = ged + (temp0*xd(i)+x(i)*SUM(arg1d(:))-temp1*SUM(arg2d(:)))/&
         &       temp
         ge = ge + temp1
      END DO
      temp1 = SUM(n)
      ged = r*(t*ge*SUM(nd)+temp1*(ge*td+t*ged))
      ge = r*(temp1*(t*ge))
   end subroutine EXCESS_GIBBS_D

   SUBROUTINE EXCESS_GIBBS_B(model, n, nb, t, tb, ge, geb)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr) :: nb(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr) :: tb
      REAL(pr) :: ge
      REAL(pr) :: geb
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: xb(SIZE(n)), gb(SIZE(n), SIZE(n)), taub(SIZE(n), SIZE(n))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg1b
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      REAL(pr), DIMENSION(SIZE(n)) :: arg2b
      REAL(pr) :: temp
      REAL(pr) :: tempb
      REAL(pr) :: temp0
      REAL(pr) :: tempb0
      INTEGER :: ad_to
      x = n/SUM(n)
      tau = model%a(:, :) + model%b(:, :)/t
      g = EXP(-(model%c*tau))
      ge = 0
      DO i=1,SIZE(n)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         ge = ge + x(i)*SUM(arg1(:))/SUM(arg2(:))
      END DO
      CALL PUSHINTEGER4(i - 1)
      nb = 0.0_8
      nb = t*ge*r*geb
      tempb0 = SUM(n)*r*geb
      geb = t*tempb0
      tb = ge*tempb0
      taub = 0.0_8
      gb = 0.0_8
      xb = 0.0_8
      CALL POPINTEGER4(ad_to)
      DO i=ad_to,1,-1
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         arg1b = 0.0_8
         arg2b = 0.0_8
         temp = SUM(arg2(:))
         temp0 = SUM(arg1(:))
         tempb = geb/temp
         xb(i) = xb(i) + temp0*tempb
         arg1b = x(i)*tempb
         arg2b = -(x(i)*temp0*tempb/temp)
         xb = xb + g(:, i)*arg2b + tau(:, i)*g(:, i)*arg1b
         gb(:, i) = gb(:, i) + x*arg2b + x*tau(:, i)*arg1b
         taub(:, i) = taub(:, i) + x*g(:, i)*arg1b
      END DO
      taub = taub - model%c*EXP(-(model%c*tau))*gb
      tb = tb - SUM(model%b*taub)/t**2
      temp = SUM(n)
      nb = nb + xb/temp - SUM(n*xb)/temp**2
      geb = 0.0_8
   end subroutine EXCESS_GIBBS_B

   SUBROUTINE EXCESS_GIBBS(model, n, t, ge)
      IMPLICIT NONE
      CLASS(NRTL) :: model
      REAL(pr), INTENT(IN) :: n(:)
      REAL(pr), INTENT(IN) :: t
      REAL(pr), INTENT(OUT) :: ge
      REAL(pr) :: x(SIZE(n)), g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
      REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
      &   n))
      REAL(pr) :: down
      INTEGER :: i, j
      INTRINSIC SUM
      INTRINSIC EXP
      INTRINSIC SIZE
      REAL(pr), DIMENSION(SIZE(n)) :: arg1
      REAL(pr), DIMENSION(SIZE(n)) :: arg2
      x = n/SUM(n)
      tau = model%a(:, :) + model%b(:, :)/t
      g = EXP(-(model%c*tau))
      ge = 0
      DO i=1,SIZE(n)
         arg1(:) = x(:)*tau(:, i)*g(:, i)
         arg2(:) = x(:)*g(:, i)
         ge = ge + x(i)*SUM(arg1(:))/SUM(arg2(:))
      END DO
      ge = SUM(n)*r*t*ge
   end subroutine EXCESS_GIBBS

end module yaeos_models_ge_NRTL

