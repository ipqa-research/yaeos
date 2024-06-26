!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 10 Nov 2023 18:24
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 10 Nov 2023 18:24
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 10 Nov 2023 18:24
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 10 Nov 2023 18:24
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (develop) - 10 Nov 2023 18:24
!
MODULE NRTL_MOD
  USE YAEOS__TAPENADE_GE_API, ONLY : gemodeltapenade
  IMPLICIT NONE
  INTEGER, PARAMETER :: pr=8
  REAL(pr), PARAMETER :: r=0.08314472
! Point to the generated procedures procedures 
! and replate `type :: NRTL` with `type, extends(GeModelTapenade) :: NRTL`
  TYPE NRTL
      REAL(pr), ALLOCATABLE :: a(:, :)
      REAL(pr), ALLOCATABLE :: b(:, :)
      REAL(pr), ALLOCATABLE :: c(:, :)
  END TYPE NRTL

CONTAINS
  FUNCTION SETUP(a_mat, b_mat, c_mat) RESULT (model)
    IMPLICIT NONE
! A setup function for the NRTL model.
! Sets the attributes with provided values.
    REAL(pr), INTENT(IN) :: a_mat(:, :)
    REAL(pr), INTENT(IN) :: b_mat(:, :)
    REAL(pr), INTENT(IN) :: c_mat(:, :)
    CLASS(NRTL) :: model
    model%a = a_mat
    model%b = b_mat
    model%c = c_mat
  END FUNCTION SETUP

!  Differentiation of excess_gibbs_d_d in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ged gedd ged0 ge
!   with respect to varying inputs: t
!   RW status of diff variables: ged:out t:in gedd:out ged0:out
!                ge:out
!  Differentiation of excess_gibbs_d in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ged ge
!   with respect to varying inputs: t
!   RW status of diff variables: ged:out t:in ge:out
!  Differentiation of excess_gibbs in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ge
!   with respect to varying inputs: n t
!   RW status of diff variables: n:in t:in ge:out
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
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
    REAL(pr) :: gd1(SIZE(n), SIZE(n)), taud1(SIZE(n), SIZE(n))
    REAL(pr) :: gd0(SIZE(n), SIZE(n)), taud0(SIZE(n), SIZE(n))
    REAL(pr) :: gd0d(SIZE(n), SIZE(n)), taud0d(SIZE(n), SIZE(n))
    REAL(pr) :: gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
    REAL(pr) :: gdd0(SIZE(n), SIZE(n)), taudd0(SIZE(n), SIZE(n))
    REAL(pr) :: gdd(SIZE(n), SIZE(n)), taudd(SIZE(n), SIZE(n))
    REAL(pr) :: gddd(SIZE(n), SIZE(n)), tauddd(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
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
      temp4d = nd(:)*taud1(:, i) + n(:)*taudd0(:, i)
      temp4 = nd(:)*tau(:, i) + n(:)*taud(:, i)
      temp9 = nd(:)*taud0(:, i) + n(:)*taudd(:, i)
      arg1ddd(:) = gd0(:, i)*temp4d + temp4*gd0d(:, i) + temp9*gd1(:, i)&
&       + g(:, i)*(nd(:)*taud0d(:, i)+n(:)*tauddd(:, i)) + n(:)*(taud0(:&
&       , i)*gdd0(:, i)+gd(:, i)*taud0d(:, i)+gdd(:, i)*taud1(:, i)+tau(&
&       :, i)*gddd(:, i))
      arg1dd(:) = temp4*gd0(:, i) + g(:, i)*temp9 + n(:)*(gd(:, i)*taud0&
&       (:, i)+tau(:, i)*gdd(:, i))
      arg1dd0(:) = temp4*gd1(:, i) + g(:, i)*temp4d + n(:)*(gd(:, i)*&
&       taud1(:, i)+tau(:, i)*gdd0(:, i))
      arg1d(:) = g(:, i)*temp4 + n(:)*(tau(:, i)*gd(:, i))
      arg1d0d(:) = n(:)*(taud0(:, i)*gd1(:, i)+g(:, i)*taud0d(:, i)+gd0(&
&       :, i)*taud1(:, i)+tau(:, i)*gd0d(:, i))
      arg1d0(:) = n(:)*(g(:, i)*taud0(:, i)+tau(:, i)*gd0(:, i))
      arg1d1(:) = n(:)*(g(:, i)*taud1(:, i)+tau(:, i)*gd1(:, i))
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2ddd(:) = nd(:)*gd0d(:, i) + n(:)*gddd(:, i)
      arg2dd(:) = nd(:)*gd0(:, i) + n(:)*gdd(:, i)
      arg2dd0(:) = nd(:)*gd1(:, i) + n(:)*gdd0(:, i)
      arg2d(:) = g(:, i)*nd(:) + n(:)*gd(:, i)
      arg2d0d(:) = n(:)*gd0d(:, i)
      arg2d0(:) = n(:)*gd0(:, i)
      arg2d1(:) = n(:)*gd1(:, i)
      arg2(:) = n(:)*g(:, i)
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
      temp1dd = n(i)*(temp0dd-(tempd*temp0d0+temp0*tempdd-temp10*tempd0)&
&       /temp-temp11*tempd0)/temp
      temp1d = n(i)*temp11
      temp1d0 = n(i)*(temp0d0-temp0*tempd0/temp)/temp
      temp1 = n(i)*temp0/temp
      temp5d = SUM(arg2dd0(:))
      temp5 = SUM(arg2d(:))
      temp11 = (nd(i)*temp0+n(i)*SUM(arg1d(:))-temp1*temp5)/temp
      temp6d = (nd(i)*temp0d0+n(i)*SUM(arg1dd0(:))-temp5*temp1d0-temp1*&
&       temp5d-temp11*tempd0)/temp
      temp6 = temp11
      temp11 = SUM(arg2dd(:))
      temp10 = (nd(i)*temp0d+n(i)*SUM(arg1dd(:))-temp5*temp1d-temp1*&
&       temp11-temp6*tempd)/temp
      geddd = geddd + (nd(i)*temp0dd+n(i)*SUM(arg1ddd(:))-temp1d*temp5d-&
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
    geddd = r*(td*ged0d+td0*gedd0+gedd*td1+t*geddd)
    gedd = r*(td*ged0+ged*td0+t*gedd)
    gedd0 = r*(td*ged1+ged*td1+t*gedd0)
    ged = r*(ge*td+t*ged)
    ged0d = r*(td0*ged1+ged0*td1+t*ged0d)
    ged0 = r*(ge*td0+t*ged0)
    ged1 = r*(ge*td1+t*ged1)
    ge = r*t*ge
  END SUBROUTINE EXCESS_GIBBS_D_D_D

!  Differentiation of excess_gibbs_d in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ged ge
!   with respect to varying inputs: t
!   RW status of diff variables: ged:out t:in ge:out
!  Differentiation of excess_gibbs in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ge
!   with respect to varying inputs: n t
!   RW status of diff variables: n:in t:in ge:out
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
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
    REAL(pr) :: gd0(SIZE(n), SIZE(n)), taud0(SIZE(n), SIZE(n))
    REAL(pr) :: gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
    REAL(pr) :: gdd(SIZE(n), SIZE(n)), taudd(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
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
      temp4 = nd(:)*tau(:, i) + n(:)*taud(:, i)
      arg1dd(:) = temp4*gd0(:, i) + g(:, i)*(nd(:)*taud0(:, i)+n(:)*&
&       taudd(:, i)) + n(:)*(gd(:, i)*taud0(:, i)+tau(:, i)*gdd(:, i))
      arg1d(:) = g(:, i)*temp4 + n(:)*(tau(:, i)*gd(:, i))
      arg1d0(:) = n(:)*(g(:, i)*taud0(:, i)+tau(:, i)*gd0(:, i))
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2dd(:) = nd(:)*gd0(:, i) + n(:)*gdd(:, i)
      arg2d(:) = g(:, i)*nd(:) + n(:)*gd(:, i)
      arg2d0(:) = n(:)*gd0(:, i)
      arg2(:) = n(:)*g(:, i)
      tempd = SUM(arg2d0(:))
      temp = SUM(arg2(:))
      temp0d = SUM(arg1d0(:))
      temp0 = SUM(arg1(:))
      temp1d = n(i)*(temp0d-temp0*tempd/temp)/temp
      temp1 = n(i)*temp0/temp
      temp5 = SUM(arg2d(:))
      temp6 = (nd(i)*temp0+n(i)*SUM(arg1d(:))-temp1*temp5)/temp
      gedd = gedd + (nd(i)*temp0d+n(i)*SUM(arg1dd(:))-temp5*temp1d-temp1&
&       *SUM(arg2dd(:))-temp6*tempd)/temp
      ged = ged + temp6
      ged0 = ged0 + temp1d
      ge = ge + temp1
    END DO
    gedd = r*(td*ged0+ged*td0+t*gedd)
    ged = r*(ge*td+t*ged)
    ged0 = r*(ge*td0+t*ged0)
    ge = r*t*ge
  END SUBROUTINE EXCESS_GIBBS_D_D

!  Differentiation of excess_gibbs_d in reverse (adjoint) mode (with options noISIZE):
!   gradient     of useful results: nd n ged t ge td
!   with respect to varying inputs: nd n ged t ge td
!   RW status of diff variables: nd:incr n:incr ged:in-zero t:incr
!                ge:in-zero td:incr
!  Differentiation of excess_gibbs in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ge
!   with respect to varying inputs: n t
!   RW status of diff variables: n:in t:in ge:out
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
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
    REAL(pr) :: gb(SIZE(n), SIZE(n)), taub(SIZE(n), SIZE(n))
    REAL(pr) :: gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
    REAL(pr) :: gdb(SIZE(n), SIZE(n)), taudb(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
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
    REAL(pr) :: tempb0
    REAL(pr), DIMENSION(SIZE(n)) :: tempb1
    REAL(pr) :: temp2
    REAL(pr) :: temp3
    REAL(pr) :: tempb2
    INTEGER :: ad_to
    EXTERNAL PUSHREAL8ARRAY
    EXTERNAL PUSHREAL8
    EXTERNAL PUSHINTEGER4
    EXTERNAL POPINTEGER4
    EXTERNAL POPREAL8
    EXTERNAL POPREAL8ARRAY
    INTEGER :: arg10
    taud = -(model%b(:, :)*td/t**2)
    tau = model%a(:, :) + model%b(:, :)/t
    gd = -(EXP(-(model%c*tau))*model%c*taud)
    g = EXP(-(model%c*tau))
    ge = 0
    ged = 0.0_8
    DO i=1,SIZE(n)
      arg10 = SIZE(n)
      CALL PUSHREAL8ARRAY(arg1d, arg10)
      arg1d(:) = g(:, i)*(tau(:, i)*nd(:)+n(:)*taud(:, i)) + n(:)*tau(:&
&       , i)*gd(:, i)
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2d(:) = g(:, i)*nd(:) + n(:)*gd(:, i)
      arg2(:) = n(:)*g(:, i)
      CALL PUSHREAL8(temp)
      temp = SUM(arg2(:))
      CALL PUSHREAL8(temp0)
      temp0 = SUM(arg1(:))
      temp1 = n(i)*temp0/temp
      ged = ged + (temp0*nd(i)+n(i)*SUM(arg1d(:))-temp1*SUM(arg2d(:)))/&
&       temp
      ge = ge + temp1
    END DO
    CALL PUSHINTEGER4(i - 1)
    tempb2 = r*gedb
    tb = tb + ge*r*geb + ged*tempb2
    geb = t*r*geb + td*tempb2
    gedb = t*tempb2
    tdb = tdb + ge*tempb2
    taudb = 0.0_8
    taub = 0.0_8
    gb = 0.0_8
    gdb = 0.0_8
    CALL POPINTEGER4(ad_to)
    DO i=ad_to,1,-1
      tempb0 = gedb/temp
      arg2d(:) = g(:, i)*nd(:) + n(:)*gd(:, i)
      temp3 = SUM(arg2d(:))
      temp1 = n(i)*temp0/temp
      temp1b = geb - temp3*tempb0
      arg1db = 0.0_8
      arg2db = 0.0_8
      temp2 = SUM(arg1d(:))
      temp0b = nd(i)*tempb0
      ndb(i) = ndb(i) + temp0*tempb0
      nb(i) = nb(i) + temp2*tempb0 + temp0*temp1b/temp
      arg1db = n(i)*tempb0
      arg2db = -(temp1*tempb0)
      tempb = -((temp0*nd(i)+n(i)*temp2-temp1*temp3)*tempb0/temp)
      tempb0 = n(i)*temp1b/temp
      temp0b = temp0b + tempb0
      tempb = tempb - temp0*tempb0/temp
      arg1b = 0.0_8
      CALL POPREAL8(temp0)
      arg1b = temp0b
      arg2b = 0.0_8
      CALL POPREAL8(temp)
      arg2b = tempb
      gb(:, i) = gb(:, i) + n*arg2b + nd*arg2db + n*tau(:, i)*arg1b + (&
&       tau(:, i)*nd+n*taud(:, i))*arg1db
      gdb(:, i) = gdb(:, i) + n*arg2db + n*tau(:, i)*arg1db
      arg10 = SIZE(n)
      CALL POPREAL8ARRAY(arg1d, arg10)
      tempb1 = g(:, i)*arg1db
      nb = nb + g(:, i)*arg2b + gd(:, i)*arg2db + tau(:, i)*g(:, i)*&
&       arg1b + tau(:, i)*gd(:, i)*arg1db + taud(:, i)*tempb1
      ndb = ndb + g(:, i)*arg2db + tau(:, i)*tempb1
      taub(:, i) = taub(:, i) + n*g(:, i)*arg1b + n*gd(:, i)*arg1db + nd&
&       *tempb1
      taudb(:, i) = taudb(:, i) + n*tempb1
    END DO
    taub = taub + model%c**2*EXP(-(model%c*tau))*taud*gdb - model%c*EXP(&
&     -(model%c*tau))*gb
    taudb = taudb - EXP(-(model%c*tau))*model%c*gdb
    tempb0 = -(SUM(model%b*taudb)/t**2)
    tb = tb - SUM(model%b*taub)/t**2 - 2*td*tempb0/t
    tdb = tdb + tempb0
    gedb = 0.0_8
    geb = 0.0_8
  END SUBROUTINE EXCESS_GIBBS_D_B

!  Differentiation of excess_gibbs in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: ge
!   with respect to varying inputs: n t
!   RW status of diff variables: n:in t:in ge:out
  SUBROUTINE EXCESS_GIBBS_D(model, n, nd, t, td, ge, ged)
    IMPLICIT NONE
    CLASS(NRTL) :: model
    REAL(pr), INTENT(IN) :: n(:)
    REAL(pr), INTENT(IN) :: nd(:)
    REAL(pr), INTENT(IN) :: t
    REAL(pr), INTENT(IN) :: td
    REAL(pr), INTENT(OUT) :: ge
    REAL(pr), INTENT(OUT) :: ged
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
    REAL(pr) :: gd(SIZE(n), SIZE(n)), taud(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
    INTRINSIC SIZE
    REAL(pr), DIMENSION(SIZE(n)) :: arg1
    REAL(pr), DIMENSION(SIZE(n)) :: arg1d
    REAL(pr), DIMENSION(SIZE(n)) :: arg2
    REAL(pr), DIMENSION(SIZE(n)) :: arg2d
    REAL(pr) :: temp
    REAL(pr) :: temp0
    REAL(pr) :: temp1
    taud = -(model%b(:, :)*td/t**2)
    tau = model%a(:, :) + model%b(:, :)/t
    gd = -(EXP(-(model%c*tau))*model%c*taud)
    g = EXP(-(model%c*tau))
    ge = 0
    ged = 0.0_8
    DO i=1,SIZE(n)
      arg1d(:) = g(:, i)*(tau(:, i)*nd(:)+n(:)*taud(:, i)) + n(:)*tau(:&
&       , i)*gd(:, i)
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2d(:) = g(:, i)*nd(:) + n(:)*gd(:, i)
      arg2(:) = n(:)*g(:, i)
      temp = SUM(arg2(:))
      temp0 = SUM(arg1(:))
      temp1 = n(i)*temp0/temp
      ged = ged + (temp0*nd(i)+n(i)*SUM(arg1d(:))-temp1*SUM(arg2d(:)))/&
&       temp
      ge = ge + temp1
    END DO
    ged = r*(ge*td+t*ged)
    ge = r*t*ge
  END SUBROUTINE EXCESS_GIBBS_D

!  Differentiation of excess_gibbs in reverse (adjoint) mode (with options noISIZE):
!   gradient     of useful results: ge
!   with respect to varying inputs: n t ge
!   RW status of diff variables: n:out t:out ge:in-zero
  SUBROUTINE EXCESS_GIBBS_B(model, n, nb, t, tb, ge, geb)
    IMPLICIT NONE
    CLASS(NRTL) :: model
    REAL(pr), INTENT(IN) :: n(:)
    REAL(pr) :: nb(:)
    REAL(pr), INTENT(IN) :: t
    REAL(pr) :: tb
    REAL(pr) :: ge
    REAL(pr) :: geb
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
    REAL(pr) :: gb(SIZE(n), SIZE(n)), taub(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
    INTRINSIC SIZE
    REAL(pr), DIMENSION(SIZE(n)) :: arg1
    REAL(pr), DIMENSION(SIZE(n)) :: arg1b
    REAL(pr), DIMENSION(SIZE(n)) :: arg2
    REAL(pr), DIMENSION(SIZE(n)) :: arg2b
    REAL(pr) :: temp
    REAL(pr) :: tempb
    REAL(pr) :: temp0
    INTEGER :: ad_to
    tau = model%a(:, :) + model%b(:, :)/t
    g = EXP(-(model%c*tau))
    ge = 0
    DO i=1,SIZE(n)
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2(:) = n(:)*g(:, i)
      ge = ge + n(i)*SUM(arg1(:))/SUM(arg2(:))
    END DO
    CALL PUSHINTEGER4(i - 1)
    tb = ge*r*geb
    geb = t*r*geb
    nb = 0.0_8
    taub = 0.0_8
    gb = 0.0_8
    CALL POPINTEGER4(ad_to)
    DO i=ad_to,1,-1
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2(:) = n(:)*g(:, i)
      arg1b = 0.0_8
      arg2b = 0.0_8
      temp = SUM(arg2(:))
      temp0 = SUM(arg1(:))
      tempb = geb/temp
      nb(i) = nb(i) + temp0*tempb
      arg1b = n(i)*tempb
      arg2b = -(n(i)*temp0*tempb/temp)
      nb = nb + g(:, i)*arg2b + tau(:, i)*g(:, i)*arg1b
      gb(:, i) = gb(:, i) + n*arg2b + n*tau(:, i)*arg1b
      taub(:, i) = taub(:, i) + n*g(:, i)*arg1b
    END DO
    taub = taub - model%c*EXP(-(model%c*tau))*gb
    tb = tb - SUM(model%b*taub)/t**2
    geb = 0.0_8
  END SUBROUTINE EXCESS_GIBBS_B

  SUBROUTINE EXCESS_GIBBS(model, n, t, ge)
    IMPLICIT NONE
    CLASS(NRTL) :: model
    REAL(pr), INTENT(IN) :: n(:)
    REAL(pr), INTENT(IN) :: t
    REAL(pr), INTENT(OUT) :: ge
    REAL(pr) :: g(SIZE(n), SIZE(n)), tau(SIZE(n), SIZE(n))
! Due to a tapenade bug, the model attributes are defined as variables,
! though they aren't used. This makes it possible to define sizes for
! generated auxiliar variables.
    REAL(pr) :: a(SIZE(n), SIZE(n)), b(SIZE(n), SIZE(n)), c(SIZE(n), SIZE(&
&   n))
    REAL(pr) :: down
    INTEGER :: i, j
    INTRINSIC EXP
    INTRINSIC SUM
    INTRINSIC SIZE
    REAL(pr), DIMENSION(SIZE(n)) :: arg1
    REAL(pr), DIMENSION(SIZE(n)) :: arg2
    tau = model%a(:, :) + model%b(:, :)/t
    g = EXP(-(model%c*tau))
    ge = 0
    DO i=1,SIZE(n)
      arg1(:) = n(:)*tau(:, i)*g(:, i)
      arg2(:) = n(:)*g(:, i)
      ge = ge + n(i)*SUM(arg1(:))/SUM(arg2(:))
    END DO
    ge = r*t*ge
  END SUBROUTINE EXCESS_GIBBS

END MODULE NRTL_MOD

