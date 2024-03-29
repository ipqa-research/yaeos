MODULE autodiff_tapenade_pr76
  USE YAEOS_TAPENADE_AR_API, ONLY : armodeltapenade
  IMPLICIT NONE
  REAL(8), ALLOCATABLE :: kij(:, :), lij(:, :)
  REAL(8), ALLOCATABLE :: ac(:), b(:), k(:)
  REAL(8), ALLOCATABLE :: tc(:), pc(:), w(:)
  REAL(8), PARAMETER :: r=0.08314472
  REAL(8), PARAMETER :: del1=1.+SQRT(2.)
  REAL(8), PARAMETER :: del2=1.-SQRT(2.)
  TYPE(ARMODELTAPENADE) :: model

CONTAINS
  SUBROUTINE SETUP(tc_in, pc_in, w_in, kij_in, lij_in)
    IMPLICIT NONE
    REAL(8) :: tc_in(:)
    REAL(8) :: pc_in(:)
    REAL(8) :: w_in(:)
    REAL(8) :: kij_in(:, :)
    REAL(8) :: lij_in(:, :)
    
    tc = tc_in
    pc = pc_in
    w = w_in
    ac = 0.45723553*r**2*tc**2/pc
    b = 0.07779607*r*tc/pc
    k = 0.37464 + 1.54226*w - 0.26993*w**2
    kij = kij_in
    lij = lij_in
    model%ar => ar
    model%ar_d => ar_d
    model%ar_b => ar_b
    model%ar_d_b => ar_d_b
    model%ar_d_d => ar_d_d
  END SUBROUTINE SETUP

  SUBROUTINE AR_D_D_D(n, nd, v, vd1, vd0, vd, t, td1, td0, td, arval, &
&   arvald1, arvald0, arvald0d, arvald, arvaldd0, arvaldd, arvalddd)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8), INTENT(IN) :: vd1, td1
    REAL(8), INTENT(IN) :: vd0, td0
    REAL(8), INTENT(IN) :: nd(:), vd, td
    REAL(8), INTENT(OUT) :: arval
    REAL(8), INTENT(OUT) :: arvald1
    REAL(8), INTENT(OUT) :: arvald0
    REAL(8), INTENT(OUT) :: arvald0d
    REAL(8), INTENT(OUT) :: arvald
    REAL(8), INTENT(OUT) :: arvaldd0
    REAL(8), INTENT(OUT) :: arvaldd
    REAL(8), INTENT(OUT) :: arvalddd
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: amixd1, ad1(SIZE(n))
    REAL(8) :: amixd0, ad0(SIZE(n))
    REAL(8) :: amixd0d, ad0d(SIZE(n))
    REAL(8) :: amixd, ad(SIZE(n)), nijd
    REAL(8) :: amixdd0, add0(SIZE(n))
    REAL(8) :: amixdd, add(SIZE(n))
    REAL(8) :: amixddd, addd(SIZE(n))
    REAL(8) :: bmix
    REAL(8) :: bmixd
    REAL(8) :: b_v
    REAL(8) :: b_vd1
    REAL(8) :: b_vd0
    REAL(8) :: b_vd0d
    REAL(8) :: b_vd
    REAL(8) :: b_vdd0
    REAL(8) :: b_vdd
    REAL(8) :: b_vddd
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d0d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1dd0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1dd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1ddd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d0d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2dd0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2dd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2ddd
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n)) :: arg10d1
    REAL(8), DIMENSION(SIZE(n)) :: arg10d0
    REAL(8), DIMENSION(SIZE(n)) :: arg10d0d
    REAL(8), DIMENSION(SIZE(n)) :: arg10d
    REAL(8), DIMENSION(SIZE(n)) :: arg10dd0
    REAL(8), DIMENSION(SIZE(n)) :: arg10dd
    REAL(8), DIMENSION(SIZE(n)) :: arg10ddd
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11d
    REAL(8) :: arg12
    REAL(8) :: arg12d1
    REAL(8) :: arg12d0
    REAL(8) :: arg12d0d
    REAL(8) :: arg12d
    REAL(8) :: arg12dd0
    REAL(8) :: arg12dd
    REAL(8) :: arg12ddd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempd0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempdd
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0d0
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0d
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0dd
    REAL(8) :: temp1
    REAL(8) :: temp1d0
    REAL(8) :: temp1d
    REAL(8) :: temp1dd
    REAL(8) :: temp2
    REAL(8) :: temp3
    REAL(8) :: temp3d0
    REAL(8) :: temp3d
    REAL(8) :: temp4
    REAL(8) :: temp4d0
    REAL(8) :: temp4d
    REAL(8) :: temp4dd
    REAL(8) :: temp5
    REAL(8) :: temp5d0
    REAL(8) :: temp5d
    REAL(8) :: temp5dd
    REAL(8) :: temp6
    REAL(8) :: temp6d0
    REAL(8) :: temp6d
    REAL(8) :: temp6dd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp7
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp7d
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp8
    REAL(8) :: temp9
    REAL(8) :: temp9d
    REAL(8) :: temp10
    REAL(8) :: temp10d
    REAL(8) :: temp11
    REAL(8) :: temp11d
    REAL(8) :: temp12
    REAL(8) :: temp12d
    REAL(8) :: temp13
    REAL(8) :: temp14
    REAL(8) :: temp14d
    REAL(8) :: temp15
    REAL(8) :: temp15d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp16
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp17
    REAL(8) :: temp18
    REAL(8) :: temp19
    REAL(8) :: temp20
    REAL(8) :: temp21
    REAL(8) :: temp22
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask0
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask1
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask2
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask3
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask4
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask5
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask6
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask7
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask8
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask9
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask10
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask11
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask12
    arg1d(:) = td/tc
    arg1d0(:) = td0/tc
    arg1d1(:) = td1/tc
    arg1(:) = t/tc
    temp16 = SQRT(arg1(:))
    mask(:) = arg1(:) .EQ. 0.0
    WHERE (mask(:)) 
      temp7d = 0.0_8
    ELSEWHERE
      temp7d = arg1d1(:)/(2.0*temp16)
    END WHERE
    temp7 = temp16
    mask0(:) = arg1(:) .EQ. 0.0
    WHERE (mask0(:)) tempd = 0.0_8
    tempdd = 0.0_8
    mask1(:) = .NOT.arg1(:) .EQ. 0.0
    WHERE (mask1(:)) 
      temp16 = arg1d0(:)/(2.0*temp7)
      tempdd = -(temp16*temp7d/temp7)
      tempd = temp16
    END WHERE
    tempd0 = temp7d
    temp = temp7
    mask2(:) = arg1(:) .EQ. 0.0
    WHERE (mask2(:)) result1d = 0.0_8
    result1dd = 0.0_8
    mask3(:) = .NOT.arg1(:) .EQ. 0.0
    WHERE (mask3(:)) 
      temp16 = arg1d(:)/(2.0*temp)
      temp7d = -(temp16*tempd0/temp)
      temp7 = temp16
    END WHERE
    result1ddd = 0.0_8
    mask4(:) = .NOT.arg1(:) .EQ. 0.0
    WHERE (mask4(:)) 
      temp16 = temp7*tempd/temp
      result1ddd = -((tempd*temp7d+temp7*tempdd-temp16*tempd0)/temp)
      result1dd = -temp16
    END WHERE
    result1dd0 = 0.0_8
    mask5(:) = .NOT.arg1(:) .EQ. 0.0
    WHERE (mask5(:)) 
      result1dd0 = temp7d
      result1d = temp7
    END WHERE
    result1d0d = tempdd
    result1d0 = tempd
    result1d1 = tempd0
    result1 = temp
    temp8 = 2*ac*k
    arg2ddd(:) = -(temp8*((k*(1.0-result1)+1.0)*result1ddd-result1dd*k*&
&     result1d1-k*(result1d0*result1dd0+result1d*result1d0d)))
    arg2dd(:) = -(temp8*((k*(1.0-result1)+1.0)*result1dd-result1d*k*&
&     result1d0))
    arg2dd0(:) = -(temp8*((k*(1.0-result1)+1.0)*result1dd0-result1d*k*&
&     result1d1))
    arg2d(:) = -(temp8*((k*(1.0-result1)+1.0)*result1d))
    temp17 = 2*ac*k
    arg2d0d(:) = -(temp17*((k*(1.0-result1)+1.0)*result1d0d-result1d0*k*&
&     result1d1))
    arg2d0(:) = -(temp17*((k*(1.0-result1)+1.0)*result1d0))
    arg2d1(:) = -(ac*2*(k*(1.0-result1)+1.0)*k*result1d1)
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    temp16 = SQRT(arg2(:))
    mask6(:) = arg2(:) .EQ. 0.0
    WHERE (mask6(:)) 
      temp7d = 0.0_8
    ELSEWHERE
      temp7d = arg2d1(:)/(2.0*temp16)
    END WHERE
    temp7 = temp16
    mask7(:) = arg2(:) .EQ. 0.0
    WHERE (mask7(:)) temp0d = 0.0_8
    temp0dd = 0.0_8
    mask8(:) = .NOT.arg2(:) .EQ. 0.0
    WHERE (mask8(:)) 
      temp16 = arg2d0(:)/(2.0*temp7)
      temp0dd = (arg2d0d(:)-temp16*2.0*temp7d)/(2.0*temp7)
      temp0d = temp16
    END WHERE
    temp0d0 = temp7d
    temp0 = temp7
    mask9(:) = arg2(:) .EQ. 0.0
    WHERE (mask9(:)) ad = 0.0_8
    add = 0.0_8
    mask10(:) = .NOT.arg2(:) .EQ. 0.0
    WHERE (mask10(:)) 
      temp16 = arg2d(:)/(2.0*temp0)
      temp7d = (arg2dd0(:)-temp16*2.0*temp0d0)/(2.0*temp0)
      temp7 = temp16
    END WHERE
    addd = 0.0_8
    mask11(:) = .NOT.arg2(:) .EQ. 0.0
    WHERE (mask11(:)) 
      temp16 = (arg2dd(:)-2.0*temp7*temp0d)/(2.0*temp0)
      addd = (arg2ddd(:)-2.0*(temp0d*temp7d+temp7*temp0dd)-temp16*2.0*&
&       temp0d0)/(2.0*temp0)
      add = temp16
    END WHERE
    add0 = 0.0_8
    mask12(:) = .NOT.arg2(:) .EQ. 0.0
    WHERE (mask12(:)) 
      add0 = temp7d
      ad = temp7
    END WHERE
    ad0d = temp0dd
    ad0 = temp0d
    ad1 = temp0d0
    a = temp0
    amix = 0.0
    bmix = 0.0
    bmixd = 0.0_8
    amixd = 0.0_8
    amixd0 = 0.0_8
    amixdd = 0.0_8
    amixddd = 0.0_8
    amixd1 = 0.0_8
    amixd0d = 0.0_8
    amixdd0 = 0.0_8
    DO i=1,SIZE(n)-1
      DO j=i+1,SIZE(n)
        nijd = n(j)*nd(i) + n(i)*nd(j)
        nij = n(i)*n(j)
        temp1 = 2*(-kij(i, j)+1)
        temp9d = nijd*ad1(i) + nij*add0(i)
        temp9 = nijd*a(i) + nij*ad(i)
        temp18 = nijd*ad0(i) + nij*add(i)
        amixddd = amixddd + temp1*(ad0(j)*temp9d+temp9*ad0d(j)+temp18*&
&         ad1(j)+a(j)*(nijd*ad0d(i)+nij*addd(i))+nij*(ad0(i)*add0(j)+ad(&
&         j)*ad0d(i)+add(j)*ad1(i)+a(i)*addd(j)))
        amixdd = amixdd + temp1*(temp9*ad0(j)+a(j)*temp18+nij*(ad(j)*ad0&
&         (i)+a(i)*add(j)))
        amixdd0 = amixdd0 + temp1*(temp9*ad1(j)+a(j)*temp9d+nij*(ad(j)*&
&         ad1(i)+a(i)*add0(j)))
        amixd = amixd + temp1*(a(j)*temp9+nij*(a(i)*ad(j)))
        amixd0d = amixd0d + temp1*nij*(ad0(i)*ad1(j)+a(j)*ad0d(i)+ad0(j)&
&         *ad1(i)+a(i)*ad0d(j))
        amixd0 = amixd0 + temp1*nij*(a(j)*ad0(i)+a(i)*ad0(j))
        amixd1 = amixd1 + temp1*nij*(a(j)*ad1(i)+a(i)*ad1(j))
        amix = amix + temp1*(nij*a(i)*a(j))
        temp1 = b(i) + b(j)
        bmixd = bmixd + temp1*(1-lij(i, j))*nijd
        bmix = bmix + temp1*((1-lij(i, j))*nij)
      END DO
    END DO
    arg10ddd(:) = 2**2*n*nd*(ad0*ad1+a*ad0d) + n**2*2*(ad0*add0+ad*ad0d+&
&     add*ad1+a*addd)
    arg10dd(:) = n*2**2*nd*a*ad0 + n**2*2*(ad*ad0+a*add)
    arg10dd0(:) = n*2**2*nd*a*ad1 + n**2*2*(ad*ad1+a*add0)
    arg10d(:) = a**2*2*n*nd + n**2*2*a*ad
    arg10d0d(:) = n**2*2*(ad0*ad1+a*ad0d)
    arg10d0(:) = n**2*2*a*ad0
    arg10d1(:) = n**2*2*a*ad1
    arg10(:) = n**2*a**2
    amixddd = amixddd + SUM(arg10ddd(:))
    amixdd = amixdd + SUM(arg10dd(:))
    amixdd0 = amixdd0 + SUM(arg10dd0(:))
    amixd = amixd + SUM(arg10d(:))
    amixd0d = amixd0d + SUM(arg10d0d(:))
    amixd0 = amixd0 + SUM(arg10d0(:))
    amixd1 = amixd1 + SUM(arg10d1(:))
    amix = amix + SUM(arg10(:))
    arg11d(:) = b*2*n*nd
    arg11(:) = n**2*b
    bmixd = bmixd + SUM(arg11d(:))
    bmix = bmix + SUM(arg11(:))
    temp1 = SUM(n)
    bmixd = (bmixd-bmix*SUM(nd)/temp1)/temp1
    bmix = bmix/temp1
    temp18 = bmix*vd/v
    temp9d = -(temp18*vd1/v)
    temp9 = temp18
    temp18 = (bmixd-temp9)/v
    temp10d = (-temp9d-temp18*vd1)/v
    temp10 = temp18
    temp18 = (temp9/v-temp10)/v
    b_vddd = vd0*((temp9d-temp9*vd1/v)/v-temp10d-temp18*vd1)/v
    b_vdd = vd0*temp18
    b_vdd0 = temp10d
    b_vd = temp10
    temp18 = bmix*vd0/(v*v)
    b_vd0d = temp18*2*vd1/v
    b_vd0 = -temp18
    b_vd1 = -(bmix*vd1/v**2)
    b_v = bmix/v
    temp18 = (del1*b_v+1.0)/(del2*b_v+1.0)
    temp10d = (del1-temp18*del2)*b_vd1/(del2*b_v+1.0)
    temp10 = temp18
    temp18 = b_vd0/(del2*b_v+1.0)
    temp1dd = (del1-del2*temp10)*(b_vd0d-temp18*del2*b_vd1)/(del2*b_v+&
&     1.0) - temp18*del2*temp10d
    temp1d = (del1-del2*temp10)*temp18
    temp1d0 = temp10d
    temp1 = temp10
    temp18 = b_vd/(del2*b_v+1.0)
    temp10d = (b_vdd0-temp18*del2*b_vd1)/(del2*b_v+1.0)
    temp10 = temp18
    temp18 = b_vdd - del2*temp10*b_vd0
    temp19 = (del1-del2*temp1)/(del2*b_v+1.0)
    arg12ddd = temp18*(-(del2*temp1d0)-temp19*del2*b_vd1)/(del2*b_v+1.0)&
&     + temp19*(b_vddd-del2*(b_vd0*temp10d+temp10*b_vd0d)) - del2*(&
&     temp1d*temp10d+temp10*temp1dd)
    arg12dd = temp19*temp18 - del2*(temp10*temp1d)
    arg12dd0 = (del1-del2*temp1)*temp10d - temp10*del2*temp1d0
    arg12d = (del1-del2*temp1)*temp10
    arg12d0d = temp1dd
    arg12d0 = temp1d
    arg12d1 = temp1d0
    arg12 = temp1
    temp19 = b_vd0/(-b_v+1.0)
    temp1dd = -((b_vd0d+temp19*b_vd1)/(1.0-b_v))
    temp1d = -temp19
    temp1d0 = -(b_vd1/(1.0-b_v))
    temp1 = LOG(-b_v + 1.0)
    temp2 = SUM(n)
    temp3d = r*bmix*(del1-del2)*td0
    temp3d0 = r*bmix*(del1-del2)*td1
    temp3 = r*(del1-del2)*t*bmix
    temp19 = (amixd0-temp3d*amix/temp3)/temp3
    temp4dd = (amixd0d-temp3d*(amixd1-amix*temp3d0/temp3)/temp3-temp19*&
&     temp3d0)/temp3
    temp4d = temp19
    temp4d0 = (amixd1-amix*temp3d0/temp3)/temp3
    temp4 = amix/temp3
    temp5dd = (arg12d0d-arg12d0*arg12d1/arg12)/arg12
    temp5d = arg12d0/arg12
    temp5d0 = arg12d1/arg12
    temp5 = LOG(arg12)
    temp6dd = -(temp2*temp1dd) - temp5d*temp4d0 - temp4*temp5dd - temp4d&
&     *temp5d0 - temp5*temp4dd
    temp6d = -(temp2*temp1d) - temp4*temp5d - temp5*temp4d
    temp6d0 = -(temp2*temp1d0) - temp4*temp5d0 - temp5*temp4d0
    temp6 = -(temp2*temp1) - temp5*temp4
    temp10d = (temp5d0-temp5*temp3d0/temp3)/temp3
    temp10 = temp5/temp3
    temp9d = bmixd*td1
    temp9 = bmix*td + bmixd*t
    temp11d = amixdd0 - r*(del1-del2)*(temp9*temp4d0+temp4*temp9d)
    temp11 = amixd - r*(del1-del2)*temp4*temp9
    temp19 = temp4*arg12d/arg12
    temp12d = (arg12d*temp4d0+temp4*arg12dd0-temp19*arg12d1)/arg12
    temp12 = temp19
    temp13 = SUM(nd)
    temp19 = b_vd/(-b_v+1.0)
    temp14d = (b_vdd0+temp19*b_vd1)/(1.0-b_v)
    temp14 = temp19
    temp15d = temp2*temp14d - temp13*temp1d0 - temp12d - temp10*temp11d &
&     - temp11*temp10d
    temp15 = temp2*temp14 - temp13*temp1 - temp12 - temp11*temp10
    temp19 = temp11/temp3
    temp18 = amixdd - r*(del1-del2)*(temp9*temp4d+bmixd*td0*temp4)
    temp20 = (arg12d*temp4d+temp4*arg12dd-temp12*arg12d0)/arg12
    temp21 = (b_vdd+temp14*b_vd0)/(-b_v+1.0)
    temp22 = temp2*temp21 - temp13*temp1d - temp20 - temp10*temp18 - (&
&     temp5d-temp3d*temp10)*temp19
    arvalddd = r*(td0*temp15d+temp22*td1+t*(temp2*(b_vddd+b_vd0*temp14d+&
&     temp14*b_vd0d+temp21*b_vd1)/(1.0-b_v)-temp13*temp1dd-(temp4d*&
&     arg12dd0+arg12d*temp4dd+arg12dd*temp4d0+temp4*arg12ddd-arg12d0*&
&     temp12d-temp12*arg12d0d-temp20*arg12d1)/arg12-temp18*temp10d-&
&     temp10*(amixddd-r*(del1-del2)*(temp4d*temp9d+temp9*temp4dd+bmixd*&
&     td0*temp4d0))-temp19*(temp5dd-temp3d*temp10d)-(temp5d-temp3d*&
&     temp10)*(temp11d-temp19*temp3d0)/temp3)+td*temp6dd)
    arvaldd = r*(td0*temp15+t*temp22+td*temp6d)
    arvaldd0 = r*(temp15*td1+t*temp15d+td*temp6d0)
    arvald = r*(t*temp15+td*temp6)
    arvald0d = r*(temp6d*td1+t*temp6dd+td0*temp6d0)
    arvald0 = r*(t*temp6d+temp6*td0)
    arvald1 = r*(t*temp6d0+temp6*td1)
    arval = r*(temp6*t)
  END SUBROUTINE AR_D_D_D

!  Differentiation of ar_d in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: arval arvald
!   with respect to varying inputs: t v
!   RW status of diff variables: t:in v:in arval:out arvald:out
!  Differentiation of ar in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: arval
!   with respect to varying inputs: n t v
!   RW status of diff variables: n:in t:in v:in arval:out
  SUBROUTINE AR_D_D(n, nd, v, vd0, vd, t, td0, td, arval, arvald0, &
&   arvald, arvaldd)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8), INTENT(IN) :: vd0, td0
    REAL(8), INTENT(IN) :: nd(:), vd, td
    REAL(8), INTENT(OUT) :: arval
    REAL(8), INTENT(OUT) :: arvald0
    REAL(8), INTENT(OUT) :: arvald
    REAL(8), INTENT(OUT) :: arvaldd
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: amixd0, ad0(SIZE(n))
    REAL(8) :: amixd, ad(SIZE(n)), nijd
    REAL(8) :: amixdd, add(SIZE(n))
    REAL(8) :: bmix
    REAL(8) :: bmixd
    REAL(8) :: b_v
    REAL(8) :: b_vd0
    REAL(8) :: b_vd
    REAL(8) :: b_vdd
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1dd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2dd
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n)) :: arg10d0
    REAL(8), DIMENSION(SIZE(n)) :: arg10d
    REAL(8), DIMENSION(SIZE(n)) :: arg10dd
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11d
    REAL(8) :: arg12
    REAL(8) :: arg12d0
    REAL(8) :: arg12d
    REAL(8) :: arg12dd
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempd
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0d
    REAL(8) :: temp1
    REAL(8) :: temp1d
    REAL(8) :: temp2
    REAL(8) :: temp3
    REAL(8) :: temp3d
    REAL(8) :: temp4
    REAL(8) :: temp4d
    REAL(8) :: temp5
    REAL(8) :: temp5d
    REAL(8) :: temp6
    REAL(8) :: temp6d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp7
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp8
    REAL(8) :: temp9
    REAL(8) :: temp10
    REAL(8) :: temp11
    REAL(8) :: temp12
    REAL(8) :: temp13
    REAL(8) :: temp14
    REAL(8) :: temp15
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask0
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask1
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask2
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask3
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask4
    arg1d(:) = td/tc
    arg1d0(:) = td0/tc
    arg1(:) = t/tc
    temp7 = SQRT(arg1(:))
    mask(:) = arg1(:) .EQ. 0.0
    WHERE (mask(:)) 
      tempd = 0.0_8
    ELSEWHERE
      tempd = arg1d0(:)/(2.0*temp7)
    END WHERE
    temp = temp7
    mask0(:) = arg1(:) .EQ. 0.0
    WHERE (mask0(:)) result1d = 0.0_8
    result1dd = 0.0_8
    mask1(:) = .NOT.arg1(:) .EQ. 0.0
    WHERE (mask1(:)) 
      temp7 = arg1d(:)/(2.0*temp)
      result1dd = -(temp7*tempd/temp)
      result1d = temp7
    END WHERE
    result1d0 = tempd
    result1 = temp
    temp8 = 2*ac*k
    arg2dd(:) = -(temp8*((k*(1.0-result1)+1.0)*result1dd-result1d*k*&
&     result1d0))
    arg2d(:) = -(temp8*((k*(1.0-result1)+1.0)*result1d))
    arg2d0(:) = -(ac*2*(k*(1.0-result1)+1.0)*k*result1d0)
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    temp7 = SQRT(arg2(:))
    mask2(:) = arg2(:) .EQ. 0.0
    WHERE (mask2(:)) 
      temp0d = 0.0_8
    ELSEWHERE
      temp0d = arg2d0(:)/(2.0*temp7)
    END WHERE
    temp0 = temp7
    mask3(:) = arg2(:) .EQ. 0.0
    WHERE (mask3(:)) ad = 0.0_8
    add = 0.0_8
    mask4(:) = .NOT.arg2(:) .EQ. 0.0
    WHERE (mask4(:)) 
      temp7 = arg2d(:)/(2.0*temp0)
      add = (arg2dd(:)-temp7*2.0*temp0d)/(2.0*temp0)
      ad = temp7
    END WHERE
    ad0 = temp0d
    a = temp0
    amix = 0.0
    bmix = 0.0
    bmixd = 0.0_8
    amixd = 0.0_8
    amixd0 = 0.0_8
    amixdd = 0.0_8
    DO i=1,SIZE(n)-1
      DO j=i+1,SIZE(n)
        nijd = n(j)*nd(i) + n(i)*nd(j)
        nij = n(i)*n(j)
        temp1 = 2*(-kij(i, j)+1)
        temp9 = nijd*a(i) + nij*ad(i)
        amixdd = amixdd + temp1*(temp9*ad0(j)+a(j)*(nijd*ad0(i)+nij*add(&
&         i))+nij*(ad(j)*ad0(i)+a(i)*add(j)))
        amixd = amixd + temp1*(a(j)*temp9+nij*(a(i)*ad(j)))
        amixd0 = amixd0 + temp1*nij*(a(j)*ad0(i)+a(i)*ad0(j))
        amix = amix + temp1*(nij*a(i)*a(j))
        temp1 = b(i) + b(j)
        bmixd = bmixd + temp1*(1-lij(i, j))*nijd
        bmix = bmix + temp1*((1-lij(i, j))*nij)
      END DO
    END DO
    arg10dd(:) = n*2**2*nd*a*ad0 + n**2*2*(ad*ad0+a*add)
    arg10d(:) = a**2*2*n*nd + n**2*2*a*ad
    arg10d0(:) = n**2*2*a*ad0
    arg10(:) = n**2*a**2
    amixdd = amixdd + SUM(arg10dd(:))
    amixd = amixd + SUM(arg10d(:))
    amixd0 = amixd0 + SUM(arg10d0(:))
    amix = amix + SUM(arg10(:))
    arg11d(:) = b*2*n*nd
    arg11(:) = n**2*b
    bmixd = bmixd + SUM(arg11d(:))
    bmix = bmix + SUM(arg11(:))
    temp1 = SUM(n)
    bmixd = (bmixd-bmix*SUM(nd)/temp1)/temp1
    bmix = bmix/temp1
    temp9 = bmix*vd/v
    temp10 = (bmixd-temp9)/v
    b_vdd = (temp9/v-temp10)*vd0/v
    b_vd = temp10
    b_vd0 = -(bmix*vd0/v**2)
    b_v = bmix/v
    temp10 = (del1*b_v+1.0)/(del2*b_v+1.0)
    temp1d = (del1-temp10*del2)*b_vd0/(del2*b_v+1.0)
    temp1 = temp10
    temp10 = b_vd/(del2*b_v+1.0)
    arg12dd = (del1-del2*temp1)*(b_vdd-temp10*del2*b_vd0)/(del2*b_v+1.0)&
&     - temp10*del2*temp1d
    arg12d = (del1-del2*temp1)*temp10
    arg12d0 = temp1d
    arg12 = temp1
    temp1d = -(b_vd0/(1.0-b_v))
    temp1 = LOG(-b_v + 1.0)
    temp2 = SUM(n)
    temp3d = r*bmix*(del1-del2)*td0
    temp3 = r*(del1-del2)*t*bmix
    temp4d = (amixd0-amix*temp3d/temp3)/temp3
    temp4 = amix/temp3
    temp5d = arg12d0/arg12
    temp5 = LOG(arg12)
    temp6d = -(temp2*temp1d) - temp4*temp5d - temp5*temp4d
    temp6 = -(temp2*temp1) - temp5*temp4
    temp10 = temp5/temp3
    temp9 = bmix*td + bmixd*t
    temp11 = amixd - r*(del1-del2)*temp4*temp9
    temp12 = temp4*arg12d/arg12
    temp13 = SUM(nd)
    temp14 = b_vd/(-b_v+1.0)
    temp15 = temp2*temp14 - temp13*temp1 - temp12 - temp11*temp10
    arvaldd = r*(temp15*td0+t*(temp2*(b_vdd+temp14*b_vd0)/(1.0-b_v)-&
&     temp13*temp1d-(arg12d*temp4d+temp4*arg12dd-temp12*arg12d0)/arg12-&
&     temp10*(amixdd-r*(del1-del2)*(temp9*temp4d+temp4*bmixd*td0))-&
&     temp11*(temp5d-temp10*temp3d)/temp3)+td*temp6d)
    arvald = r*(t*temp15+td*temp6)
    arvald0 = r*(t*temp6d+temp6*td0)
    arval = r*(temp6*t)
  END SUBROUTINE AR_D_D

!  Differentiation of ar_d in reverse (adjoint) mode (with options noISIZE):
!   gradient     of useful results: nd n t v arval arvald vd td
!   with respect to varying inputs: nd n t v arval arvald vd td
!   RW status of diff variables: nd:incr n:incr t:incr v:incr arval:in-zero
!                arvald:in-zero vd:incr td:incr
!  Differentiation of ar in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: arval
!   with respect to varying inputs: n t v
!   RW status of diff variables: n:in t:in v:in arval:out
  SUBROUTINE AR_D_B(n, nb, nd, ndb, v, vb, vd, vdb, t, tb, td, tdb, &
&   arval, arvalb, arvald, arvaldb)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8) :: nb(:), vb, tb
    REAL(8), INTENT(IN) :: nd(:), vd, td
    REAL(8) :: ndb(:), vdb, tdb
    REAL(8) :: arval
    REAL(8) :: arvalb
    REAL(8) :: arvald
    REAL(8) :: arvaldb
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: amixb, ab(SIZE(n)), nijb
    REAL(8) :: amixd, ad(SIZE(n)), nijd
    REAL(8) :: amixdb, adb(SIZE(n)), nijdb
    REAL(8) :: bmix
    REAL(8) :: bmixb
    REAL(8) :: bmixd
    REAL(8) :: bmixdb
    REAL(8) :: b_v
    REAL(8) :: b_vb
    REAL(8) :: b_vd
    REAL(8) :: b_vdb
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1b
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1db
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1b
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1db
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2b
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2db
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n)) :: arg10b
    REAL(8), DIMENSION(SIZE(n)) :: arg10d
    REAL(8), DIMENSION(SIZE(n)) :: arg10db
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11b
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11d
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11db
    REAL(8) :: arg12
    REAL(8) :: arg12b
    REAL(8) :: arg12d
    REAL(8) :: arg12db
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempb
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0b
    REAL(8) :: temp1
    REAL(8) :: temp1b
    REAL(8) :: temp2
    REAL(8) :: temp2b
    REAL(8) :: temp3
    REAL(8) :: temp3b
    REAL(8) :: temp4
    REAL(8) :: temp4b
    REAL(8) :: temp5
    REAL(8) :: temp5b
    REAL(8) :: temp6
    REAL(8) :: temp6b
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask0
    REAL(8), DIMENSION(SIZE(tc, 1)) :: tempb0
    REAL(8) :: temp7
    REAL(8) :: tempb1
    REAL(8) :: temp8
    REAL(8) :: tempb2
    REAL(8), DIMENSION(SIZE(n, 1)) :: tempb3
    REAL(8), DIMENSION(SIZE(n)) :: tempb4
    REAL(8) :: temp9
    REAL(8) :: tempb5
    REAL(8) :: temp10
    REAL(8) :: temp11
    REAL(8) :: temp12
    REAL(8) :: temp13
    REAL(8) :: tempb6
    REAL(8) :: tempb7
    REAL(8) :: tempb8
    REAL(8) :: tempb9
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_to0
    EXTERNAL PUSHINTEGER4
    EXTERNAL PUSHREAL8
    EXTERNAL POPREAL8
    EXTERNAL POPINTEGER4
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask1
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask2
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask3
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask4
    arg1d(:) = td/tc
    arg1(:) = t/tc
    temp = SQRT(arg1(:))
    mask(:) = arg1(:) .EQ. 0.0
    WHERE (mask(:)) 
      result1d = 0.0_8
    ELSEWHERE
      result1d = arg1d(:)/(2.0*temp)
    END WHERE
    result1 = temp
    arg2d(:) = -(ac*2*(k*(1.0-result1)+1.0)*k*result1d)
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    temp0 = SQRT(arg2(:))
    mask0(:) = arg2(:) .EQ. 0.0
    WHERE (mask0(:)) 
      ad = 0.0_8
    ELSEWHERE
      ad = arg2d(:)/(2.0*temp0)
    END WHERE
    a = temp0
    amix = 0.0
    bmix = 0.0
    bmixd = 0.0_8
    amixd = 0.0_8
    DO i=1,SIZE(n)-1
      ad_from = i + 1
      DO j=ad_from,SIZE(n)
        nijd = n(j)*nd(i) + n(i)*nd(j)
        nij = n(i)*n(j)
        temp1 = 2*(-kij(i, j)+1)
        amixd = amixd + temp1*(a(j)*(a(i)*nijd+nij*ad(i))+nij*a(i)*ad(j)&
&         )
        amix = amix + temp1*(nij*a(i)*a(j))
        temp1 = b(i) + b(j)
        bmixd = bmixd + temp1*(1-lij(i, j))*nijd
        bmix = bmix + temp1*((1-lij(i, j))*nij)
      END DO
      CALL PUSHINTEGER4(j - 1)
      CALL PUSHINTEGER4(ad_from)
    END DO
    CALL PUSHINTEGER4(i - 1)
    arg10d(:) = a**2*2*n*nd + n**2*2*a*ad
    arg10(:) = n**2*a**2
    amixd = amixd + SUM(arg10d(:))
    amix = amix + SUM(arg10(:))
    arg11d(:) = b*2*n*nd
    arg11(:) = n**2*b
    bmixd = bmixd + SUM(arg11d(:))
    bmix = bmix + SUM(arg11(:))
    temp1 = SUM(n)
    CALL PUSHREAL8(bmixd)
    bmixd = (bmixd-bmix*SUM(nd)/temp1)/temp1
    CALL PUSHREAL8(bmix)
    bmix = bmix/temp1
    b_vd = (bmixd-bmix*vd/v)/v
    b_v = bmix/v
    CALL PUSHREAL8(temp1)
    temp1 = (del1*b_v+1.0)/(del2*b_v+1.0)
    arg12d = (del1-temp1*del2)*b_vd/(del2*b_v+1.0)
    arg12 = temp1
    temp1 = LOG(-b_v + 1.0)
    temp2 = SUM(n)
    temp3 = r*(del1-del2)*t*bmix
    temp4 = amix/temp3
    temp5 = LOG(arg12)
    temp6 = -(temp2*temp1) - temp5*temp4
    tempb5 = r*arvaldb
    temp13 = temp2*b_vd/(-b_v+1.0)
    temp12 = SUM(nd)
    temp11 = temp4*arg12d/arg12
    temp7 = bmix*td + t*bmixd
    temp10 = amixd - r*(del1-del2)*temp4*temp7
    temp8 = temp5/temp3
    temp6b = t*r*arvalb + td*tempb5
    tb = tb + temp6*r*arvalb + (temp13-temp1*temp12-temp11-temp10*temp8)&
&     *tempb5
    tempb6 = t*tempb5
    tempb7 = tempb6/(1.0-b_v)
    temp1b = -(temp12*tempb6) - temp2*temp6b
    ndb = ndb - temp1*tempb6
    tempb8 = -(tempb6/arg12)
    amixdb = -(temp8*tempb6)
    tempb9 = r*(del1-del2)*temp8*tempb6
    tempb2 = -(temp10*tempb6/temp3)
    temp5b = tempb2 - temp4*temp6b
    temp4b = temp7*tempb9 + arg12d*tempb8 - temp5*temp6b
    temp3b = -(temp8*tempb2) - amix*temp4b/temp3**2
    tempb1 = temp4*tempb9
    tdb = tdb + temp6*tempb5 + bmix*tempb1
    arg12db = temp4*tempb8
    arg12b = temp5b/arg12 - temp11*tempb8
    temp2b = b_vd*tempb7 - temp1*temp6b
    amixb = temp4b/temp3
    tempb5 = r*(del1-del2)*temp3b
    bmixb = td*tempb1 + t*tempb5
    tb = tb + bmixd*tempb1 + bmix*tempb5
    temp1 = (del1*b_v+1.0)/(del2*b_v+1.0)
    temp9 = b_vd/(del2*b_v+1.0)
    tempb5 = (del1-del2*temp1)*arg12db/(del2*b_v+1.0)
    b_vdb = temp2*tempb7 + tempb5
    b_vb = temp13*tempb7 - temp1b/(1.0-b_v) - del2*temp9*tempb5
    temp1b = arg12b - del2*temp9*arg12db
    CALL POPREAL8(temp1)
    tempb5 = temp1b/(del2*b_v+1.0)
    b_vb = b_vb + (del1-del2*(del1*b_v+1.0)/(del2*b_v+1.0))*tempb5
    temp8 = bmix*vd/v
    tempb5 = b_vdb/v
    bmixdb = t*tempb1 + tempb5
    tempb2 = -(tempb5/v)
    bmixb = bmixb + b_vb/v + vd*tempb2
    vb = vb - bmix*b_vb/v**2 - (bmixd-temp8)*tempb5/v - temp8*tempb2
    vdb = vdb + bmix*tempb2
    CALL POPREAL8(bmix)
    CALL POPREAL8(bmixd)
    temp7 = bmix/temp1
    temp9 = SUM(nd)
    tempb2 = bmixdb/temp1
    bmixdb = tempb2
    tempb1 = -(temp9*tempb2/temp1)
    temp1b = -(bmix*bmixb/temp1**2) - (bmixd-temp9*temp7)*tempb2/temp1 -&
&     temp7*tempb1
    bmixb = bmixb/temp1 + tempb1
    arg11b = 0.0_8
    arg11b = bmixb
    arg11db = 0.0_8
    arg11db = bmixdb
    arg10b = 0.0_8
    arg10b = amixb
    arg10db = 0.0_8
    arg10db = amixdb
    ab = 0.0_8
    adb = 0.0_8
    tempb3 = 2*a**2*arg10db
    nb = nb + temp2b + temp1b + 2*n*b*arg11b + nd*b*2*arg11db + 2*n*a**2&
&     *arg10b + 2**2*n*a*ad*arg10db + nd*tempb3
    ndb = ndb + n*b*2*arg11db - temp7*tempb2 + n*tempb3
    tempb4 = 2*n**2*arg10db
    ab = 2*a*n**2*arg10b + 2**2*a*n*nd*arg10db + ad*tempb4
    adb = a*tempb4
    CALL POPINTEGER4(ad_to0)
    DO i=ad_to0,1,-1
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO j=ad_to,ad_from,-1
        nij = n(i)*n(j)
        temp1 = b(i) + b(j)
        nijb = temp1*(1-lij(i, j))*bmixb
        nijd = n(j)*nd(i) + n(i)*nd(j)
        nijdb = temp1*(1-lij(i, j))*bmixdb
        temp1 = 2*(-kij(i, j)+1)
        tempb2 = a(j)*temp1*amixb
        ab(j) = ab(j) + nij*a(i)*temp1*amixb
        nijb = nijb + a(i)*tempb2
        ab(i) = ab(i) + nij*tempb2
        tempb1 = temp1*amixdb
        ab(j) = ab(j) + (a(i)*nijd+nij*ad(i))*tempb1
        tempb2 = a(j)*tempb1
        nijb = nijb + a(i)*ad(j)*tempb1 + ad(i)*tempb2
        ab(i) = ab(i) + nij*ad(j)*tempb1 + nijd*tempb2
        adb(j) = adb(j) + nij*a(i)*tempb1
        nijdb = nijdb + a(i)*tempb2
        adb(i) = adb(i) + nij*tempb2
        nb(i) = nb(i) + n(j)*nijb
        nb(j) = nb(j) + n(i)*nijb + nd(i)*nijdb
        ndb(i) = ndb(i) + n(j)*nijdb
        nb(i) = nb(i) + nd(j)*nijdb
        ndb(j) = ndb(j) + n(i)*nijdb
      END DO
    END DO
    temp0b = 0.0_8
    temp0b = ab
    arg2db = 0.0_8
    arg2b = 0.0_8
    result1b = 0.0_8
    result1db = 0.0_8
    tempb = 0.0_8
    arg1db = 0.0_8
    arg1b = 0.0_8
    mask1(:) = .NOT.mask0(:)
    WHERE (mask1(:)) 
      tempb0 = adb/(2.0*temp0)
      arg2db = tempb0
      temp0b = temp0b - arg2d*tempb0/temp0
    END WHERE
    tempb0 = -(ac*2*k*arg2db)
    mask2 = arg2 .EQ. 0.0
    WHERE (mask2) 
      arg2b = 0.0_8
    ELSEWHERE
      arg2b = temp0b/(2.0*SQRT(arg2))
    END WHERE
    result1b = -(k*2*(k*(1.0-result1)+1.0)*ac*arg2b) - k*result1d*tempb0
    result1db = (k*(1.0-result1)+1.0)*tempb0
    tempb = result1b
    mask3(:) = .NOT.mask(:)
    WHERE (mask3(:)) 
      tempb0 = result1db/(2.0*temp)
      arg1db = tempb0
      tempb = tempb - arg1d*tempb0/temp
    END WHERE
    mask4 = arg1 .EQ. 0.0
    WHERE (mask4) 
      arg1b = 0.0_8
    ELSEWHERE
      arg1b = tempb/(2.0*SQRT(arg1))
    END WHERE
    tb = tb + SUM(arg1b/tc)
    tdb = tdb + SUM(arg1db/tc)
    arvalb = 0.0_8
    arvaldb = 0.0_8
  END SUBROUTINE AR_D_B

!  Differentiation of ar in forward (tangent) mode (with options noISIZE):
!   variations   of useful results: arval
!   with respect to varying inputs: n t v
!   RW status of diff variables: n:in t:in v:in arval:out
  SUBROUTINE AR_D(n, nd, v, vd, t, td, arval, arvald)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8), INTENT(IN) :: nd(:), vd, td
    REAL(8), INTENT(OUT) :: arval
    REAL(8), INTENT(OUT) :: arvald
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: amixd, ad(SIZE(n)), nijd
    REAL(8) :: bmix
    REAL(8) :: bmixd
    REAL(8) :: b_v
    REAL(8) :: b_vd
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2d
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n)) :: arg10d
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11d
    REAL(8) :: arg12
    REAL(8) :: arg12d
    REAL(8), DIMENSION(SIZE(tc, 1)) :: temp
    REAL(8), DIMENSION(SIZE(ac, 1)) :: temp0
    REAL(8) :: temp1
    REAL(8) :: temp2
    REAL(8) :: temp3
    REAL(8) :: temp4
    REAL(8) :: temp5
    REAL(8) :: temp6
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask
    LOGICAL, DIMENSION(SIZE(tc, 1)) :: mask0
    arg1d(:) = td/tc
    arg1(:) = t/tc
    temp = SQRT(arg1(:))
    mask(:) = arg1(:) .EQ. 0.0
    WHERE (mask(:)) 
      result1d = 0.0_8
    ELSEWHERE
      result1d = arg1d(:)/(2.0*temp)
    END WHERE
    result1 = temp
    arg2d(:) = -(ac*2*(k*(1.0-result1)+1.0)*k*result1d)
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    temp0 = SQRT(arg2(:))
    mask0(:) = arg2(:) .EQ. 0.0
    WHERE (mask0(:)) 
      ad = 0.0_8
    ELSEWHERE
      ad = arg2d(:)/(2.0*temp0)
    END WHERE
    a = temp0
    amix = 0.0
    bmix = 0.0
    bmixd = 0.0_8
    amixd = 0.0_8
    DO i=1,SIZE(n)-1
      DO j=i+1,SIZE(n)
        nijd = n(j)*nd(i) + n(i)*nd(j)
        nij = n(i)*n(j)
        temp1 = 2*(-kij(i, j)+1)
        amixd = amixd + temp1*(a(j)*(a(i)*nijd+nij*ad(i))+nij*a(i)*ad(j)&
&         )
        amix = amix + temp1*(nij*a(i)*a(j))
        temp1 = b(i) + b(j)
        bmixd = bmixd + temp1*(1-lij(i, j))*nijd
        bmix = bmix + temp1*((1-lij(i, j))*nij)
      END DO
    END DO
    arg10d(:) = a**2*2*n*nd + n**2*2*a*ad
    arg10(:) = n**2*a**2
    amixd = amixd + SUM(arg10d(:))
    amix = amix + SUM(arg10(:))
    arg11d(:) = b*2*n*nd
    arg11(:) = n**2*b
    bmixd = bmixd + SUM(arg11d(:))
    bmix = bmix + SUM(arg11(:))
    temp1 = SUM(n)
    bmixd = (bmixd-bmix*SUM(nd)/temp1)/temp1
    bmix = bmix/temp1
    b_vd = (bmixd-bmix*vd/v)/v
    b_v = bmix/v
    temp1 = (del1*b_v+1.0)/(del2*b_v+1.0)
    arg12d = (del1-temp1*del2)*b_vd/(del2*b_v+1.0)
    arg12 = temp1
    temp1 = LOG(-b_v + 1.0)
    temp2 = SUM(n)
    temp3 = r*(del1-del2)*t*bmix
    temp4 = amix/temp3
    temp5 = LOG(arg12)
    temp6 = -(temp2*temp1) - temp5*temp4
    arvald = r*(t*(temp2*b_vd/(1.0-b_v)-temp1*SUM(nd)-temp4*arg12d/arg12&
&     -temp5*(amixd-temp4*r*(del1-del2)*(bmix*td+t*bmixd))/temp3)+temp6*&
&     td)
    arval = r*(temp6*t)
  END SUBROUTINE AR_D

!  Differentiation of ar in reverse (adjoint) mode (with options noISIZE):
!   gradient     of useful results: arval
!   with respect to varying inputs: n t v arval
!   RW status of diff variables: n:out t:out v:out arval:in-zero
  SUBROUTINE AR_B(n, nb, v, vb, t, tb, arval, arvalb)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8) :: nb(:), vb, tb
    REAL(8) :: arval
    REAL(8) :: arvalb
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: amixb, ab(SIZE(n)), nijb
    REAL(8) :: bmix
    REAL(8) :: bmixb
    REAL(8) :: b_v
    REAL(8) :: b_vb
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1b
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1b
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2b
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n)) :: arg10b
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11b
    REAL(8) :: arg12
    REAL(8) :: arg12b
    REAL(8) :: temp
    REAL(8) :: tempb
    REAL(8) :: temp0
    REAL(8) :: temp1
    REAL(8) :: temp2
    REAL(8) :: temp3
    REAL(8) :: temp4
    REAL(8) :: tempb0
    REAL(8) :: tempb1
    INTEGER :: ad_from
    INTEGER :: ad_to
    INTEGER :: ad_to0
    arg1(:) = t/tc
    result1 = SQRT(arg1(:))
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    a = SQRT(arg2(:))
    amix = 0.0
    bmix = 0.0
    DO i=1,SIZE(n)-1
      ad_from = i + 1
      DO j=ad_from,SIZE(n)
        nij = n(i)*n(j)
        amix = amix + 2*nij*(a(i)*a(j))*(1-kij(i, j))
        bmix = bmix + nij*(b(i)+b(j))*(1-lij(i, j))
      END DO
      CALL PUSHINTEGER4(j - 1)
      CALL PUSHINTEGER4(ad_from)
    END DO
    CALL PUSHINTEGER4(i - 1)
    arg10(:) = n**2*a**2
    amix = amix + SUM(arg10(:))
    arg11(:) = n**2*b
    bmix = bmix + SUM(arg11(:))
    CALL PUSHREAL8(bmix)
    bmix = bmix/SUM(n)
    b_v = bmix/v
    arg12 = (1.0+del1*b_v)/(1.0+del2*b_v)
    nb = 0.0_8
    temp0 = LOG(-b_v + 1.0)
    temp1 = SUM(n)
    temp2 = r*(del1-del2)*t*bmix
    temp3 = amix/temp2
    temp4 = LOG(arg12)
    tempb = t*r*arvalb
    nb = -(temp0*tempb)
    b_vb = temp1*tempb/(1.0-b_v)
    arg12b = -(temp3*tempb/arg12)
    tempb0 = -(temp4*tempb/temp2)
    amixb = tempb0
    tempb1 = -(r*(del1-del2)*temp3*tempb0)
    tb = (-(temp1*temp0)-temp4*temp3)*r*arvalb + bmix*tempb1
    tempb = arg12b/(del2*b_v+1.0)
    b_vb = b_vb + (del1-del2*(del1*b_v+1.0)/(del2*b_v+1.0))*tempb
    bmixb = t*tempb1 + b_vb/v
    vb = -(bmix*b_vb/v**2)
    CALL POPREAL8(bmix)
    temp = SUM(n)
    nb = nb - bmix*bmixb/temp**2
    bmixb = bmixb/temp
    arg11b = 0.0_8
    arg11b = bmixb
    arg10b = 0.0_8
    arg10b = amixb
    nb = nb + 2*n*b*arg11b + 2*n*a**2*arg10b
    ab = 0.0_8
    ab = 2*a*n**2*arg10b
    CALL POPINTEGER4(ad_to0)
    DO i=ad_to0,1,-1
      CALL POPINTEGER4(ad_from)
      CALL POPINTEGER4(ad_to)
      DO j=ad_to,ad_from,-1
        tempb = (1-kij(i, j))*2*amixb
        nij = n(i)*n(j)
        nijb = (1-lij(i, j))*(b(i)+b(j))*bmixb + a(i)*a(j)*tempb
        ab(i) = ab(i) + nij*a(j)*tempb
        ab(j) = ab(j) + nij*a(i)*tempb
        nb(i) = nb(i) + n(j)*nijb
        nb(j) = nb(j) + n(i)*nijb
      END DO
    END DO
    arg2b = 0.0_8
    WHERE (arg2 .EQ. 0.0) 
      arg2b = 0.0_8
    ELSEWHERE
      arg2b = ab/(2.0*SQRT(arg2))
    END WHERE
    result1b = 0.0_8
    result1b = -(k*2*(k*(1.0-result1)+1.0)*ac*arg2b)
    arg1b = 0.0_8
    WHERE (arg1 .EQ. 0.0) 
      arg1b = 0.0_8
    ELSEWHERE
      arg1b = result1b/(2.0*SQRT(arg1))
    END WHERE
    tb = tb + SUM(arg1b/tc)
    arvalb = 0.0_8
  END SUBROUTINE AR_B

  SUBROUTINE AR(n, v, t, arval)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:), v, t
    REAL(8), INTENT(OUT) :: arval
    REAL(8) :: amix, a(SIZE(n)), ai(SIZE(n)), z2(SIZE(n)), nij
    REAL(8) :: bmix
    REAL(8) :: b_v
    REAL(8) :: aij(SIZE(n), SIZE(n)), bij(SIZE(n), SIZE(n))
    INTEGER :: i, j
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC LOG
    INTRINSIC SIZE
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: result1
    REAL(8), DIMENSION(SIZE(tc, 1)) :: arg2
    REAL(8), DIMENSION(SIZE(n)) :: arg10
    REAL(8), DIMENSION(SIZE(n, 1)) :: arg11
    REAL(8) :: arg12
    arg1(:) = t/tc
    result1 = SQRT(arg1(:))
    arg2(:) = ac*(1.0+k*(1.0-result1))**2
    a = SQRT(arg2(:))
    amix = 0.0
    bmix = 0.0
    DO i=1,SIZE(n)-1
      DO j=i+1,SIZE(n)
        nij = n(i)*n(j)
        amix = amix + 2*nij*(a(i)*a(j))*(1-kij(i, j))
        bmix = bmix + nij*(b(i)+b(j))*(1-lij(i, j))
      END DO
    END DO
    arg10(:) = n**2*a**2
    amix = amix + SUM(arg10(:))
    arg11(:) = n**2*b
    bmix = bmix + SUM(arg11(:))
    bmix = bmix/SUM(n)
    b_v = bmix/v
    arg12 = (1.0+del1*b_v)/(1.0+del2*b_v)
    arval = (-(SUM(n)*LOG(1.0-b_v))-amix/(r*t*bmix)*1.0/(del1-del2)*LOG(&
&     arg12))*(r*t)
  END SUBROUTINE AR

  PURE FUNCTION VOLUME_INITALIZER(n, p, t) RESULT (v0)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: n(:)
    REAL(8), INTENT(IN) :: p
    REAL(8), INTENT(IN) :: t
    REAL(8) :: v0
    INTRINSIC SUM
    v0 = SUM(n*b)/SUM(b)
  END FUNCTION VOLUME_INITALIZER

END MODULE 