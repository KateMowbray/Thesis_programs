MODULE M_rkckR 

 USE M_derivsR, ONLY: DERIVS
 USE GLOBAL 
 IMPLICIT NONE

  PRIVATE
  PUBLIC :: RKCK

 CONTAINS

  SUBROUTINE RKCK (R,DRDT,U,DUDT, T,H,MU,ROUT,UOUT,RERR,T1,T2)
 !####################################################################  
 !This part basically takes one Runga-Kutta step, given 6 position 
 !variables, t, and stepsize h. It calculates a fifth and fourth order
 !solution and takes the difference as the rerr (the error)
!#################################################################

  IMPLICIT NONE
  INTEGER :: I
  REAL(num) :: H, T,MU
  REAL(num), INTENT(IN)  		:: U, DUDT, T1, T2 
  REAL(num), INTENT(OUT)		:: UOUT
  REAL(num), DIMENSION(3), INTENT(IN) 	:: R, DRDT 
  REAL(num), DIMENSION(3), INTENT(OUT) 	:: ROUT
  !REAL(num), DIMENSION(3) 		:: UGB, UC
  REAL(num), DIMENSION(4) 		:: RERR
  REAL(num), DIMENSION(3) 		:: AK2, AK3, AK4, AK5, AK6, RTEMP
  REAL(num) 				:: UK2, UK3, UK4, UK5, UK6, UTEMP
  !REAL(num) 				:: GK2, GK3, GK4, GK5, GK6, GAMMATEMP
  !REAL(num) 				:: A2, A3, A4, A5, A6, B21, B31, B32, B41, B42, B43
  !REAL(num) 				:: B51, B52, B53, B54, B61, B62, B63, B64, B65
  !REAL(num) 				:: C1, C3, C4, C6, DC1, DC3, DC4, DC5, DC6
  REAL(num), PARAMETER 			:: A2=0.2_num, A3=0.3_num, A4=0.6_num, A5=1.0_num, A6=0.875_num
  REAL(num), PARAMETER 			:: B21=0.2_num,B31=0.075_num,B32=0.225_num,B41=0.3_num,B42=-0.9_num,B43=1.2_num
  REAL(num), PARAMETER 			:: B51=-11.0_num/54.0_num,B52=2.5_num, B53=-70.0_num/27.0_num,B54=35.0_num/27.0_num
  REAL(num), PARAMETER 			:: B61=1631.0_num/55296.0_num, B62=175.0_num/512.0_num, B63=575.0_num/13824.0_num 
  REAL(num), PARAMETER			:: B64=44275.0_num/110592.0_num, B65=253.0_num/4096.0_num
  REAL(num), PARAMETER 			:: C1=37.0_num/378.0_num, C3=250.0_num/621.0_num, C4=125.0_num/594.0_num
  REAL(num), PARAMETER			:: C6=512.0_num/1771.0_num
  REAL(num), PARAMETER 			:: DC1=C1-2825.0_num/27648.0_num,DC3=C3-18575.0_num/48384.0_num 
  REAL(num), PARAMETER			:: DC4=C4-13525.0_num/55296.0_num,DC5=-277.0_num/14336.0_num,DC6=C6-0.25_num
  
  !REAL(num), DIMENSION(3)		:: B,E,Vf
  !REAL(num), DIMENSION(3)		:: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT

  
  1101 format (A, 'RTEMP1=[',ES9.2,',',ES9.2,',',ES9.2,'], UTEMP1=',ES9.2)
  1102 format (A, 'RTEMP2=[',ES9.2,',',ES9.2,',',ES9.2,'], UTEMP2=',ES9.2)
  1103 format (A, 'RTEMP3=[',ES9.2,',',ES9.2,',',ES9.2,'], UTEMP3=',ES9.2)
  1104 format (A, 'RTEMP4=[',ES9.2,',',ES9.2,',',ES9.2,'], UTEMP4=',ES9.2)
  1105 format (A, 'RTEMP5=[',ES9.2,',',ES9.2,',',ES9.2,'], UTEMP5=',ES9.2)  
  !1114 format (A, 'AK4=[',ES9.2,',',ES9.2,',',ES9.2,'], UK4=',ES9.2,', GK4=',ES9.2)
  !1124 format (A, 'B=[',ES9.2,',',ES9.2,',',ES9.2,'], E=[',ES9.2,',',ES9.2,',',ES9.2,']')
  
  DO I = 1,3                               !first step
    RTEMP(I) = R(I) + B21*H*DRDT(I)          
  ENDDO
    UTEMP = U + B21*H*DUDT  
  !PRINT 1101, 'RKQS s1:', RTEMP, UTEMP, GAMMATEMP
  CALL DERIVS (T+A2*H, RTEMP, AK2, UTEMP,UK2,MU,T1,T2)     !second step
  DO I = 1,3
    RTEMP(I) = R(I) + H *(B31*DRDT(I)+B32*AK2(I))
  ENDDO
    UTEMP = U+ H *(B31*DUDT +B32*UK2)
!    GAMMATEMP = GAMMA+ H *(B31*DGAMMADT +B32*GK2)
  !PRINT 1102, 'RKQS s2:', RTEMP, UTEMP, GAMMATEMP
  CALL DERIVS (T+A3*H, RTEMP, AK3,UTEMP,UK3,MU,T1,T2)      !third step
  DO I = 1,3
    RTEMP(I) = R(I) + H*(B41*DRDT(I)+B42*AK2(I)+B43*AK3(I))
  ENDDO
    UTEMP = U + H*(B41*DUDT+B42*UK2+B43*UK3)
!    GAMMATEMP = GAMMA + H*(B41*DGAMMADT+B42*GK2+B43*GK3)
  !PRINT 1103, 'RKQS s3:', RTEMP, UTEMP, GAMMATEMP
  CALL DERIVS (T+A4*H, RTEMP, AK4,UTEMP,UK4,MU,T1,T2)   !fourth step
  DO I = 1,3                  
   RTEMP(I) = R(I) + H*(B51*DRDT(I)+B52*AK2(I)+B53*AK3(I) + B54*AK4(I))
  ENDDO
   UTEMP = U + H*(B51*DUDT+B52*UK2+B53*UK3 + B54*UK4)
!   GAMMATEMP = GAMMA + H*(B51*DGAMMADT+B52*GK2+B53*GK3 + B54*GK4)
  !PRINT 1104, 'RKQS s4:', RTEMP, UTEMP, GAMMATEMP
  CALL DERIVS (T+A5*H, RTEMP, AK5, UTEMP, UK5,MU,T1,T2)    !fifth step
  DO I = 1,3                  
   RTEMP(I) = R(I) + H*(B61*DRDT(I)+B62*AK2(I)+B63*AK3(I) + B64*AK4(I)+B65*AK5(I))
  ENDDO
   UTEMP = U + H*(B61*DUDT+B62*UK2+B63*UK3 + B64*UK4+B65*UK5)
!   GAMMATEMP = GAMMA + H*(B61*DGAMMADT+B62*GK2+B63*GK3 + B64*GK4+B65*GK5)
  !PRINT 1105, 'RKQS s5:', RTEMP, UTEMP, GAMMATEMP
  CALL DERIVS (T+A6*H, RTEMP, AK6, UTEMP, UK6,MU,T1,T2)     !sixth step
  DO I = 1,3           !accumulate increments with weights       
   ROUT(I) = R(I) + H*(C1*DRDT(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
  ENDDO
   UOUT = U + H*(C1*DUDT+C3*UK3+C4*UK4+C6*UK6)
!   GAMMAOUT = GAMMA + H*(C1*DGAMMADT+C3*GK3+C4*GK4+C6*GK6)
   !estimate error as difference between 4th and 5th order methods
  DO I= 1,3
    RERR(I)=H*(DC1*DRDT(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*AK6(I))
  ENDDO
    RERR(4)=H*(DC1*DUDT+DC3*UK3+DC4*UK4+DC5*UK5+DC6*UK6)
!    RERR(5)=H*(DC1*DGAMMADT+DC3*GK3+DC4*GK4+DC5*GK5+DC6*GK6) 
  RETURN
  !PRINT*, 'RKQS s6 done'
 END SUBROUTINE RKCK


END MODULE M_rkckR


