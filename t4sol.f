C 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    CODE LAST CHANGED 27 MAY  1981 
C    COMMENTS LAST CHANGED 19 FEBRUARY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SOLSM1 2.2 DATED 8 AUGUST 1980 
C 
        SUBROUTINE SOLSM1(Z,MM,NN,AMAX,LIST,LENG,IFLAG) 
C 
C 
        DIMENSION Z(MM,NN),AMAX(2,LENG),LIST(2,LENG) 
C 
C    Simple user interface subroutine which calls SOLID. 
C 
        IF(MM.LT.3.OR.NN.LT.3)GO TO 2 
C 
C    Compute the range of Z values. 
C 
        ZH = Z(1,1) 
        ZL = ZH 
C 
        DO 1 MS = 1,MM 
        DO 1 NS = 1,NN 
        ZVAL = Z(MS,NS) 
        ZH = AMAX1(ZH,ZVAL) 
    1   ZL = AMIN1(ZL,ZVAL) 
C 
        ZD = ZH-ZL 
C 
C    Create artificial side lengths XL and YL. 
C 
        XL = 3.0*ZD 
        YL = XL*FLOAT(NN-1)/FLOAT(MM-1) 
C 
C    Force SOLID to invent an eye position. 
C 
        XE = -1.0 
        YE = -1.0 
        ZE = -1.0 
C 
C    Call SOLSM2 which will call SOLID. 
C 
        CALL SOLSM2(Z,MM,NN,XL,YL,XE,YE,ZE,AMAX,LIST,LENG,IFLAG) 
        RETURN 
C 
C    Error in MM or NN. 
C 
    2   IFLAG = 3 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SOLSM2 2.1 DATED 7 JULY 1980 
C 
        SUBROUTINE SOLSM2(Z,MM,NN,XL,YL,XE,YE,ZE,AMAX,LIST,LENG,IFLAG) 
C 
C 
        DIMENSION Z(MM,NN),AMAX(2,LENG),LIST(2,LENG) 
        DIMENSION ZX(1,1),ZY(1,1) 
C 
C    Simple user interface routine to SOLID. 
C 
        IF(MM.LT.3.OR.NN.LT.3)GO TO 1 
C 
C    If the user is producing a coarse plot give him some 
C    interpolation. 
C 
        MIX1 = 1 
        NIY1 = 1 
        IF(MM.LT.50)MIX1 = MAX0(3,50/MM) 
        IF(NN.LT.50)NIY1 = MAX0(3,50/NN) 
C 
C    Call SOLID to produce the plot. 
C 
        CALL SOLID(Z,ZX,ZY,MM,NN,1,1,1,1,MM,NN,XL,YL,XE,YE,ZE, 
     1    AMAX,LIST,LENG,LMAX1,IFLAG,-1.0,MIX1,NIY1,1,1, 
     2    .TRUE.,.FALSE.,.FALSE.) 
        RETURN 
C 
C    Error in the value MM or NN. 
C 
    1   IFLAG = 3 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SOLSM3 2.1 DATED 16 MARCH 1982 
C 
        SUBROUTINE SOLSM3(Z,ZX,ZY,MM,NN,XL,YL,XE,YE,ZE,AMAX,LIST, 
     1    LENG,IFLAG,KSEG) 
C 
C    Simple interface routine giving access to gradients and allowing 
C    the specification of inter node curves. 
C 
        DIMENSION Z(MM,NN),ZX(MM,NN),ZY(MM,NN),AMAX(2,LENG),LIST(2,LENG) 
C 
        CALL SOLID(Z,ZX,ZY,MM,NN,MM,NN,MM,NN,MM,NN,XL,YL,XE,YE,ZE, 
     1    AMAX,LIST,LENG,LMAX1,IFLAG,-1.0,KSEG,KSEG,1,1,.TRUE., 
     2    .TRUE.,.TRUE.) 
C 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    FUNCTION IFUZ 2.1 DATED 9 DECEMBER 1978 
C 
        FUNCTION IFUZ(VAL) 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
C 
C    This function is used to decide when the result of 
C    a floating point operation is close enough to zero to 
C    be considered zero. 
C 
C    IFUZ is set to: 
C          0 if -FUZ<VAL<FUZ 
C          1 if VAL>FUZ 
C         -1 if VAL<-FUZ 
C 
        IF(ABS(VAL).GT.FUZ)GO TO 1 
    3   IFUZ = 0 
        RETURN 
    1   IF(VAL)2,3,4 
    2   IFUZ = -1 
        RETURN 
    4   IFUZ = 1 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SEARCH 2.1  DATED 9 DECEMBER 1978 
C 
        SUBROUTINE SEARCH(AMAX,LIST,LENG,P,LI,LJ) 
C 
C 
        DIMENSION AMAX(2,LENG),LIST(2,LENG) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C 
C    This subroutine finds the segment in the C.H.P. In which P lies. 
C    LGUESS is a guess at an index of a point near P. 
C    If LGUESS is not the index of a valid point a search 
C    from the end of the profile nearest P is performed. 
C    If it is the search starts from LGUESS. 
C 
C 
C    The value of IFLAG1 on return indicates various conditions: 
C 
C       1  LJ is to the J of the point and LI to the I. 
C       2  point coincides with PI of the C.H.P. 
C          PJ set equal to PI 
C       3  point coincides with the lh end of the C.H.P. 
C       4  point coincides with the rh end of the C.H.P. 
C       5  point is off lh end 
C       6  point is off rh end 
C 
C 
C    Find the distance to the lh and rh ends of the C.H.P. 
C 
        DL = P-AMAX(1,LSTART) 
        DR = AMAX(1,LEND)-P 
C 
C    Check that P is within the correct range. 
C 
        IF(IFUZ(DL))1,2,3 
    3   IF(IFUZ(DR))4,5,6 
C 
C    If LGUESS is a valid location use it. 
C 
    6   IF(LGUESS.LE.0.OR.LGUESS.GT.LENG.OR. 
     1    LGUESS.EQ.LSTART.OR.LGUESS.EQ.LEND)GO TO 7 
        IF(LIST(1,LGUESS).GT.0)GO TO 8 
C 
C    Chose whether to start from the right or left hand end of the C.H.P. 
C 
    7   IF(DL.GT.DR)GO TO 9 
C 
        L1 = LSTART 
C 
C    Start from the left hand end. 
C 
C 
C    Move to the right. 
C 
   10   L1 = LIST(1,L1) 
        IF(L1)11,4,13 
   11   IFLAG1 = 99 
        RETURN 
   13   DL = P-AMAX(1,L1) 
C 
C    Found it? 
C 
        IF(IFUZ(DL))14,15,10 
   14   IF(I.EQ.2)GO TO 16 
C 
C 
C    Point lies between LIST(2,L1) and L1. 
C 
        LI = L1 
        LJ = LIST(2,L1) 
        IFLAG1 = 1 
        RETURN 
C 
C 
C    Point lies between L1 and LIST(2,L1) 
C 
   16   LJ = L1 
        LI = LIST(2,L1) 
        IFLAG1 = 1 
        RETURN 
C 
C    Start from right hand end. 
C 
    9   L1 = LEND 
C 
C    Move to the left 
C 
   17   L1 = LIST(2,L1) 
        IF(L1)18,1,20 
   18   IFLAG1 = 98 
        RETURN 
   20   DR = AMAX(1,L1)-P 
C 
C    Found it? 
C 
        IF(IFUZ(DR))21,15,17 
   21   IF(I.EQ.2)GO TO 22 
C 
C    Point lies between L1 and LIST(1,L1) 
C 
        LJ = L1 
        LI = LIST(1,L1) 
        IFLAG1 = 1 
        RETURN 
C 
C    Point lies between LIST(1,L1) and L1. 
C 
   22   LI = L1 
        LJ = LIST(1,L1) 
        IFLAG1 = 1 
        RETURN 
C 
C    Start from LGUESS. 
C 
    8   L1 = LGUESS 
C 
C    Decide whether to move left or right. 
C 
        DL = P-AMAX(1,L1) 
        IF(IFUZ(DL))17,15,10 
C 
C    Point is off the lh end. 
C 
    1   IFLAG1 = 5 
        LI = LSTART 
        LJ = LSTART 
        RETURN 
C 
C    Point is off the rh end. 
C 
    4   IFLAG1 = 6 
        LI = LEND 
        LJ = LEND 
        RETURN 
C 
C    Point coincides with the lh end. 
C 
    2   IFLAG1 = 3 
        LI = LSTART 
        LJ = LSTART 
        RETURN 
C 
C    Point coincides with the rh end. 
C 
    5   IFLAG1 = 4 
        LI = LEND 
        LJ = LEND 
        RETURN 
C 
C    Point coincides with L1. 
C 
   15   IFLAG1 = 2 
        LI = L1 
        LJ = L1 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE ACROSS 2.1 DATED 9 DECEMBER 1978 
C 
        SUBROUTINE ACROSS(AMAX,LENG,LI,LJ,PX,QX,IDG) 
C 
C 
        DIMENSION AMAX(2,LENG) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C 
C    This subroutine finds the intersection of the line (PI,QI),(PJ,QJ) with 
C    the segment LI to LJ of the C.H.P. and sets IFLAG1 and IDG as follows: 
C 
C    IFLAG1 
C      1  intersection in the segment LJ-LI 
C      2        "       coincides with LI 
C      3        "           "       "  LJ 
C      4        "      to I of segment 
C      5        "       " J "     " 
C      6  lines are parallel. 
C 
C    IDG 
C      1  no parallel lines or coincidences. 
C      2  PI,QI coincide with line. 
C      3  PJ,QJ coincide with line. 
C      4  parallel above 
C      5     "     below 
C      6  colinear 
C 
        IDG = 1 
C 
C    Find the gradient of the C.H.P. segment LJ-LI. 
C 
        PRI = AMAX(1,LI) 
        PRJ = AMAX(1,LJ) 
        A1 = PRI-PRJ 
C 
C    Test for a vertical segment. 
C 
        IF(IFUZ(10.0*A1))1,2,1 
C 
C    Section to deal with a vertical segment. 
C 
    2   PX = PRI 
        QX = A*PRI+B 
        IF(IFUZ(QX-AMAX(2,LJ)))50,51,52 
   52   IF(IFUZ(QX-AMAX(2,LI)))53,54,55 
   53   IFLAG1 = 1 
   56   IF(IFUZ(PX-PI).EQ.0)IDG = 2 
        IF(IFUZ(PX-PJ).EQ.0)IDG = 3 
        RETURN 
   51   IFLAG1 = 3 
        GO TO 56 
   50   IFLAG1 = 5 
        RETURN 
   54   IFLAG1 = 2 
        GO TO 56 
   55   IFLAG1 = 4 
        RETURN 
C 
C    The segment is not vertical.  Carry on with the gradient 
C    computation. 
C 
    1   QRI = AMAX(2,LI) 
        QRJ = AMAX(2,LJ) 
        A1 = (QRI-QRJ)/A1 
C 
C    Find the difference in gradients and the constant term in the 
C    equation of the C.H.P. segment. 
C 
        GD = A-A1 
        B1 = QRI-A1*PRI 
        IF(IFUZ(GD))5,6,5 
C 
C    The lines are parallel. 
C 
    6   PX = 0.0 
        QX = 0.0 
        IFLAG1 = 6 
        QRI1 = A*PRI+B 
C 
C    Test to see if the segment is above or below the new line. 
C 
        IF(IFUZ(QRI1-QRI))8,9,10 
C 
C    The lines are colinear. 
C 
    9   IDG = 6 
        IFLAG1 = 3 
        PX = AMAX(1,LJ) 
        QX = AMAX(2,LJ) 
        RETURN 
C 
C    The lines are parallel - the new line below. 
C 
    8   IDG = 5 
        RETURN 
C 
C    The lines are parallel - the new line above. 
C 
   10   IDG = 4 
        RETURN 
C 
C 
C    Find the intersection of the lines , PX,QX. 
C 
    5   PX = (B1-B)/GD 
        QX = (A*B1-A1*B)/GD 
        XM = FLOAT(J-I) 
C 
C    Find the relation of PX to PJ-PI. 
C 
        IF(IFUZ(PX-PJ))14,13,14 
C 
C    PJ coincides with PX. 
C 
   13   IDG = 3 
C 
C    PI coincides with PX. 
C 
   14   IF(IFUZ(PX-PI))16,15,16 
   15   IDG = 2 
C 
C    Find out if the intersection is to the J or the I of the segment. 
C 
   16   IF(IFUZ((PX-PRJ)*XM))17,18,19 
C 
C    Intersection coincides with the J end of the segment. 
C 
   18   IFLAG1 = 3 
        RETURN 
   19   IF(IFUZ((PRI-PX)*XM))20,21,22 
C 
C    Intersection coincides with the I end of the segment. 
C 
   21   IFLAG1 = 2 
        RETURN 
   22   IFLAG1 = 1 
        RETURN 
C 
C    Intersection is off the J end of the segment. 
C 
   17   IFLAG1 = 5 
        RETURN 
C 
C    Intersection is off the I end of the segment. 
C 
   20   IFLAG1 = 4 
        RETURN 
        END 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE ABOVE 2.1 DATED 9 DECEMBER 1978 
C 
C 
        SUBROUTINE ABOVE(AMAX,LENG,LI,LJ,P,Q) 
C 
C 
        DIMENSION AMAX(2,LENG) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C 
C    Subroutine to find out if (P,Q) is above the segment LI-LJ of 
C    the C.H.P. and sets IFLAG1 as follows: 
C 
C    IFLAG1 
C     1  above 
C     2  below 
C     3  coincides 
C     4  outside segment on J side 
C     5     "        "   "  I  " 
C 
C 
C    Find the ends of the segment. 
C 
        XM = FLOAT(J-I) 
        PRI = AMAX(1,LI) 
        PRJ = AMAX(1,LJ) 
        QRI = AMAX(2,LI) 
        QRJ = AMAX(2,LJ) 
C 
C    Find the equation of the segment. 
C 
        A1 = PRI-PRJ 
C 
C    Test for a vertical segment. 
C 
        IF(IFUZ(10.0*A1))11,10,11 
C 
C    Section to deal with a vertical segment. 
C 
   10   IF(IFUZ(P-PRJ))50,51,52 
   52   IFLAG1 = 3+J 
        RETURN 
   50   IFLAG1 = 3+I 
        RETURN 
   51   IF(IFUZ(Q-QRJ))53,54,55 
   54   IFLAG1 = 3 
        RETURN 
   53   IF(IFUZ(Q-QRI))56,54,54 
   56   IFLAG1 = 2 
        RETURN 
   55   IF(IFUZ(Q-QRI))54,54,57 
   57   IFLAG1 = 1 
        RETURN 
C 
C    Section to deal with a non-vertical segment. 
C 
   11   A1 = (QRI-QRJ)/A1 
        B1 = QRI-A1*PRI 
C 
C    Find the relationship of P to the J end. 
C 
        IF(IFUZ((P-PRJ)*XM))1,2,3 
    2   IF(IFUZ(Q-QRJ))4,5,6 
    4   IFLAG1 = 2 
        RETURN 
    5   IFLAG1 = 3 
        RETURN 
    6   IFLAG1 = 1 
        RETURN 
    1   IFLAG1 = 4 
        RETURN 
C 
C    Find the relationship of P to the I end. 
C 
    3   IF(IFUZ((PRI-P)*XM))7,8,9 
    8   IF(IFUZ(Q-QRI))4,5,6 
    7   IFLAG1 = 5 
        RETURN 
C 
C    Find the point on the segment with the value P. 
C 
    9   QR = A1*P+B1 
C 
C    Find out if it is above or below Q. 
C 
        IF(IFUZ(Q-QR))4,5,6 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SOLID 2.6 DATED 27 MAY 1981 
C 
C 
        SUBROUTINE SOLID(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,M,N, 
     1    XL,YL,XE,YE,ZE,AMAX,LIST,LENG,LMAX1,IFLAG,HAZE1, 
     2    MIX1,NIY1,MSX1,NSY1,IQUAD1,IGDX1,IGDY1) 
C 
C 
        DIMENSION Z(MM,NN),ZX(MMX,NNX),ZY(MMY,NNY),AMAX(2,LENG), 
     1    LIST(2,LENG) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD,IQUAD1,IGDX1,IGDY1 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
        DATA CFUZ/1.0E-5/ 
C 
C    This is the master subroutine in the SOLID package. 
C    The arguments are as follows (an r in the left 
C    hand margin indicates a read only argument, a u indicates 
C    an argument that must be set before the routine is called): 
C 
C     r Z      Real array dimension Z(MM,NN).  This should contain 
C              the rectangular matrix of heights to be plotted unless 
C              the surface is defined by the function ZED (see below). 
C 
C     r ZX     Real array dimension ZX(MMX,NNX). This should 
C              contain the partial derivatives of the surface 
C              with respect to x (the x direction corresponds to 
C              increasing the first subscript of Z). 
C 
C     r ZY     Real array dimension ZY(MMY,NNY).  This should 
C              contain the partial derivatives of the surface 
C              with respect to y (the y direction corresponds to 
C              increasing the second subscript of Z). 
C 
C    ur MM     Integer variable.  The true first dimension of the 
C              array Z. 
C 
C    ur NN     Integer variable.  The true second dimension of the 
C              array Z. 
C 
C    ur MMX    Integer variable.  The true first dimension of the 
C              array ZX.  This must be equal to MM or 1 (see below). 
C 
C    ur NNX    Integer variable.  The true second dimension of the 
C              array ZX.  This must equal NN or 1 (see below). 
C 
C    ur MMY    Integer variable.  The true first dimension of the 
C              array ZY.  This must equal MM or 1 (see below). 
C 
C    ur NNY    Integer variable.  The true second dimension of the 
C              array ZY.  This must equal NN or 1 (see below). 
C 
C    ur M      Integer variable. The first dimension of the arrays 
C              Z,ZX and ZY that is actually to be plotted.  If MSX1 
C              is greater than 1 M is the number of rectangular grid 
C              elements plus 1 to be plotted in the x direction. 
C              ((M-1)*MSX1+1) must be less than or equal to MM unless 
C              MM is 1 (see below). 
C 
C    ur N      Integer variable.  The second dimension of the arrays 
C              Z,ZX and ZY that is actually to be plotted.  If NSY1 
C              is greater than 1 N is the number of rectangular grid 
C              elements plus 1 to be plotted in the y direction. 
C              ((N-1)*NSY1+1) must be less than or equal to NN 
C              unless NN is 1 (see below). 
C 
C    ur XL     Real variable.  The length in the x direction of the 
C              rectangular matrix of heights to be plotted.  XL 
C              must be in the same units as Z. 
C 
C    ur YL     Real variable.  The length in the y direction of the 
C              rectangular matrix of heights to be plotted.  YL 
C              must be in the same units as Z. 
C 
C    u  XE     Real variable.  The x coordinate of the eye position 
C              from which the surface is to be viewed.  XE must be in 
C              the same units as Z. 
C 
C    u  YE     Real variable.  The y coordinate of the eye position 
C              from which the surface is to be viewed.  YE must 
C              be in the same units as Z. 
C 
C    u  ZE     Real variable.  The z coordinate of the eye position 
C              from which the surface is to be viewed.  ZE must be 
C              in the same units as Z. 
C 
C       AMAX   Real array dimension AMAX(2,LENG).  Used as workspace. 
C 
C       LIST   Integer array dimension LIST(2,LENG).  Used as 
C              workspace. 
C 
C    ur LENG   Integer variable. The second dimension of the 
C              arrays AMAX and LIST.  The value of LENG needed 
C              for a given surface will depend upon how finely 
C              it is defined by the height matrix and any 
C              interpolation and upon how rough it is.  In 
C              general a value of 
C 
C                     LENG = 2*(M*MIX1+N*NIY1) 
C 
C              will prove sufficient. 
C 
C       LMAX1  Integer variable.  On successsful completion 
C              of the plot LMAX1 will contain the maximum 
C              amount of the arrays AMAX and LIST actually 
C              used.  If this is substantially less than 
C              LENG the value of LENG may be reduced on 
C              subsequent calls to SOLID for plotting the 
C              same surface to a value of about, say 
C 
C                 LENG = INT(1.5*FLOAT(LMAX1)) 
C 
C 
C       IFLAG  Integer variable.  On successful completion of the 
C              plot IFLAG will contain 0.  Non-zero values 
C              indicate various failures: 
C 
C                  IFLAG            ERROR 
C 
C                    1              Out of space in the C.H.P. 
C                                   (LENG too small) 
C                    2              Illegal eye coordinate(s) 
C                    3              Illegal value in argument(s) 
C 
C 
C 
C    ur HAZE1  Real variable.  To be used for future enhancements 
C              Should be set to -1.0. 
C 
C    ur MIX1   Integer variable.  The number of straight line 
C              segments to interpolate between two grid points 
C              in the x direction to give smooth curves.  MIX1 
C              should be greater than or equal to 1. 
C 
C    ur NIY1   Integer variable.  The number of straight line 
C              segments to interpolate between two grid points 
C              in the y direction to give smooth curves.  NIY1 
C              should be greater than or equal to 1. 
C 
C    ur MSX1   Integer variable. The number of grid points to be plotted 
C              in the x direction to form an edge of a rectangular grid 
C              element.  Usually this will be 1.  Setting MSX1 to a 
C              value greater than 1 allows the user to define the 
C              surface more finely than the rectangular grid elements 
C              that will make up his plot.  If heights are being 
C              supplied in the array Z some of this array, of necesity, 
C              will be unused if MSX1 is greater than 1. 
C 
C    ur NSY1   Integer variable. The number of grid points to be plotted 
C              in the y direction to form an edge of a rectangular grid 
C              element.  Usually this will be 1.  Setting NSY1 to a 
C              value greater than 1 allows users to define the surface 
C              more finely than the rectangular grid elements that will 
C              make up his plot.  If heights are being supplied in the 
C              array Z some of this array will be unused if 
C              NSY1 is greater than 1. 
C 
C    ur IQUAD  Logical variable.  Controls the type of 
C              interpolation used in the drawing of smooth 
C              curves.  .TRUE. gives seamed quadratic 
C              interpolation, .FALSE. gives cubic interpolation. 
C              If no interpolation is needed (MIX1.EQ.1.AND.NIY1.EQ.1) 
C              then IQUAD can be set to either value. 
C 
C    ur IGDX1  Logical variable.  If the user is supplying SOLID 
C              with partial derivatives of the surface in the array 
C              ZX or from the function ZEDX he should set this 
C              variable .TRUE.  Othewise it must be set .FALSE. 
C              (see below). 
C 
C    ur IGDY1  Logical variable.  If the user is supplying SOLID 
C              with partial derivatives of the surface in the array 
C              ZY or from the function ZEDY he should set this 
C              variable .TRUE.  Othewise it should be set .FALSE. 
C              (see below). 
C 
C 
C 
C 
C    The eye position (XE,YE,ZE) must be such that the eye is 
C    not above the rectangle 
C 
C                    |x| < 0.5*XL 
C                    |y| < 0.5*YL 
C 
C    The z coordinate, ZE, should, in general, be positive.  However, 
C    if ZE is negative an automatic viewing position will be 
C    created using it and the sign of XE and YE and the lengths XL 
C    and YL as follows: 
C 
C      1)  The viewing position will be looking at a corner. 
C 
C      2)  XE will be set to SIGN(3.0*XL,XE) 
C          YE  "    "  "  "  SIGN(2.0*YL,YE) 
C          ZE  "    "  "  "  ABS(ZE)*0.3*(ABS(XE)+ABS(YE)) 
C 
C    The effect of this is to produce a reasonable picture for 
C    most surfaces if (XE,YE,ZE) is set to (-1.0,-1.0,-1.0) 
C    for example.  The eye position generated is returned in 
C    XE, YE and ZE, so that subsequent calls of SOLID will 
C    produce a view from the same position if these variables 
C    are left unaltered by the user's program. 
C 
C    The reasons for dimensioning the arrays Z, ZX and ZY separately 
C    are as follows: 
C 
C      1)  If MIX1 (or NIY1) is 1 the user wants straight lines 
C          plotted between the grid points in the x (or y) direction 
C          and the partial derivative of the surface with respect 
C          to x (or y) used in the fitting of smooth curves between 
C          grid points is not needed, so the array space used to 
C          store the information can be saved by dimensioning 
C          ZX to ZX(1,1) (or ZY to ZY(1,1)). 
C 
C      2)  If the user wants smooth curves between grid points but 
C          can not compute the gradient(s) he should dimension the 
C          gradient array(s) to (1,1) and set MIX1 and NIY1 
C          appropriately to  a value or values greater than 1. 
C          The routine will then fit gradients to the height 
C          matrix using parabolas and the smooth curves will be drawn. 
C 
C      3)  If MM and NN are set to 1 the routine assumes that the 
C          user wishes to plot a surface that he has written in the 
C          function 
C 
C                      FUNCTION ZED(MS,NS) 
C 
C          The two integers MS and NS are indexed between 1 to M and 
C          1 to N by the SOLID routine. ZED should return a height value 
C          in the same way as one might be set in the array Z. 
C          If the user decides upon this course of action and wishes 
C          to provide the SOLID routine with gradients for the fitting 
C          of smooth curves through the grid points he can also 
C          supply two functions for the computation of the partial 
C          derivatives of the surface with respect to x and y: 
C 
C                      FUNCTION ZEDX(MS,NS) 
C                      FUNCTION ZEDY(MS,NS) 
C 
C          If he does not, and sets IGDX1 and IGDY1 to .FALSE., SOLID 
C          will fit parabolas to the function ZED in order to compute 
C          gradients.  This will mean that ZED will be referenced a 
C          large number of times, and the user should therefor take 
C          steps to ensure that it is as efficient as possible. 
C 
C 
C 
C 
C 
C 
C    Check for erroneous subroutine arguments. 
C 
        IF(MM.LE.0.OR.NN.LE.0.OR.MMX.LE.0.OR.NNX.LE.0.OR.MMY. 
     1    LE.0.OR.NNY.LE.0.OR.LENG.LE.0)GO TO 9999 
        IF(M.LE.2.OR.N.LE.2.)GO TO 9999 
        IF(XL.LE.0.0.OR.YL.LE.0.0)GO TO 9999 
        IF(MIX1.LE.0.OR.NIY1.LE.0)GO TO 9999 
        IF(MSX1.LE.0.OR.NSY1.LE.0)GO TO 9999 
C 
C    To be removed..... 
C 
        IF(HAZE1.GE.0.0)GO TO 9999 
C 
C    Record the size of the plot grid in MTOP and NTOP. 
C 
        MTOP = M 
        NTOP = N 
C 
C    Set up the logical variables in the common block. 
C 
        IGRADX = .TRUE. 
        IF(IGDX1) IGRADX = .FALSE. 
        IGRADY = .TRUE. 
        IF(IGDY1) IGRADY = .FALSE. 
        IHAZE = .FALSE. 
        IFUN = .FALSE. 
        IF(MM.EQ.1.AND.NN.EQ.1) IFUN = .TRUE. 
        IFLAG1 = 0 
        LDJ = -1 
        M2 = MM 
        IF(IFUN)M2 = (M-1)*MSX1+1 
        IF(.NOT.IFUN.AND.(((M-1)*MSX1+1).GT.MM.OR.((N-1)*NSY1+1).GT.NN)) 
     1    GO TO 9999 
        IQUAD = IQUAD1 
C 
C    Find the minimum z value (the surface is plotted so that all z 
C    values appear positive, the minimum one lying in the plane z = 0.0) 
C    also find the maximum z value so that the plot window needed 
C    may be subsequently calculated. 
C 
C 
C    Set up the values in the common block so that the routine ZED9 
C    considers the surface to be being observed by an eye facing 
C    the corner Z(1,1). 
C 
        K(1) = 0 
        K(2) = 1 
        K(3) = -1 
        K(4) = 1 
        K(5) = 0 
        K(6) = 0 
        JGRAD = 1 
        ZMIN = 0.0 
C 
C    Set up the variables needed for interpolation. 
C 
        MSX = MSX1 
        NSY = NSY1 
        MIX2 = MIX1/2 
        NIY2 = NIY1/2 
        NODX = MOD(MIX1,2) 
        NODY = MOD(NIY1,2) 
        MLIM = (M-1)*MSX 
        NLIM = (N-1)*NSY 
        MT1 = MLIM+1 
        NT1 = NLIM+1 
        XINCB = XL/FLOAT(MLIM) 
        YINCB = YL/FLOAT(NLIM) 
        XINCS = XINCB/FLOAT(MIX1) 
        YINCS = YINCB/FLOAT(NIY1) 
        XHAL = 0.5*XINCB 
        YHAL = 0.5*YINCB 
        RXHAL = 1.0/XHAL 
        RYHAL = 1.0/YHAL 
C 
C    Set up the x and y coordinates of a point one distance increment 
C    (XINCB and YINCB) away from the corner Z(1,1). 
C 
        X0 = -XL*0.5-XINCB 
        Y0 = -YL*0.5-YINCB 
        MS = 1 
        NS = 1 
C 
C    Find the first height. 
C 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        ZTEM = ZZ 
        ZMAX = ZZ 
        XDIR = 1.0 
        YDIR = 1.0 
C 
C    Is interpolation in the y direction neeeded? 
C 
        IF(NIY2.LE.0)GO TO 84 
C 
C    Yes.  Raster scan in the y direction for the maximum 
C    and minimum z values. 
C 
        DO 77 MF = 1,MT1,MSX 
        MS = MF 
        Y = Y0+YINCB 
        NF = 1 
        NS = 1 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
        DO 77 NF = 1,NLIM 
        NS = NF+1 
        CALL YPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,Y1,Y2,YDIR) 
        YOFF = -YHAL+YINCS 
        DO 78 NIT = 1,NIY2 
        ZZ = POLYA*(YOFF**3)+POLYB1*YOFF*YOFF+POLYC*YOFF+POLYD 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
   78   YOFF = YOFF+YINCS 
        NIT2 = NIY2+NODY 
        DO 79 NIT = 1,NIT2 
        ZZ = POLYA*(YOFF**3)+POLYB2*YOFF*YOFF+POLYC*YOFF+POLYD 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
   79   YOFF = YOFF+YINCS 
   77   Y = Y2 
C 
C      Is interpolation in the x direction required? 
C 
        IF(MIX2.LE.0)GO TO 90 
        GO TO 86 
   84   IF(MIX2.LE.0)GO TO 85 
C 
C    Yes.  Raster scan in the x dirction for the maximum 
C    and minimum z values. 
C 
   86   DO 87 NF = 1,NT1,NSY 
        NS = NF 
        X = X0+XINCB 
        MF = 1 
        MS = 1 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
        DO 87 MF = 1,MLIM 
        MS = MF+1 
        CALL XPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,X1,X2,XDIR) 
        XOFF = -XHAL+XINCS 
        DO 88 NIT = 1,MIX2 
        ZZ = POLYA*(XOFF**3)+POLYB1*XOFF*XOFF+POLYC*XOFF+POLYD 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
   88   XOFF = XOFF+XINCS 
        NIT2 = MIX2+NODX 
        DO 89 NIT = 1,NIT2 
        ZZ = POLYA*(XOFF**3)+POLYB2*XOFF*XOFF+POLYC*XOFF+POLYD 
        ZMAX = AMAX1(ZZ,ZMAX) 
        ZTEM = AMIN1(ZZ,ZTEM) 
   89   XOFF = XOFF+XINCS 
   87   X = X2 
        GO TO 90 
C 
C    No interpolation is being carried out.  Scan the height 
C    matrix or function for the extreme values. 
C 
   85   IF(IFUN)GO TO 91 
C 
C    Scan the matrix. 
C 
        DO 82 MS = 1,M 
        DO 82 NS = 1,N 
        ZMAX = AMAX1(ZMAX,Z(MS,NS)) 
   82   ZTEM = AMIN1(ZTEM,Z(MS,NS)) 
        GO TO 90 
C 
C    Scan the function. 
C 
   91   DO 93 MS = 1,M 
        DO 93 NS = 1,N 
        ZZ = ZED(MS,NS) 
        ZMAX = AMAX1(ZMAX,ZZ) 
   93   ZTEM = AMIN1(ZTEM,ZZ) 
C 
C    Transfer the minimum z value to ZMIN. 
C 
   90   ZMIN = ZTEM 
        ZD = ZMAX-ZMIN 
C 
C    Find out if the user wants automatic viewing and create 
C    an eye position if need be. 
C 
   80   IF(ZE)74,73,73 
   74   XE = SIGN(3.0*XL,XE) 
        YE = SIGN(2.0*YL,YE) 
        ZE = ABS(ZE)*0.3*(ABS(XE)+ABS(YE)) 
   73   XL2 = XL*0.5 
        YL2 = YL*0.5 
        ZE1 = ZE 
C 
C    Find the eye position and use it to decide upon the rotation 
C    of the height matrix or function.  The values stored 
C    in the arrays K and JGRAD are used by the routine ZED9 
C    to rotate the entries in the arrays Z, ZX and ZY or the 
C    functions ZED, ZEDX and ZEDY. 
C 
        IF(XE.LT.-XL2.AND.YE.LT.YL2)GO TO 1 
        IF(XE.LT.XL2.AND.YE.GT.YL2)GO TO 2 
        IF(XE.GT.XL2.AND.YE.GT.-YL2)GO TO 3 
        IF(XE.GT.-XL2.AND.YE.LT.-YL2)GO TO 4 
        GO TO 5 
C 
C 
C    Eye is bottom left. 
C    Some of the rotation dependent values in the 
C    common block have already been set to deal with this eventuality. 
C 
    1   XE1 = XE 
        YE1 = YE 
        MIX = MIX1 
        NIY = NIY1 
        M1 = M 
        N1 = N 
        XL1 = XL 
        YL1 = YL 
        CALL PSPSET(XE1,YE1,ZE,0.0,0.0,ZE,0.0,0.0,ZE*1.1) 
        IF(IFLAG1.GT.50)GO TO 5 
        GO TO 6 
C 
C    Eye is top left. 
C 
    2   K(1) = -1 
        K(2) = 0 
        K(3) = N 
        K(4) = 0 
        K(5) = 1 
        K(6) = 0 
        M1 = N 
        N1 = M 
        MSX = NSY1 
        NSY = MSX1 
        MIX = NIY1 
        NIY = MIX1 
        JGRAD = 2 
        XE1 = -YE 
        YE1 = XE 
        XL1 = YL 
        YL1 = XL 
        CALL PSPSET(XE1,YE1,ZE,0.0,0.0,ZE,0.0,0.0,ZE*1.1) 
        IF(IFLAG1.GT.50)GO TO 5 
        GO TO 6 
C 
C    Eye is top right. 
C 
    3   K(1) = 0 
        K(2) = -1 
        K(3) = N 
        K(4) = -1 
        K(5) = 0 
        K(6) = (M-1)*MSX1+2 
        M1 = M 
        N1 = N 
        MSX = MSX1 
        NSY = NSY1 
        MIX = MIX1 
        NIY = NIY1 
        JGRAD = 3 
        XE1 = -XE 
        YE1 = -YE 
        XL1 = XL 
        YL1 = YL 
        CALL PSPSET(XE1,YE1,ZE,0.0,0.0,ZE,0.0,0.0,ZE*1.1) 
        IF(IFLAG1.GT.50)GO TO 5 
        GO TO 6 
C 
C    Eye is bottom right. 
C 
    4   K(1) = 1 
        K(2) = 0 
        K(3) = -1 
        K(4) = 0 
        K(5) = -1 
        K(6) = (M-1)*MSX1+2 
        M1 = N 
        N1 = M 
        MSX = NSY1 
        NSY = MSX1 
        MIX = NIY1 
        NIY = MIX1 
        JGRAD = 4 
        XE1 = YE 
        YE1 = -XE 
        XL1 = YL 
        YL1 = XL 
        CALL PSPSET(XE1,YE1,ZE,0.0,0.0,ZE,0.0,0.0,ZE*1.1) 
        IF(IFLAG1.GT.50)GO TO 5 
C 
C    The relationship between the eye position and the surface has now 
C    been established.  From this point on the routines will pretend to 
C    be looking at the surface from a position to the left (and 
C    possibly below)Z(1,1) (or ZED(1,1)).  The surface will be rotated 
C    by the routine ZED9 so that this gives the required picture. 
C 
C    Check that sufficient space exists in the C.H.P. 
C 
    6   IF(LENG.LE.N1*MIX+5)GO TO 9998 
C 
C    Set up the lower left corner (X0,Y0) and the axis increments. 
C 
        MLIM = (M1-1)*MSX 
        NLIM = (N1-1)*NSY 
        MT1 = MLIM+1 
        NT1 = NLIM+1 
        XINCB = XL1/FLOAT(MLIM) 
        YINCB = YL1/FLOAT(NLIM) 
        XINCS = XINCB/FLOAT(MIX) 
        YINCS = YINCB/FLOAT(NIY) 
        XHAL = 0.5*XINCB 
        YHAL = 0.5*YINCB 
        RXHAL = 1.0/XHAL 
        RYHAL = 1.0/YHAL 
        MIX2 = MIX/2 
        NIY2 = NIY/2 
        NODX = MOD(MIX,2) 
        NODY = MOD(NIY,2) 
        X0 = -XL1*0.5-XINCB 
        Y0 = -YL1*0.5-YINCB 
C 
C    Find the enclosing rectangle around the perspective plot by 
C    inspecting the perspective transformation of the minimum 
C    enclosing cuboid around the surface.  Also compute the squared 
C    distances of the nearest and furthest corners of that 
C    cuboid to the eye position. 
C 
        MS = 1 
        NS = 1 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        X = X0+XINCB 
        Y = Y0+YINCB 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 5 
        IF(IFLAG1.GT.50)GO TO 5 
        PMIN = P 
        QMIN = Q 
        PMAX = P 
        QMAX = Q 
        D0 = (XE1-X)**2+(YE1-Y)**2+(ZE1-ZZ)**2 
        DD = D0 
        DO 110 KK = 1,3,2 
        KK1 = KK-2 
        X = FLOAT(KK1)*(X0+XINCB) 
        DO 110 KKK = 1,3,2 
        KK2 = KKK-2 
        Y = FLOAT(KK2)*(Y0+YINCB) 
        DO 110 KKKK = 1,2 
        KK3 = KKKK-1 
        ZZ = FLOAT(KK3)*ZD 
        D = (XE1-X)**2+(YE1-Y)**2+(ZE1-ZZ)**2 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 5 
        IF(IFLAG1.GT.50)GO TO 5 
        PMIN = AMIN1(P,PMIN) 
        PMAX = AMAX1(P,PMAX) 
        QMIN = AMIN1(Q,QMIN) 
        QMAX = AMAX1(Q,QMAX) 
        D0 = AMIN1(D0,D) 
  110   DD = AMAX1(DD,D) 
C 
C 
C    If aerial perspective is required set it up. 
C    (This will be a future enhancement). 
C 
C       IF(HAZE1.GT.0.0)CALL HAZSET(HAZE1) 
C 
C    Start the plot. 
C 
        PD = 0.01*(PMAX-PMIN) 
        QD = 0.01*(QMAX-QMIN) 
        CALL PLTON(PMIN-PD,PMAX+PD,QMIN-QD,QMAX+QD) 
C 
C    Decide what is meant by zero. 
C 
        FUZ = CFUZ*SQRT((PMAX-PMIN)**2+(QMAX-QMIN)**2) 
C 
C 
C    The first point to be plotted is at the maximum y 
C    coordinate and the minimum x coordinate. 
C 
        X = X0+XINCB 
        Y = -(Y0+YINCB) 
        ZZ = 0.0 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 5 
        CALL HDRAW(P,Q,0) 
        POLD = P 
        QOLD = Q 
        PJ = P 
        QJ = Q 
        NS = NT1 
        MS = 1 
        MF = 1 
        NF = NT1 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 5 
C 
C    Set up an artificial C.H.P. to build upon. 
C 
        AMAX(1,1) = PMIN-PD 
        AMAX(2,1) = Q 
        LIST(1,1) = 2 
        LIST(2,1) = 0 
        AMAX(1,2) = P 
        AMAX(2,2) = Q 
        LIST(1,2) = 3 
        LIST(2,2) = 1 
        LSTART = 1 
        LDJ = 2 
        AMAX(1,3) = P 
        AMAX(2,3) = QMIN-QD 
        AMAX(1,4) = PMAX+PD 
        AMAX(2,4) = QMIN-QD 
        LIST(1,3) = 4 
        LIST(2,3) = 2 
        LIST(1,4) = 0 
        LIST(2,4) = 3 
        LFREE = 5 
        DO 798 L = LFREE,LENG 
  798   LIST(1,L) = -L-1 
        LIN = 4 
        LMAX = 4 
        LEND = 4 
        LSTART = 1 
        LGUESS = 1 
        I = 1 
        J = 2 
        CALL HDRAW(P,Q,1) 
C 
C 
C    Scan along the x = X0 face drawing it and modifying 
C    the C.H.P. 
C 
        DO 440 NT = 2,N1 
        DO 441 NTT = 1,NSY 
        NS = NS-1 
  441   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        ZTEM = ZZ 
        ZZ = 0.0 
        CALL PSP(P1,Q1) 
        IF(IFLAG1.GT.50)GO TO 5 
        CALL HDRAW(P1,Q1,1) 
        CALL HDRAW(POLD,QOLD,1) 
        CALL HDRAW(PJ,QJ,0) 
        ZZ = ZTEM 
        POLD = P1 
  440   QOLD = Q1 
C 
C    Is the eye facing a corner? 
C 
        IF(IFUZ(Y0+YINCB-YE1).LE.0)GO TO 75 
C 
C    Yes. Eye is facing a corner. Plot the y = Y0 face. 
C 
        Y = Y0+YINCB 
        X = X0+XINCB 
        IF(LENG.LE.M1*MIX+N1*NIY+5)GO TO 999 
        NS = 1 
        MS = 1 
        DO 550 MT = 2,M1 
        DO 551 MTT = 1,MSX 
        MS = MS+1 
  551   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        ZREC = ZZ 
        ZZ = 0.0 
        CALL PSP(P1,Q1) 
        IF(IFLAG1.GT.50)GO TO 5 
        CALL HDRAW(P1,Q1,1) 
        CALL HDRAW(POLD,QOLD,1) 
        CALL HDRAW(PJ,QJ,0) 
        ZZ = ZREC 
        POLD = P1 
  550   QOLD = Q1 
C 
        NSM = -20 
        NS0 = 1 
        NS1 = NSY+1 
        MS = 1 
        NS = NS1 
C 
C    Go and tidy up the rh end of the C.H.P. and 
C    then plot the surface. 
C 
        GO TO 250 
C 
C    No.  The eye is facing the x = X0 face. 
C 
C 
C    Find out which pair of grid lines the eye is looking 
C    between.  Check if the eye is exactly lined up with 
C    one of them. 
C 
   75   MS = 1 
        E = FLOAT(N1-1)*(YE1-(Y0+YINCB))/YL1+1.0 
        NSM = -10 
        NS0 = MIN0(INT(E),N1-1) 
        NS1 = NS0+1 
        IF(IFUZ(E-FLOAT(NS0)).EQ.0)NSM = NS0 
        IF(IFUZ(E-FLOAT(NS1)).NE.0) GO TO 672 
        NSM = NS1 
        IF(NS1.EQ.N1)GO TO 673 
        NS0 = NS1 
        NS1 = NS1+1 
        GO TO 672 
C 
C    Special case.  The eye is facing along an edge. 
C 
  673   NS0 = (NS0-1)*NSY+1 
        NS1 = (NS1-1)*NSY+1 
        IF(NSM.GT.0)NSM = (NSM-1)*NSY+1 
        NS = NS0 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 300 
  672   NS0 = (NS0-1)*NSY+1 
        NS1 = (NS1-1)*NSY+1 
        IF(NSM.GT.0)NSM = (NSM-1)*NSY+1 
        NS = NS1 
C 
C    Move invisibly to the first point to be plotted. 
C 
  250   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Tidy up the rh end of the C.H.P. 
C 
        L = LIST(1,LDJ) 
        LIST(1,L) = 0 
        AMAX(1,L) = PMAX+PD 
        AMAX(2,L) = QJ 
        LIST(1,LEND) = -LFREE 
        LFREE = LEND 
        LEND = L 
        LDJ = -1 
C 
C    Move one step in the +x direction. 
C 
  200   I = 1 
        J = 2 
        DO 700 MT = 1,MSX 
        MS = MS+1 
  700   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Move one step in the -y direction. 
C 
        DO 701 NT = 1,NSY 
        NS = NS-1 
  701   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Check if a change of direction is needed. 
C 
        IF(MS.EQ.MT1) GO TO 202 
        IF(NS.EQ.1) GO TO 203 
        IF(NS.EQ.NS0)GO TO 204 
        GO TO 200 
C 
C    Change direction. 
C 
  202   IF(NS.EQ.NS0)GO TO 207 
        NSA = NS0+NS1-NS 
        IF(NSA.LE.1)GO TO 205 
        NS = NSA 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 500 
  207   IF(NS.NE.1)GO TO 208 
        IF(NSM.EQ.NS.OR.NSM.EQ.-20)GO TO 205 
        DO 702 MT = 1,MSX 
        MS = MS-1 
  702   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
        GO TO 205 
  208   IF(NS.EQ.NSM)GO TO 500 
        GO TO 501 
C 
  205   NS = NS+NSY 
        IF(NS.GE.NT1)GO TO 99 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 400 
C 
  203   IF(NSM.EQ.NS.OR.NSM.EQ.-20)GO TO 209 
        DO 703 MT = 1,MSX 
        MS = MS-1 
  703   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        I = 2 
        J = 1 
        DO 704 MT = 1,MSX 
        MS = MS+1 
  704   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 400 
  209   MS = MS+MSX 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 400 
C 
  204   IF(NS.EQ.NSM) GO TO 500 
        GO TO 501 
C 
C    Move one step in the +x direction. 
C 
  300   I = 2 
        J = 1 
        DO 705 MT = 1,MSX 
        MS = MS+1 
  705   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Move one step in the +y direction. 
C 
        DO 706 NT = 1,NSY 
        NS = NS+1 
  706   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Check if a change of direction is needed. 
C 
        IF(MS.EQ.MT1)GO TO 302 
        IF(NS.EQ.NT1) GO TO 303 
        IF(NS.EQ.NS1) GO TO 401 
        IF(NS.EQ.NSM) GO TO 400 
        GO TO 300 
C 
C    Change direction. 
C 
  302   IF(NS.EQ.NS1)GO TO 307 
        IF(NS.EQ.NSM) GO TO 306 
        NSA = NS0+NS1-NS 
        IF(NSA.GE.NT1)GO TO 305 
        NS = NSA 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 400 
  307   IF(NS.NE.NT1)GO TO 401 
        IF(NSM.EQ.NS)GO TO 305 
        DO 707 MT = 1,MSX 
        MS = MS-1 
  707   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
C 
  305   NS = NS-NSY 
        IF(NS.LE.1) GO TO 99 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 500 
C 
  303   IF(NSM.EQ.NS.OR.NSM.EQ.-20) GO TO 309 
        DO 708 MT = 1,MSX 
        MS = MS-1 
  708   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        I = 1 
        J = 2 
        DO 709 MT = 1,MSX 
        MS = MS+1 
  709   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 500 
  309   MS = MS+MSX 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 500 
C 
  306   IF(NS.EQ.N1) GO TO 303 
C 
C    Move one step in the +y direction. 
C 
  400   I = 2 
        J = 1 
        DO 710 NT = 1,NSY 
        NS = NS+1 
  710   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        IF(NS.EQ.NT1)GO TO 402 
C 
C    Move one step in the -x diirection. 
C 
  401   I = 2 
        J = 1 
        DO 711 MT = 1,MSX 
        MS = MS-1 
  711   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        IF(MS.EQ.1) GO TO 403 
        GO TO 400 
C 
C    Change direction. 
C 
  402   DO 712 MT = 1,MSX 
        MS = MS-1 
  712   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        IF(MS.EQ.MT1) GO TO 405 
        GO TO 200 
C 
  403   NS = NS+NSY 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 200 
C 
  405   NS = NS0+NS1-NS 
        IF(NS.LE.1) GO TO 99 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
C 
C    Move one step in the -y direction. 
C 
  500   I = 1 
        J = 2 
        DO 713 NT = 1,NSY 
        NS = NS-1 
  713   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        IF(NS.EQ.1)GO TO 502 
C 
C    Move one step in the -x direction. 
C 
  501   I = 1 
        J = 2 
        DO 714 MT = 1,MSX 
        MS = MS-1 
  714   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        IF(MS.EQ.1) GO TO 503 
        GO TO 500 
C 
C    Change direction. 
C 
  502   IF(NSM.EQ.1.OR.NSM.EQ.-20) GO TO 300 
        DO 715 MT = 1,MSX 
        MS = MS-1 
  715   CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,1) 
        IF(IFLAG1.GT.50)GO TO 9997 
        MS = MS+MSX 
        IF(MS.EQ.MT1)GO TO 505 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 300 
C 
  503   NS = NS-NSY 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 300 
C 
  505   NS = NS0+NS1-NS 
        IF(NS.GE.NT1) GO TO 99 
        CALL SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,0) 
        IF(IFLAG1.GT.50)GO TO 9997 
        GO TO 400 
C 
C    Plot has finished successfully. 
C    Record maximum ammount of C.H.P. used. 
C 
   99   LMAX1 = LMAX 
C 
C    If the eye was facing along a grid line plot that 
C    line. 
C 
        IF(NSM.LE.0) GO TO 98 
        MS = 1 
        NS = NSM 
        X = X0+XINCB 
        Y = Y0+FLOAT(NS)*YINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 9997 
        PSTART = P 
        QSTART = Q 
        QEND = Q 
C 
C    Is x interpolation needed? 
C 
        IF(MIX2.LE.0)GO TO 30 
C 
C    Yes.  Scan along the interpolant looking for maximum Q. 
C 
        XDIR = 1.0 
        DO 32 MF = 1,MLIM 
        MS = MF+1 
        CALL XPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,X1,X2,XDIR) 
        XOFF = -XHAL+XINCS 
        X = X1+XINCS 
        DO 33 NIT = 1,MIX2 
        ZZ = POLYA*(XOFF**3)+POLYB1*XOFF*XOFF+POLYC*XOFF+POLYD 
        X = X+XINCS 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 9997 
        QEND = AMAX1(QEND,Q) 
   33   XOFF = XOFF+XINCS 
        NIT2 = MIX2+NODX 
        DO 34 NIT = 1,NIT2 
        ZZ = POLYA*(XOFF**3)+POLYB2*XOFF*XOFF+POLYC*XOFF+POLYD 
        X = X+XINCS 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 9997 
        QEND = AMAX1(QEND,Q) 
   34   XOFF = XOFF+XINCS 
   32   X = X2 
        GO TO 31 
C 
C    No.  Scan the height matrix or function looking for 
C    maximum Q. 
C 
   30   DO 97 MS = 2,MT1 
        X = X0+FLOAT(MS)*XINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)GO TO 9997 
   97   QEND = AMAX1(QEND,Q) 
   31   IF(IFUZ(QEND-QSTART).EQ.0)GO TO  98 
        CALL HDRAW(PSTART,QSTART,0) 
        CALL HDRAW(PSTART,QEND,1) 
   98   IFLAG = 0 
        CALL PLTOFF 
        RETURN 
C 
C    Out of space in the C.H.P. 
C 
  999   CALL PLTOFF 
 9998   IFLAG = 1 
        LMAX1 = LMAX 
        RETURN 
C 
C    Eye position illegal. 
C 
    5   IFLAG = 2 
        RETURN 
C 
C    Illegal value of subroutine argument. 
C 
 9999   IFLAG = 3 
        RETURN 
C 
C    Error detected after a call to PLTON. 
C 
 9997   CALL PLTOFF 
        IF(IFLAG1.LE.92)GO TO 5 
        IF(IFLAG1.EQ.999)GO TO 9998 
        IFLAG = IFLAG1 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE SEGP 2.1 DATED 20 FEBRUARY 1980 
C 
        SUBROUTINE SEGP(AMAX,LIST,LENG,Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY, 
     1    IPLT) 
C 
C 
        DIMENSION AMAX(2,LENG),LIST(2,LENG),Z(MM,NN),ZX(MMX,NNX), 
     1    ZY(MMY,NNY) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C 
C    This subroutine plots a line between two points on the grid 
C    covering the surface. 
C 
        IF(IPLT.NE.0)GO TO 1 
C 
C    The line is invisible - just move to the point. 
C 
        X = X0+FLOAT(MS)*XINCB 
        Y = Y0+FLOAT(NS)*YINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        CALL HDRAW(P,Q,0) 
        MF = MS 
        NF = NS 
        PJ = P 
        QJ = Q 
        RETURN 
C 
C    The line is visible. 
C 
    1   XDIR = FLOAT(MS-MF) 
C 
C    Moving in the x direction? 
C 
        IF(XDIR)2,50,2 
C 
C    Yes.  Is interpolation required? 
C 
    2   IF(MIX2.LE.0)GO TO 3 
C 
C    It is.  Compute the polynomial(s) to be used for 
C    interpolation. 
C 
        CALL XPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,X1,X2,XDIR) 
        X = X1+XINCS1 
        XOFF = -XHAL1+XINCS1 
C 
C    Draw the first half of the interpolated curve. 
C 
        DO 4 NIT = 1,MIX2 
        ZZ = POLYA*(XOFF**3)+POLYB1*XOFF*XOFF+POLYC*XOFF+POLYD 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        X = X+XINCS1 
    4   XOFF = XOFF+XINCS1 
        NIT2 = MIX2+NODX 
C 
C    Draw the second half of the interpolated curve. 
C 
        DO 6 NIT = 1,NIT2 
        ZZ = POLYA*(XOFF**3)+POLYB2*XOFF*XOFF+POLYC*XOFF+POLYD 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        X = X+XINCS1 
    6   XOFF = XOFF+XINCS1 
        X = X2 
        MF = MS 
        RETURN 
C 
C    Plot in the x direction without interpolation. 
C 
    3   X = X0+FLOAT(MS)*XINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        MF = MS 
        RETURN 
   50   YDIR = FLOAT(NS-NF) 
C 
C    Moving in the y direction? 
C 
        IF(YDIR)20,500,20 
C 
C    Yes.  Check if interpolation is needed. 
C 
   20   IF(NIY2.LE.0)GO TO 30 
C 
C    It is.  Compute the polynomial(s) to be used. 
C 
        CALL YPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,Y1,Y2,YDIR) 
C 
        Y = Y1+YINCS1 
        YOFF = -YHAL1+YINCS1 
C 
C    Draw the first half of the interpolated curve. 
C 
        DO 40  NIT = 1,NIY2 
        ZZ = POLYA*(YOFF**3)+POLYB1*YOFF*YOFF+POLYC*YOFF+POLYD 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        Y = Y+YINCS1 
   40   YOFF = YOFF+YINCS1 
        NIT2 = NIY2+NODY 
C 
C 
C    Draw the second half of the interpolated curve. 
C 
        DO 60 NIT = 1,NIT2 
        ZZ = POLYA*(YOFF**3)+POLYB2*YOFF*YOFF+POLYC*YOFF+POLYD 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        Y = Y+YINCS1 
   60   YOFF = YOFF+YINCS1 
        Y = Y2 
        NF = NS 
        RETURN 
C 
C    Move in the y direction without interpolation. 
C 
   30   Y = Y0+FLOAT(NS)*YINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        CALL PSP(P,Q) 
        IF(IFLAG1.GT.50)RETURN 
        PI = P 
        QI = Q 
        CALL RUN(AMAX,LIST,LENG) 
        IF(IFLAG1.GT.50)RETURN 
        PJ = PI 
        QJ = QI 
        NF = NS 
  500   RETURN 
C 
C 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE ZED9 2.3 DATED 21 MAY 1980 
C 
        SUBROUTINE ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
C 
C 
        DIMENSION Z(MM,NN),ZX(MMX,NNX),ZY(MMY,NNY) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD,NLP 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C    This subroutine rotates the entries in Z or the function ZED 
C    and calculates a height from one of them.  It also calculates 
C    the partial derivatives of the surface with respect to x and y. 
C 
C 
C    Use the entries in the array K to find the correct 
C    entry on the grid from MS and NS. 
C 
        L = (K(1)*MS+K(2)*NS+K(3))*M2+K(4)*MS+K(5)*NS+K(6) 
        MS1 = MOD(L-1,M2)+1 
        NS1 = 1+(L-1)/M2 
C 
C    User supplied functions or matrix of heights? 
C 
        IF(IFUN)GO TO 1000 
C 
C    Matrix of heights - extract a z value. 
C 
        ZZ = Z(MS1,NS1)-ZMIN 
C 
C    Use the entries in JGRAD to set up the sign 
C    and direction of the partial derivatives if 
C    needed. 
C 
        NLP = .FALSE. 
C 
C    If no x interpolation is needed consider y. 
C 
        IF(MIX2.LE.0)GO TO 20 
C 
C    Has the user supplied x gradients? 
C 
        IF(IGRADX)GO TO 1 
C 
C    Yes.  Pick the appropriate one. 
C 
        GO TO (2,3,4,5),JGRAD 
    2   GX = ZX(MS1,NS1) 
        GO TO 20 
    3   GX = -ZY(MS1,NS1) 
        GO TO 20 
    4   GX = -ZX(MS1,NS1) 
        GO TO 20 
    5   GX = ZY(MS1,NS1) 
        GO TO 20 
C 
C    No.  Fit parabolas in the x and y directions. 
C 
    1   IF(MS1.EQ.1)GO TO 6 
        IF(MS1.EQ.MTOP)GO TO 7 
        MS2 = MS1-1 
        MS3 = MS1+1 
        GXX = (Z(MS3,NS1)-Z(MS2,NS1))*0.5/XINCB 
        GO TO 8 
    6   GXX = (Z(2,NS1)-Z(1,NS1))/XINCB 
        GO TO 8 
    7   GXX = (Z(MTOP,NS1)-Z(MTOP-1,NS1))/XINCB 
    8   IF(NS1.EQ.1)GO TO 9 
        IF(NS1.EQ.NTOP)GO TO 10 
        NS2 = NS1-1 
        NS3 = NS1+1 
        GY = (Z(MS1,NS3)-Z(MS1,NS2))*0.5/YINCB 
        GO TO 11 
    9   GY = (Z(MS1,2)-Z(MS1,1))/YINCB 
        GO TO 11 
   10   GY = (Z(MS1,NTOP)-Z(MS1,NTOP-1))/YINCB 
   11   IF(NLP)GO TO 21 
C 
C    Pick the appropriate parabola gradient. 
C 
        GO TO (12,13,14,15),JGRAD 
   12   GX = GXX 
        GO TO 20 
   13   GX = -GY 
        GO TO 20 
   14   GX = -GXX 
        GO TO 20 
   15   GX = GY 
C 
C    If no y interpolation is needed return. 
C 
   20   IF(NIY2.LE.0)RETURN 
        NLP = .TRUE. 
C 
C    Has the user supplied y gradients? 
C 
        IF(IGRADY)GO TO 22 
C 
C    Yes.  Pick the appropriate one. 
C 
        GO TO (32,33,34,35),JGRAD 
   32   GY = ZY(MS1,NS1) 
        RETURN 
   33   GY = ZX(MS1,NS1) 
        RETURN 
   34   GY = -ZY(MS1,NS1) 
        RETURN 
   35   GY = -ZX(MS1,NS1) 
   36   RETURN 
C 
C    No.  If parabolas have already been fitted use them. 
C    Otherwise go back and fit them. 
C 
   22   IF(MIX2.LE.0)GO TO 1 
C 
C    Pick the appropriate parabola gradient. 
C 
   21   GO TO (36,37,38,39),JGRAD 
   37   GY = GXX 
        RETURN 
   38   GY = -GY 
        RETURN 
   39   GY = -GXX 
        RETURN 
C 
C    User supplied functions - extract a z value. 
C 
 1000   ZZ = ZED(MS1,NS1)-ZMIN 
C 
C    Use the entries in JGRAD to set up the sign 
C    and direction of the partial derivatives if 
C    needed. 
C 
        NLP = .FALSE. 
C 
C    If no x interpolation is needed consider y. 
C 
        IF(MIX2.LE.0)GO TO 120 
C 
C    Has the user supplied x gradients? 
C 
        IF(IGRADX)GO TO 101 
C 
C    Yes.  Pick the appropriate one. 
C 
        GO TO (102,103,104,105),JGRAD 
  102   GX = ZEDX(MS1,NS1) 
        GO TO 120 
  103   GX = -ZEDY(MS1,NS1) 
        GO TO 120 
  104   GX = -ZEDX(MS1,NS1) 
        GO TO 120 
  105   GX = ZEDY(MS1,NS1) 
        GO TO 120 
C 
C    No.  Fit parabolas in the x and y directions. 
C 
  101   IF(MS1.EQ.1)GO TO 106 
        IF(MS1.EQ.MTOP)GO TO 107 
        MS2 = MS1-1 
        MS3 = MS1+1 
        GXX = (ZED(MS3,NS1)-ZED(MS2,NS1))*0.5/XINCB 
        GO TO 108 
  106   GXX = (ZED(2,NS1)-ZED(1,NS1))/XINCB 
        GO TO 108 
  107   GXX = (ZED(MTOP,NS1)-ZED(MTOP-1,NS1))/XINCB 
  108   IF(NS1.EQ.1)GO TO 109 
        IF(NS1.EQ.NTOP)GO TO 110 
        NS2 = NS1-1 
        NS3 = NS1+1 
        GY = (ZED(MS1,NS3)-ZED(MS1,NS2))*0.5/YINCB 
        GO TO 111 
  109   GY = (ZED(MS1,2)-ZED(MS1,1))/YINCB 
        GO TO 111 
  110   GY = (ZED(MS1,NTOP)-ZED(MS1,NTOP-1))/YINCB 
  111   IF(NLP)GO TO 121 
C 
C    Pick the appropriate parabola gradient. 
C 
        GO TO (112,113,114,115),JGRAD 
  112   GX = GXX 
        GO TO 120 
  113   GX = -GY 
        GO TO 120 
  114   GX = -GXX 
        GO TO 120 
  115   GX = GY 
C 
C    If no y interpolation is needed return. 
C 
  120   IF(NIY2.LE.0)RETURN 
        NLP = .TRUE. 
C 
C    Has the user supplied y gradients? 
C 
        IF(IGRADY)GO TO 122 
C 
C    Yes.  Pick the appropriate one. 
C 
        GO TO (132,133,134,135),JGRAD 
  132   GY = ZEDY(MS1,NS1) 
        RETURN 
  133   GY = ZEDX(MS1,NS1) 
        RETURN 
  134   GY = -ZEDY(MS1,NS1) 
        RETURN 
  135   GY = -ZEDX(MS1,NS1) 
  136   RETURN 
C 
C    No.  If parabolas have already been fitted use them. 
C    Otherwise go back and fit them. 
C 
  122   IF(MIX2.LE.0)GO TO 101 
C 
C    Pick the appropriate parabola gradient. 
C 
  121   GO TO (136,137,138,139),JGRAD 
  137   GY = GXX 
        RETURN 
  138   GY = -GY 
        RETURN 
  139   GY = -GXX 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE PSP 2.1 DATED 9 DECEMBER 1978 
C 
        SUBROUTINE PSP(P,Q) 
C 
C 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C    This subroutine transforms a point (X,Y,ZZ) in 
C    three dimensions to a point (P,Q) in the image 
C    plane.  (P,Q) is the true perspective image of 
C    the point (X,Y,ZZ). 
C 
        R = W1*X+W2*Y+W3*ZZ-WE 
        IF(R.GT.0.0) GO TO 1 
        P = 0.0 
        Q = 0.0 
        IFLAG1 = 90 
        RETURN 
    1   P = (U1*X+U2*Y+U3*ZZ-UE)/R 
        Q = (V1*X+V2*Y+V3*ZZ-VE)/R 
        IFLAG1 = 0 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE PSPSET 2.1 DATED 9 DECEMBER 1978 
C 
        SUBROUTINE PSPSET(E1,E2,E3,C1,C2,C3,T1,T2,T3) 
C 
C 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C    This subroutine sets up the perspective transformation 
C    used by the routine PSP. 
C 
C 
C 
C    E gives the eye position 
C    C gives the position of a point required to be in the 
C      centre of the picture 
C    T gives the position of a point required to be above 
C      the centre of the picture 
C 
        W1 = C1-E1 
        W2 = C2-E2 
        W3 = C3-E3 
        R = SQRT(W1*W1+W2*W2+W3*W3) 
        IF(R.GT.0.0)GO TO 1 
        IFLAG1 = 91 
        RETURN 
    1   W1 = W1/R 
        W2 = W2/R 
        W3 = W3/R 
        WE = W1*E1+W2*E2+W3*E3 
        U1 = W2*(T3-E3)-W3*(T2-E2) 
        U2 = W3*(T1-E1)-W1*(T3-E3) 
        U3 = W1*(T2-E2)-W2*(T1-E1) 
        R = SQRT(U1*U1+U2*U2+U3*U3) 
        IF(R.GT.0.0)GO TO 2 
        IFLAG1 = 92 
        RETURN 
    2   U1 = U1/R 
        U2 = U2/R 
        U3 = U3/R 
        UE = U1*E1+U2*E2+U3*E3 
        V1 = U2*W3-U3*W2 
        V2 = U3*W1-U1*W3 
        V3 = U1*W2-U2*W1 
        VE = V1*E1+V2*E2+V3*E3 
        IFLAG1 = 0 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE RUN 2.2 DATED 18 NOVEMBER 1979 
C 
        SUBROUTINE RUN(AMAX,LIST,LENG) 
C 
C 
        DIMENSION AMAX(2,LENG),LIST(2,LENG) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD,LAST 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C 
C     This subroutine plots a line from PJ,QJ (where it assumes the 
C     pen to be) to PI,QI.  It hides any of that line that lies 
C     below the C.H.P. 
C 
        LAST = .FALSE. 
C 
C    Compute equation of line to be plotted: 
C       A*P + B = Q 
C 
        A = PI-PJ 
        IF(IFUZ(A*10.0))2,1000,2 
    2   A = (QI-QJ)/A 
        B = QI-A*PI 
C 
C    Find segment of C.H.P. in which PI,QI lies 
C 
        P = PI 
        CALL SEARCH(AMAX,LIST,LENG,P,LI2,LJ2) 
        IF(IFLAG1.GT.50)RETURN 
        LJ0 = LIST(J,LI2) 
        LGUESS = LI2 
        IF(IFLAG1.GT.2)GO TO 1001 
C 
C    If C.H.P. contains a discontinuity plotting 
C    must start from it. 
C 
        IF(LDJ)100,100,3 
C 
C    The segment of the C.H.P. under consideration is LJJ - LII. 
C    Do not bother with the vertical part of the discontinuity. 
C 
    3   LJJ = LIST(I,LDJ) 
        LJ1 = LDJ 
        LI1 = LDJ 
        LI0 = LIST(I,LDJ) 
    4   LII = LIST(I,LJJ) 
C 
C    Does the new line cross LJJ-LII of the C.H.P.? 
C 
        CALL ACROSS(AMAX,LENG,LII,LJJ,PX,QX,IDG) 
C 
C    Is LJJ the start of the last C.H.P. segment to be considered? 
C 
        IF(LJJ.EQ.LJ0)GO TO 50 
C 
C    Test for crossing. 
C 
        IF(IFLAG1-3)6,6,5 
C 
C    None found.  Go on to the next segment of the C.H.P. 
C 
    5   LJJ = LII 
        GO TO 4 
C 
C    Crossing found.  Draw a vector to the point of intersection. 
C 
C 
C    Check if intersection really is in the new line. 
C 
    6   IF(IFUZ(FLOAT(J-I)*(PI-PX)))400,81,81 
   81   CALL HDRAW(PX,QX,1) 
C 
C    Insert the intersection in the C.H.P. 
C 
        LL = LIST(I,LDJ) 
        IF(IFLAG1-2)9,8,7 
C 
C    Intersection coincides with LJJ. 
C 
    7   LIST(I,LDJ) = LJJ 
        LIST(J,LJJ) = LDJ 
        LE = LJJ 
        LJ1 = LJJ 
C 
C    Wipe out the C.H.P. between LL and LE, returning entries found to 
C    the free chain. 
C 
   10   IF(LL.EQ.LE)GO TO 11 
        LX = LIST(I,LL) 
        LIST(1,LL) = -LFREE 
        LFREE = LL 
        LIN = LIN-1 
        LL = LX 
        GO TO 10 
   11   LDJ = -1 
        IFLAG1 = 0 
C 
C    C.H.P. has been changed and a new line plotted to a intersection. 
C    Is any of the new line left? 
C 
        IF(LAST)RETURN 
C 
C    Some new line left.  Set PJ and QJ to the intersection point 
C    and jump to the part of the routine dealing with running under 
C    the C.H.P. 
C 
        LI1 = LJ1 
        LI0 = LIST(I,LJ1) 
        PJ = PX 
        QJ = QX 
        GO TO 101 
C 
C    Intersection coincides with LII. 
C 
    8   LIST(I,LDJ) = LII 
        LIST(J,LII) = LDJ 
        LE = LII 
        LJ1 = LII 
C 
C    Wipe out dead section of the C.H.P. 
C 
        GO TO 10 
C 
C    Intersection lies between LJJ and LII. Make an entry in the C.H.P. 
C    at it. 
C 
    9   LIN = LIN+1 
        LMAX = MAX0(LIN,LMAX) 
        IF(LIN-LENG)12,12,13 
C 
C    C.H.P. overflow. 
C 
   13   IFLAG1 = 999 
        RETURN 
   12   LF = LFREE 
        LFREE = -LIST(1,LFREE) 
        LIST(I,LDJ) = LF 
        LIST(J,LF) = LDJ 
        LIST(I,LF) = LII 
        LIST(J,LII) = LF 
        AMAX(1,LF) = PX 
        AMAX(2,LF) = QX 
        LE = LII 
        LJ1 = LF 
C 
C    Wipe out dead section of the C.H.P. 
C 
        GO TO 10 
C 
C    In the segment of the C.H.P. where the end of the 
C    new line lies.  Set last .TRUE. 
C 
   50   LAST = .TRUE. 
        IF(IFLAG1.GE.4)GO TO 400 
        GO TO 6 
C 
C    No discontinuity.  Find the segment of the C.H.P. where the new 
C    line starts. 
C 
  100   P = PJ 
        CALL SEARCH(AMAX,LIST,LENG,P,LI1,LJ1) 
        IF(IFLAG1.GT.50)RETURN 
        IF(IFLAG1.GT.2)GO TO 1001 
C 
C    LI0 is the I end of the segment of the C.H.P. in which 
C    PJ lies (LJ1-LI0). 
C 
        LI0 = LIST(I,LJ1) 
        P = PJ 
        Q = QJ 
C 
C    Find out if PJ,QJ is (above), below or coiincides 
C    with LJ1-LI0. 
C 
        CALL ABOVE(AMAX,LENG,LI0,LJ1,P,Q) 
        GO TO (759,759,102,1003,1004),IFLAG1 
C 
C    Start of the new line is below the C.H.P.  Start tracking 
C    along the C.H.P.... 
C 
  759   LJJ = LJ1 
  754   LII = LIST(I,LJJ) 
C 
C    ..Checking if the new line crosses it. 
C 
        CALL ACROSS(AMAX,LENG,LII,LJJ,PX,QX,IDG) 
C 
C    If crossing occurs in the first segment of the C.H.P. 
C    check that it does so strictly to the I side of PJ. 
C 
        IF(LJJ.NE.LJ1.OR.IFLAG1.GT.3)GO TO 772 
        IF(IFUZ(FLOAT(J-I)*(PX-PJ)))773,773,772 
  773   IFLAG1 = 5 
  772   IF(LJJ.EQ.LJ0)GO TO 350 
        IF(IFLAG1-3)756,756,755 
  755   LJJ = LII 
        GO TO 754 
C 
C    Start of new line coincides with C.H.P.  Find its Q 
C    value at LI0. 
C 
  102   Q = AMAX(1,LI0)*A+B 
C 
C    Is that Q below, coincident or above the C.H.P. at LI0? 
C 
        IF(IFUZ(Q-AMAX(2,LI0)))101,101,103 
C 
C    The new line starts off visible.  Track along the C.H.P... 
C 
  103   LJJ = LJ1 
C 
C    Check if the new line is shorter than the first seg of the C.H.P. 
C 
        IF(LJJ.EQ.LJ0)GO TO 150 
        LJJ = LIST(I,LJJ) 
  104   LII = LIST(I,LJJ) 
C 
C    .. Looking for crossings. 
C 
        CALL ACROSS(AMAX,LENG,LII,LJJ,PX,QX,IDG) 
        IF(LJJ.EQ.LJ0)GO TO 150 
C 
C    Check for intersection. 
C 
        IF(IFLAG1-3)106,106,105 
C 
C    None found.  Consider next segment. 
C 
  105   LJJ = LII 
        GO TO 104 
C 
C    Intersection found. 
C    Plot to it. 
C 
C 
C    Check if intersection really is in the new line. 
C 
  106   IF(IFUZ(FLOAT(J-I)*(PI-PX)))400,181,181 
  181   CALL HDRAW(PX,QX,1) 
C 
C    Insert the intersection in the C.H.P. 
C 
        LL = LIST(I,LJ1) 
        IF(IFLAG1-2)109,108,107 
C 
C    Intersection coincides with LJJ. 
C 
  107   LIST(I,LJ1) = LJJ 
        LIST(J,LJJ) = LJ1 
        LE = LJJ 
        LJ1 = LJJ 
C 
C    Wipe out C.H.P. between LL and LE, returning any 
C    entries found to the free chain. 
C 
  110   IF(LL.EQ.LE)GO TO 111 
        LX = LIST(I,LL) 
        LIST(1,LL) = -LFREE 
        LFREE = LL 
        LIN = LIN-1 
        LL = LX 
        GO TO 110 
C 
C    C.H.P. has been modified.  Is any of the new line left? 
C 
  111   IFLAG1 = 0 
        IF(LAST)RETURN 
        PJ = PX 
        QJ = QX 
        LI1 = LJ1 
        LI0 = LIST(I,LJ1) 
        GO TO 101 
C 
C    Intersection coincides with LII. 
C 
  108   LIST(I,LJ1) = LII 
        LIST(J,LII) = LJ1 
        LE = LII 
        LJ1 = LII 
        GO TO 110 
C 
C    Intersection lies between LJJ and LII.  Make an 
C    extra entry in the C.H.P. at it. 
C 
  109   LIN = LIN+1 
        LMAX = MAX0(LIN,LMAX) 
        IF(LIN-LENG)112,112,13 
  112   LF = LFREE 
        LFREE = -LIST(1,LF) 
        LIST(I,LJ1) = LF 
        LIST(J,LF) = LJ1 
        LIST(I,LF) = LII 
        LIST(J,LII) = LF 
        AMAX(1,LF) = PX 
        AMAX(2,LF) = QX 
        LE = LII 
        LJ1 = LF 
C 
C    Wipe out dead section of the C.H.P. 
C 
        GO TO 110 
  150   LAST = .TRUE. 
C 
C    If a discontinuity has been generated go to 400. 
C 
        IF(IFLAG1.GE.4.OR.LJ1.EQ.LJ0)GO TO 400 
        GO TO 106 
C 
C    Discontinuity needs to be drawn and inserted in the C.H.P. 
C 
  400   CALL HDRAW(PI,QI,1) 
C 
C    Start C.H.P. modification. 
C 
        LL = LIST(I,LJ1) 
C 
C    Check if there is enough free space in the C.H.P. 
C 
        IF(LIN+2-LENG)401,401,13 
C 
C    LJ1 is the start of the new line. 
C 
  401   LF1 = LFREE 
        LFREE = -LIST(1,LF1) 
        LIST(I,LJ1) = LF1 
        LIST(J,LF1) = LJ1 
        AMAX(1,LF1) = PI 
        AMAX(2,LF1) = QI 
        LDJ = LF1 
        LIN = LIN+1 
        LMAX = MAX0(LIN,LMAX) 
C 
C    Drop a perpendicular from the discontinuity. 
C 
C 
C    Check if the bottom of the discontinuity coincides 
C    with a point in the C.H.P. 
C 
        IF(LJ2.EQ.LI2)GO TO 402 
C 
C    It does not. 
C 
        LF = LFREE 
        LFREE = -LIST(1,LF) 
        LIST(I,LF1) = LF 
        LIST(J,LF) = LF1 
        AMAX(1,LF) = PI 
        AMAX(2,LF) = A1*PI+B1 
        LIN = LIN+1 
        LMAX = MAX0(LIN,LMAX) 
        LIST(I,LF) = LI2 
        LIST(J,LI2) = LF 
        LE = LI2 
        GO TO 110 
C 
C    It does. 
C 
  402   LE = LJ2 
        LIST(I,LF1) = LE 
        LIST(J,LE) = LF1 
        GO TO 110 
C 
C    The new line starts off invisibly.  Track along the C.H.P..... 
C 
  101   LJJ = LJ1 
C 
C    Check if new line is shorter than first segment 
C    of the C.H.P. 
C 
        IF(LJJ.EQ.LJ0)GO TO 389 
        LJJ = LIST(I,LJJ) 
  302   LII = LIST(I,LJJ) 
C 
C    .. Looking for intersections. 
C 
        CALL ACROSS(AMAX,LENG,LII,LJJ,PX,QX,IDG) 
C 
C    Check for end of new line. 
C 
        IF(LJJ.EQ.LJ0)GO TO 350 
C 
C    Has an intersection occured? 
C 
        IF(IFLAG1-3)756,756,305 
C 
C    No - look at the next segment of the C.H.P. 
C 
  305   LJJ = LII 
        GO TO 302 
C 
C    Crossing has occured.  Move to it and insert it in the C.H.P. 
C 
  756   CALL HDRAW(PX,QX,0) 
        PJ = PX 
        QJ = QX 
        IF(IFLAG1-2)736,737,738 
C 
C    Intersection lies between LJJ and LII. 
C 
  736   LIN = LIN+1 
        LMAX = MAX0(LMAX,LIN) 
        IF(LIN-LENG)735,735,13 
  735   LF = LFREE 
        LFREE = -LIST(1,LF) 
        LIST(I,LJJ) = LF 
        LIST(J,LF) = LJJ 
        LIST(I,LF) = LII 
        LIST(J,LII) = LF 
        AMAX(1,LF) = PX 
        AMAX(2,LF) = QX 
        LJJ = LF 
        GO TO 738 
C 
C    Intersection coincides with LII. 
C 
  737   LJJ = LII 
C 
C    Intersection coincides with LJJ. 
C 
  738   LJ1 = LJJ 
        LI1 = LJJ 
        LI0 = LIST(I,LJJ) 
C 
C    Is part of the new line left? 
C 
        P = PI 
        Q = QI 
        CALL ABOVE(AMAX,LENG,LI0,LJ1,P,Q) 
        IF(IFLAG1.EQ.1.OR.IFLAG1.EQ.5)GO TO 103 
C 
C    No - return 
C 
        IFLAG1 = 0 
        RETURN 
C 
C    The end of the new line.  Does it cross the 
C    last segment? 
C 
  350   IF(IFLAG1.LT.4)GO TO 388 
C 
C    No - return 
C 
  389   IFLAG1 = 0 
        RETURN 
C 
C    Yes - it crosses. 
C 
  388   P = PI 
        Q = QI 
        IFL1 = IFLAG1 
C 
C    Check if it really does. 
C 
        CALL ABOVE(AMAX,LENG,LII,LJJ,P,Q) 
        IF(IFLAG1.EQ.2)GO TO 389 
C 
C    Yes - move to point of intersection. 
C 
        CALL HDRAW(PX,QX,0) 
C 
C    Insert that point in the C.H.P. 
C 
        IF(IFL1-2)360,389,362 
  362   IF(IDG.EQ.2)GO TO 389 
C 
C    A bit of the line is visible. Create the appropriate discontinuity. 
C 
        CALL HDRAW(PI,QI,1) 
        LIN = LIN+1 
        LMAX = MAX0(LIN,LMAX) 
        IF(LIN-LENG)700,700,13 
  700   LF = LFREE 
        LFREE = -LIST(1,LF) 
        LIST(I,LJJ) = LF 
        LIST(J,LF) = LJJ 
        AMAX(1,LF) = PI 
        AMAX(2,LF) = QI 
        LDJ = LF 
C 
C    Drop a perpendicular. 
C 
        IF(LJ2.EQ.LI2)GO TO 749 
        LIN = LIN+1 
        LMAX = MAX0(LMAX,LIN) 
        IF(LIN-LENG)701,701,13 
  701   LF1 = LFREE 
        LFREE = -LIST(1,LF1) 
        LIST(I,LF) = LF1 
        LIST(J,LF1) = LF 
        LIST(I,LF1) = LII 
        LIST(J,LII) = LF1 
        AMAX(1,LF1) = PI 
        AMAX(2,LF1) = A1*PI+B1 
        IFLAG1 = 0 
        RETURN 
C 
C    Perpendicular coincides with LII. 
C 
  749   LIST(I,LF) = LII 
        LIST(J,LII) = LF 
        IFLAG1 = 0 
        RETURN 
C 
C    Put the intersection in the C.H.P. 
C 
  360   LIN = LIN+1 
        LMAX = MAX0(LMAX,LIN) 
        IF(LIN-LENG)702,702,13 
  702   LF1 = LFREE 
        LFREE = -LIST(1,LF1) 
        AMAX(1,LF1) = PX 
        AMAX(2,LF1) = QX 
        LIST(I,LJJ) = LF1 
        LIST(J,LF1) = LJJ 
        LIST(I,LF1) = LII 
        LIST(J,LII) = LF1 
        LJJ = LF1 
        GO TO 362 
C 
C    Error returns. 
C 
 1000   IFLAG1 = 97 
        RETURN 
 1001   IFLAG1 = 96 
        RETURN 
 1002   IFLAG1 = 95 
        RETURN 
 1003   IFLAG1 = 94 
        RETURN 
 1004   IFLAG1 = 93 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE YPOL 2.1 DATED 29 APRIL 1980 
C 
        SUBROUTINE YPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,Y1,Y2,YDIR) 
C 
C 
        DIMENSION Z(MM,NN),ZX(MMX,NNX),ZY(MMY,NNY) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C    This subroutine computes polynomial(s) to interpolate 
C    between two adjacent points in the y direction. 
C 
C 
C    Compute the value and gradient at the first point. 
C 
C 
        G1 = GY 
        Z1 = ZZ 
        Y1 = Y 
C 
C    Compute the value and gradient at the second point. 
C 
        Y2 = Y0+FLOAT(NS)*YINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        G2 = GY 
        Z2 = ZZ 
        YINCS1 = YINCS*YDIR 
        YHAL1 = YDIR*YHAL 
        RYHAL1 = RYHAL*YDIR 
C 
C    Is a seamed quadratic required? 
C 
        IF(IQUAD)GO TO 1 
C 
C    No.  Fit a cubic. 
C 
        POLYA = 0.25*RYHAL1*RYHAL1*(G1+G2-RYHAL1*(Z2-Z1)) 
        POLYB1 = 0.25*RYHAL1*(G2-G1) 
        POLYB2 = POLYB1 
        POLYC = 0.25*(3.0*RYHAL1*(Z2-Z1)-G1-G2) 
        POLYD = 0.5*(Z1+Z2-0.5*YHAL1*(G2-G1)) 
        RETURN 
C 
C    Yes.  Fit two quadratics. 
C 
    1   POLYA = 0.0 
        POLYB1 = 0.5*RYHAL1*(RYHAL1*(Z2-Z1)-1.5*G1-0.5*G2) 
        POLYB2 = 0.5*RYHAL1*(RYHAL1*(Z1-Z2)+0.5*G1+1.5*G2) 
        POLYC = RYHAL1*(Z2-Z1)-0.5*(G1+G2) 
        POLYD = 0.5*(Z1+Z2)+0.25*YHAL1*(G1-G2) 
C 
        RETURN 
        END 
C 
C 
C 
C*********************************************************************** 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE XPOL 2.1  DATED 29 APRIL 1980 
C 
        SUBROUTINE XPOL(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY,X1,X2,XDIR) 
C 
C 
        DIMENSION Z(MM,NN),ZX(MMX,NNX),ZY(MMY,NNY) 
        LOGICAL IFUN,IGRADX,IGRADY,IHAZE,IQUAD 
C 
        COMMON/SOLBLK/U1,U2,U3,UE,V1,V2,V3,VE,W1,W2,W3,WE, 
     1    FUZ,D0,DD,X,Y,ZZ,XINCB,YINCB,K(6),IFUN,IFLAG1,IHAZE,ZMIN, 
     2    I,J,LFREE,X0,Y0,LIN,LMAX,LGUESS,LDJ,LSTART,LEND, 
     3    XE1,YE1,ZE1,PI,QI,PJ,QJ,MS,NS,M2,A,B,A1,B1,NODX,NODY,JGRAD, 
     4    XINCS,YINCS,GX,GY,XHAL,YHAL,RXHAL,RYHAL,MF,NF,MIX2,NIY2, 
     5    POLYA,POLYB1,POLYB2,POLYC,POLYD,IQUAD,XINCS1,XHAL1,YINCS1, 
     6    YHAL1,MTOP,NTOP,IGRADX,IGRADY,MSX,NSY 
C 
C    This subroutine computes polynomial(s) to interpolate 
C    between two adjacent grid points in the x direction. 
C 
C 
C    Compute the values and gradients at the first point. 
C 
        G1 = GX 
        Z1 = ZZ 
        X1 = X 
C 
C    Compute the values and gradients at the second point. 
C 
        X2 = X0+FLOAT(MS)*XINCB 
        CALL ZED9(Z,ZX,ZY,MM,NN,MMX,NNX,MMY,NNY) 
        G2 = GX 
        Z2 = ZZ 
        XINCS1 = XINCS*XDIR 
        XHAL1 = XDIR*XHAL 
        RXHAL1 = RXHAL*XDIR 
C 
C    Is a seamed quadratic required? 
C 
        IF(IQUAD)GO TO 1 
C 
C    No.  Fit a cubic. 
C 
        POLYA = 0.25*RXHAL1*RXHAL1*(G1+G2-RXHAL1*(Z2-Z1)) 
        POLYB1 = 0.25*RXHAL1*(G2-G1) 
        POLYB2 = POLYB1 
        POLYC = 0.25*(3.0*RXHAL1*(Z2-Z1)-G1-G2) 
        POLYD = 0.5*(Z1+Z2-0.5*XHAL1*(G2-G1)) 
        RETURN 
C 
C    Yes.  Fit two quadratics. 
C 
    1   POLYA = 0.0 
        POLYB1 = 0.5*RXHAL1*(RXHAL1*(Z2-Z1)-1.5*G1-0.5*G2) 
        POLYB2 = 0.5*RXHAL1*(RXHAL1*(Z1-Z2)+0.5*G1+1.5*G2) 
        POLYC = RXHAL1*(Z2-Z1)-0.5*(G1+G2) 
        POLYD = 0.5*(Z1+Z2)+0.25*XHAL1*(G1-G2) 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN SOLID 2  CREATED 9 DECEMBER 1978 
C    SUBROUTINE HDRAW 2.1  DATED 2 FEBRUARY 1981 
C 
        SUBROUTINE HDRAW(X,Y,IPLT) 
C 
C    This subroutine interfaces SOLID graphics to the graphics 
C    specification for TILE 4.  If the TILE 4 graphics 
C    interface has been written no further work is needed to 
C    make SOLID 2 run. 
C 
        IF(IPLT.NE.0) GO TO 1 
        CALL PLTMOV(X,Y) 
        RETURN 
    1   CALL PLTLIN(X,Y) 
        RETURN 
        END 
C 
C*********************************************************************** 
