C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    PROGRAM TESTC 4.3 DATED 26 MARCH 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C 
C    This program tests the TILE 4 routines 
C        SCATRC 
C        TILE4M 
C        GRLD 
C        RECT1G 
C    and the routines that they call. 
C 
C    It constructs a tessellation of a random point pattern 
C    on a rectangle, establishes the values of some test 
C    functions at the data points and then reconstructs the 
C    functions over a smaller rectangle using the C1 natural 
C    neighbour interpolant.  The values of the original 
C    functions and their reconstructed values are written 
C    to a file without carriage control characters for 
C    later plotting.  This can be supressed by setting the variable 
C    NSTORE zero or negative.  Statistics on the accuracy of the 
C    reconstructions are written to another file with carriage 
C    control characters. 
C 
	SUBROUTINE NNINT(NBP,XX,YY,VAL,XMIN2,YMIN2,XMAX2,YMAX2,
	1  NCELL,XCELL,YCELL,VCELL) !,NearestWell)
C2      INTEGER*2 L 
        EXTERNAL UNI4M 
        DIMENSION CN(3,4),JADDR(4),PT(2,NBP),NADDR(NBP),L(10000000), 
     1    VAL(NBP),SBAREA(50),PTOFF(2,50),VALOFF(50),GRAD(2,NBP), 
     2    DELSBA(2,50),Z(25,25),ZX(25,25),ZY(25,25) 
     
     
      Integer NCELL
      DIMENSION VSHAPE(NBP),VV(NBP)
	DIMENSION XX(NBP),YY(NBP)	
	DIMENSION XCELL(NCELL),YCELL(NCELL),VCELL(NCELL)
	Dimension NearestWell(NCELL)
C 
        DATA JCNS,NPTS,LTOP,KNBMAX,MM,NN/4,5,10000000,50,20,20/ 
!        DATA XMIN2,XMAX2,YMIN2,YMAX2/248.99,299.01,248.99,299.01/ 
        DATA MSEED/97129/ 
        DATA NFMAX/8/ 
        DATA EPSCN,EPSPT/0.0,0.0/ 
        DATA NSTART/0/ 
        DATA NWRITE,NSTORE/12,13/ 

C DEFINE THE RANGE OF THE RECTANGLE DOMAIN
!	XMIN2= 402584.100928644
!	XMAX2= 439918.002897194
!	YMIN2= 3562000.3572702
!	YMAX2= 3654114.0285765
		
 !X   000000       000000     
 !Y   00000        00000  

C STARING TILE4	  	  
	  NPTS=1
        MSD = MSEED 
        CALL SCATRC(CN,JCNS,PT,NPTS,IFLAG,UNI4M,MSD,XMIN2,XMAX2, 
     1    YMIN2,YMAX2) 
	  NPTS=NBP
        IF(IFLAG.EQ.0) GOTO 1 
        WRITE(NWRITE,9000)IFLAG 
        STOP 612 
	  
!
1	CONTINUE

	DO I=1,NBP
	PT(1,I)=XX(I) !PPT(1,I)
	PT(2,I)=YY(I) !PPT(2,I)
	ENDDO

       CALL TILE4M(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB) 
        IF(IFLAG.EQ.0) GOTO 2 
        WRITE(NWRITE,9000)IFLAG 
!        STOP 613 
C 
C    Loop through the NFMAX functions. 
C 
2      Continue    
!    2   DO 7  NFUNC = 1,NFMAX 
C 
C    Set the function values at the data sites in VAL. 
C 
!        DO 3  NPT = 1,NPTS 
!    3   VAL(NPT) = FUNC(PT(1,NPT),PT(2,NPT),NFUNC) 
!	DO I=1,NBP
!	VAL(I)=VV(I)
!	ENDDO

C 
C    Call GRLD to estimate the gradients at the data sites. 
C 
!        CALL GRLD(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VAL,SBAREA, 
!     1    PTOFF,VALOFF,KNBMAX,GRA`D) 
!        IF(IFLAG.EQ.0) GOTO 4 
!        WRITE(NWRITE,9000)IFLAG 
!        STOP 614 
C 
C    Call RECT1G to compute the C1 interpolant for the function 
C    at the nodes of the MM*NN grid on the smaller rectangle. 
C 
	VSHAPE=1.0
!	OPEN(UNIT=16, FILE='ZONE.DATA')
!	WRITE(16,*) '2D VORONOI ZONES: I, J, NEAREST POINT' 
	DO I=1,NCELL
	XPT=XCELL(I)
	YPT=YCELL(I)

    4   CALL NNBR0(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE,VAL,Z0
     3   ,VSHAPE)

	NearestWell(I)=NINDEX
	VCELL(I)= Z0
!	WRITE(16,*) I,J,NINDEX
	IF (NINDEX.LT.0) THEN
	WRITE(*,*)
	ENDIF

	ENDDO	
	
        IF(IFLAG.EQ.0) GOTO 5 
        WRITE(NWRITE,9000)IFLAG 
        write(NWRITE,*) XPT,YPT
        STOP 615 
C 
C    Record the original function and its reconstruction on NSTORE 
C    if required and compute error statistics. 
C 

C 
C    Run complete. 
C 
  5       Continue
! 5       WRITE(NWRITE,9002) 
!        STOP 
	RETURN
 8000   FORMAT(3(I13,2X)) 
 8001   FORMAT(4(E13.6,2X)) 
 9000   FORMAT(///1H ,19(1H*)/1H ,16H ERROR: IFLAG = ,I3/1H ,19(1H*)///) 
 9001   FORMAT(///1H ,10HFUNCTION: ,I3/1H ,16HRMS ERROR:      ,E13.6/ 
     1    1H ,16HABSOLUTE ERROR: ,E13.6/1H ,16HMAXIMUM ERROR:  ,E13.6 
     2    ///) 
! 9002   FORMAT(////1H ,22HTESTC 4.3 RUN COMPLETE//) 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    FUNCTION FUNC 4.1 DATED 8 JANUARY 1981 
C 
        FUNCTION FUNC(X,Y,NFUNC) 
C 
C    This routine supplies the value of one of eight 
C    test functions at the point (X,Y). 
C 
        IF(NFUNC.LE.0.OR.NFUNC.GT.8) STOP 616 
        GOTO (1,2,3,4,5,6,7,8),NFUNC 
C 
C    Function 1: Cliff edge. 
C 
    1   FUNC = 0.0 
        IF(X.GT.0.5) FUNC = 1.0 
        RETURN 
C 
C    Function 2: V-shaped valley. 
C 
    2   FUNC = ABS(X-Y) 
        RETURN 
C 
C    Function 3: Spherical quadratic (parabolic function). 
C 
    3   FUNC = (X-0.25)**2+(Y-0.25)**2 
        RETURN 
C 
C    Function 4: Non-spherical quadratic (hyperbolic paraboloid). 
C 
    4   FUNC = (X-0.5)**2-(Y-0.5)**2 
        RETURN 
C 
C    Function 5: Spherical bivariate normal density function. 
C 
    5   FUNC = EXP(-8.0*((X-0.75)**2+(Y-0.75)**2)) 
        RETURN 
C 
C    Function 6: Sum of two normal densities to give a bimodal 
C                density function. 
C 
    6   FUNC = EXP(-16.0*((X-0.75)**2+(Y-0.25)**2)) 
        FUNC = FUNC+EXP(-16.0*((X-0.25)**2+(Y-0.75)**2)) 
        RETURN 
C 
C    Function 7: Fan of ripples. 
C 
    7   FUNC = SIN(6.283185*X*(Y+1.0)) 
        RETURN 
C 
C    Function 8: Concentric circular ripples. 
C 
    8   FUNC = COS(12.566371*SQRT((X-0.25)**2+(Y-0.25)**2)) 
        RETURN 
        END 
