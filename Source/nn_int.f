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

!/
!This subroutine NNINT was modified by Frank Tsai for the main program Lithology.90
!/

	SUBROUTINE NNINT(NBP,XX,YY,VAL,XMIN2,YMIN2,XMAX2,YMAX2,
	1  NCELL,XCELL,YCELL,VCELL)
C       INTEGER*2 L 
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
        DATA MSEED/97129/ 
        DATA NFMAX/8/ 
        DATA EPSCN,EPSPT/0.0,0.0/ 
        DATA NSTART/0/ 
        DATA NWRITE,NSTORE/12,13/ 

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
!        WRITE(NWRITE,9000)IFLAG 
2      Continue    

	VSHAPE=1.0

	DO I=1,NCELL
	XPT=XCELL(I)
	YPT=YCELL(I)

    4   CALL NNBR0(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE,VAL,Z0
     3   ,VSHAPE)

	NearestWell(I)=NINDEX
	VCELL(I)= Z0
	IF (NINDEX.LT.0) THEN
	WRITE(*,*)
	ENDIF

	ENDDO	
	
        IF(IFLAG.EQ.0) GOTO 5 
        WRITE(NWRITE,9000)IFLAG 
        write(NWRITE,*) XPT,YPT
        STOP 615 
C 
C    Run complete. 
C 
  5       Continue

	RETURN
 8000   FORMAT(3(I13,2X)) 
 8001   FORMAT(4(E13.6,2X)) 
 9000   FORMAT(///1H ,19(1H*)/1H ,16H ERROR: IFLAG = ,I3/1H ,19(1H*)///) 
 9001   FORMAT(///1H ,10HFUNCTION: ,I3/1H ,16HRMS ERROR:      ,E13.6/ 
     1    1H ,16HABSOLUTE ERROR: ,E13.6/1H ,16HMAXIMUM ERROR:  ,E13.6 
     2    ///) 
        END 
C 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    GENERAL ROUTINES LAST MODIFIED 16 MAY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    FUNCTION UNI4M 4.4 DATED 22 DECEMBER 1980 
C 
        FUNCTION UNI4M(MSEED) 
C 
C    Returns a pseudorandom number from the uniform distribution 
C    on (0,1).  The value of MSEED is used to generate the number 
C    and then modified.  This allows repeated calls to the function to 
C    generate a sequence of different random numbers, all depending on 
C    the initial value of MSEED. 
C 
C    The generator used is: 
C 
C      MSEED=MOD(24298*MSEED+99991,199017) 
C 
C    broken up so as to be safe in 24 bit integers. 
C    This is the generator used on Texas Instruments TI58/59 
C    calculators. 
C 
C    The user will usually supply UNI4M via the argument PRNG to some 
C    of the tile routines. 
C 
C    This generator is a convenient and simple one which is provided 
C    so as to make programs self-contained and portable.  Its 
C    statistical properties are respectable, but by no means as good 
C    as those of more sophisticated generators, and its cycle length 
C    is relatively short.  It should NOT normally be used for 
C    statistical simulation, for which purpose a generator such as 
C    that provided in the NAG library is to be preferred.  This 
C    generator is supplied with TILE 4 so as to facilitate 
C    the provision of test and demonstration programs. 
C 
C    Note: beware the fact that this function has a side-effect. 
C    Since it is used via EXTERNAL declarations, even the most 
C    arrogant optimising compiler should be wary, but you never 
C    know. 
C 
        MSEED = MOD(23*MOD(32*MOD(33*MOD(MSEED,199017),199017),199017)+ 
     1      10*MSEED+99991,199017) 
        UNI4M = (FLOAT(MSEED)+0.5)/199017.0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SHOVEL 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SHOVEL(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,IGARB,NPRINT) 
C 
C    Prints out the current status of the TILE4 database on stream 
C    NPRINT with carriage control characters. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
        WRITE(NPRINT,9001) LTOP,LBASE,LPTR,IGARB 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        WRITE(NPRINT,9002) JCNS,(L(LL), LL = LLO,LHI) 
        WRITE(NPRINT,9003) 
        DO 1  JCN = 1,JCNS 
        IF(JADDR(JCN).EQ.0) GOTO 2 
        LLO = JADDR(JCN)+2 
        LHI = LLO-1+L(LLO-1) 
        WRITE(NPRINT,9004) JCN,CN(1,JCN),CN(2,JCN),CN(3,JCN), 
     1    (L(LL), LL=LLO,LHI) 
        GO TO 1 
    2   WRITE(NPRINT,9005) JCN,CN(1,JCN),CN(2,JCN),CN(3,JCN) 
    1   CONTINUE 
        WRITE(NPRINT,9006) NPTS,NPTSIN,NFREE,NSTART 
        DO 3  NPT = 1,NPTS 
        LLL = NADDR(NPT) 
        IF(LLL) 4,701,5 
  701   STOP 747 
    5   LLO = LLL+2 
        LHI = L(LLO-1)+LLO-1 
        WRITE(NPRINT,9007) NPT,PT(1,NPT),PT(2,NPT), 
     1    (L(LL), LL = LLO,LHI) 
        GOTO 3 
    4   LLL = -LLL 
        WRITE(NPRINT,9008) NPT,LLL 
    3   CONTINUE 
        RETURN 
 9001   FORMAT(1H //1H ,14HFORTRAN TILE 4//1H ,10HHEAP SIZE ,I6/ 
     1    1H ,18HBASE OF WORKSPACE ,I4/1H ,19HBASE OF FREE SPACE ,I6/ 
     2    1H ,37HNUMBER OF FORCED GARBAGE COLLECTIONS ,I4) 
 9002   FORMAT(//1H ,22HNUMBER OF CONSTRAINTS ,I4//1H , 
     1    15HBOUNDARY LIST: ,10I6,(/1H ,15X,10I6)) 
 9003   FORMAT(//1H ,27HCONSTRAINT CONTIGUITY LISTS/) 
 9004   FORMAT(1H ,I4,3H : ,3F10.4,3H : ,10I6,(/1H ,40X,10I6)) 
 9005   FORMAT(1H ,I4,3H : ,3F10.4,12H : REDUNDANT) 
 9006   FORMAT(//1H ,26HNUMBER OF POINT LOCATIONS ,I6/ 
     1    1H ,26HNUMBER OF ACCEPTED POINTS ,I6/ 
     2    1H ,26HFREE LOCATION CHAIN ENTRY ,I6/ 
     3    1H ,20HWALK STARTING POINT ,6X,I6/// 
     4    1H ,22HPOINT CONTIGUITY LISTS/) 
 9007   FORMAT(1H ,I6,3H : ,2F10.4,3H : ,10I6,(/1H ,32X,10I6)) 
 9008   FORMAT(1H ,I6,20H : FREE - POINTS TO ,I6) 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RANPT 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE RANPT(VTX,KVTX,PROP,PRNG,MSEED,XPT,YPT) 
C 
C    The KVTX vertices of a convex polygon are held in (either) cyclic 
C    order in (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX, and the cumulative 
C    proportional subareas of a triangular subdivision are held in 
C    PROP(K) for K = 3,...,KVTX.  PROP can be loaded by a call to 
C    SUBDIV.  PRNG should be a U(0,1) pseudorandom number generator 
C    returning a single-precision real value and using and updating 
C    an integer seed supplied as MSEED.  The routine RANPT returns as 
C    (XPT,YPT) the coordinates of a uniform pseudorandom point in 
C    the polygon.  Normally three calls to PRNG are made. 
C 
        EXTERNAL PRNG 
        DIMENSION VTX(2,KVTX),PROP(KVTX) 
        IF(KVTX.LT.3) STOP 750 
        DO 1  J = 1,10 
        R1 = PRNG(MSEED) 
        R2 = PRNG(MSEED) 
        R3 = PRNG(MSEED) 
        IF(R1+R2-1.0) 2,1,3 
    3   R1 = 1.0-R1 
        R2 = 1.0-R2 
    2   DO 4  K = 3,KVTX 
        IF(PROP(K).GT.R3) GOTO 5 
    4   CONTINUE 
        STOP 751 
    5   XPT = R1*VTX(1,K-1)+R2*VTX(1,K)+(1.0-R1-R2)*VTX(1,1) 
        YPT = R1*VTX(2,K-1)+R2*VTX(2,K)+(1.0-R1-R2)*VTX(2,1) 
        RETURN 
    1   CONTINUE 
        STOP 752 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SCATRC 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SCATRC(CN,JCNS,PT,NPTS,IFLAG,PRNG,MSEED, 
     1    XMIN,XMAX,YMIN,YMAX) 
C 
C    Sets up entries in CN to specify as the window the rectangle 
C 
C       (XMIN,XMAX) * (YMIN,YMAX) 
C 
C    and then uses the external pseudorandom number generator 
C    PRNG, seeded by MSEED which it, and hence SCATRC, modifies, 
C    to set up entries in PT for NPTS points in random positions. 
C    IFLAG is zero on successful return.  Nonzero values indicate 
C    error return, as follows. 
C 
C       2  JCNS less than 4 
C       3  XMIN.GE.XMAX or YMIN.GE.YMAX 
C 
        EXTERNAL PRNG 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
C 
C    Check dimensionality. 
C 
        IF(NPTS.LE.0) STOP 753 
C 
C    Call SETREC to establish the window. 
C 
        CALL SETREC(CN,JCNS,IFLAG,XMIN,XMAX,YMIN,YMAX) 
C 
C    Check IFLAG on return from SETREC. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Set up the points using PRNG. 
C 
    1   DO 2  N = 1,NPTS 
        PT(1,N) = XMIN+(XMAX-XMIN)*PRNG(MSEED) 
    2   PT(2,N) = YMIN+(YMAX-YMIN)*PRNG(MSEED) 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SETREC 4.1 DATED 5 AUGUST 1980 
C 
        SUBROUTINE SETREC(CN,JCNS,IFLAG,XMIN,XMAX,YMIN,YMAX) 
C 
C    Sets up entries in CN to specify as the window the rectangle 
C 
C       (XMIN,XMAX) * (YMIN,YMAX). 
C 
C    IFLAG is zero on successful return.  Nonzero values indicate 
C    errors, as follows. 
C 
C       2   JCNS less than 4 
C       3   XMIN.GE.XMAX or YMIN.GE.YMAX 
C 
        DIMENSION CN(3,JCNS) 
C 
C    Check dimensionality. 
C 
        IF(JCNS.GE.4) GOTO 1 
        IFLAG = 2 
        RETURN 
C 
C    Check window specification. 
C 
    1   IF(XMIN.LT.XMAX.AND.YMIN.LT.YMAX) GOTO 2 
        IFLAG = 3 
        RETURN 
C 
C    Set up the window. 
C 
    2   CN(1,1) = 1.0 
        CN(2,1) = 0.0 
        CN(3,1) = -XMAX 
        CN(1,2) = 0.0 
        CN(2,2) = -1.0 
        CN(3,2) = YMIN 
        CN(1,3) = -1.0 
        CN(2,3) = 0.0 
        CN(3,3) = XMIN 
        CN(1,4) = 0.0 
        CN(2,4) = 1.0 
        CN(3,4) = -YMAX 
        IF(JCNS.EQ.4) GOTO 3 
        DO 4  J = 5,JCNS 
        CN(1,J) = 1.0 
        CN(2,J) = 0.0 
    4   CN(3,J) = -(XMAX+1.0) 
C 
C    Set IFLAG to zero and return. 
C 
    3   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TDUMP 4.2 DATED 16 FEBRUARY 1981 
C 
        SUBROUTINE TDUMP(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,NDUMP) 
C 
C    This subroutine dumps the entire tessellation data structure 
C    to device NDUMP using unformatted write statements.  It may 
C    then be read back subsequently using subroutine TLOAD (q.v.). 
C    The user is recommended to perform a garbage collection by a 
C    call to the routine GARBAJ immediately prior to a call to 
C    TDUMP in order to minimise the quantity of output. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Dump the simple INTEGERs. 
C 
        WRITE(NDUMP) JCNS,NPTS,NFREE,NSTART,NPTSIN,LTOP,LPTR,LBASE 
C 
C    Dump the simple REALs. 
C 
        WRITE(NDUMP) EPSCN,EPSPT 
C 
C    Dump the INTEGER arrays other than L. 
C 
        WRITE(NDUMP) (JADDR(J),J = 1,JCNS),(NADDR(N),N = 1,NPTS) 
C 
C    Dump L up to LPTR-1. 
C 
        LPT1 = LPTR-1 
        WRITE(NDUMP) (L(LL),LL = 1,LPT1) 
C 
C    Dump the REAL arrays. 
C 
        WRITE(NDUMP) ((CN(I,J),I = 1,3),J = 1,JCNS),(PT(1,N),PT(2,N), 
     1    N = 1,NPTS) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TLOAD 4.3 DATED 16 FEBRUARY 1981 
C 
        SUBROUTINE TLOAD(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,NLOAD) 
C 
C    This subroutine loads the entire tessellation data structure 
C    from device NLOAD using unformatted read statements. 
C    The routine checks that the dimensions of the 
C    arrays CN, JADDR, PT, NADDR are the same 
C    as they were when the file on device NLOAD was created 
C    using routine TDUMP (q.v.) and that the array L is large enough. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Load the simple INTEGERs and check the array dimensions. 
C 
        READ(NLOAD) JCNS1,NPTS1,NFREE,NSTART,NPTSIN,LTOP1,LPTR,LBASE 
        IF(JCNS1.NE.JCNS.OR.NPTS1.NE.NPTS.OR.LPTR.GT.LTOP) STOP 754 
C 
C    Load the simple REALs. 
C 
        READ(NLOAD) EPSCN,EPSPT 
C 
C    Load the INTEGER arrays other than L. 
C 
        READ(NLOAD) (JADDR(J),J = 1,JCNS),(NADDR(N),N = 1,NPTS) 
C 
C    Load L up to LPTR-1 
C 
        LPT1 = LPTR-1 
        READ(NLOAD) (L(LL),LL = 1,LPT1) 
C 
C    Load the REAL arrays. 
C 
        READ(NLOAD) ((CN(I,J),I = 1,3),J = 1,JCNS),(PT(1,N),PT(2,N), 
     1    N = 1,NPTS) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CNLD 4.2  DATED 16 MAY 1981 
C 
        SUBROUTINE CNLD(CN,JCNS,VTX,KVTX) 
C 
C    Subroutine to compute the window constraints from the coordinates 
C    of the window vertices. 
C 
C    The window vertex coordinates should be loaded by the user into 
C    (VTX(1,K),VTX(2,K)) K = 1......KVTX in clockwise order 
C    round the window.  On return the array CN will contain the proper 
C    description of the window.  The only error detected by the routine 
C    is the setting of KVTX to a greater value than JCNS.  This causes 
C    a trap to a numbered STOP statement.  The user should take care 
C    that the vertices that he specifies do form a convex window.  If 
C    they do not this will usually be trapped when WINLD is subsequently 
C    called, but this is not inevitable, as the intersecting half 
C    planes defined by a non-convex set of vertices may still form 
C    a valid, though incorrect, window. 
C 
        DIMENSION CN(3,JCNS),VTX(2,KVTX) 
C 
C    Check for a DIMENSION error. 
C 
        IF(KVTX.GT.JCNS.OR.KVTX.LT.3) STOP 755 
C 
C    Set up the first vertex and keep track of the minimum and maximum 
C    vertex x coordinates. 
C 
        XOLD = VTX(1,1) 
        YOLD = VTX(2,1) 
        XMX = XOLD 
        XMN = XMX 
C 
C    Scan through the vertices setting up CN. 
C 
        DO 1 K = 1,KVTX 
        IF(K.NE.KVTX) GOTO 2 
        XNEW = VTX(1,1) 
        YNEW = VTX(2,1) 
        GOTO 3 
    2   K1 = K+1 
        XNEW = VTX(1,K1) 
        YNEW = VTX(2,K1) 
        IF(XNEW.GT.XMX) XMX = XNEW 
        IF(XNEW.LT.XMN) XMN = XNEW 
    3   A = YOLD-YNEW 
        B = XNEW-XOLD 
C 
C    Normalise the gradient terms. 
C 
        SCA = 1.0/SQRT(A*A+B*B) 
        A = A*SCA 
        B = B*SCA 
        CN(3,K) = -A*XOLD-B*YOLD 
        CN(1,K) = A 
        CN(2,K) = B 
        XOLD = XNEW 
    1   YOLD = YNEW 
C 
C    Fill any extra CN entries with a harmless constraint. 
C 
        IF(JCNS.EQ.KVTX) RETURN 
        XMX = XMN-2.0*XMX 
        KVT1 = KVTX+1 
        DO 4 K = KVT1,JCNS 
        CN(1,K) = 1.0 
        CN(2,K) = 0.0 
    4   CN(3,K) = XMX 
        RETURN 
        END 
C 
C 
C*********************************************************************** 
C 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    MANIPULATION ROUTINES LAST MODIFIED 28 JULY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILE 4.6 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TILE(CN,JCNS,JADDR,PT,NPTS,NADDR,L,LTOP, 
     1    LPTR,IFLAG,IGARB) 
C 
C    A simple master routine for constructing a tessellation 
C    within an arbitrary window. 
C 
C    The subroutine is compatible with that offered as the 
C    master subroutine in TILE 003.  The calling program is 
C    responsible for entering the constraints into CN and the 
C    points into PT, and for setting the dimensions JCNS,NPTS, 
C    LTOP; JCNS serves also as the number of constraints and 
C    NPTS as the number of points.  Information is returned 
C    through the standard data structure in the heap L (which 
C    is garbaged immediately before return) accessed through 
C    the address vectors JADDR,NADDR.  LPTR is the base of the 
C    free area in L, IGARB the number of garbage collections 
C    excluding the final one.  IFLAG is zero on successful 
C    return.  A nonzero value of IFLAG indicates an error 
C    condition, as follows. 
C 
C       1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C       2  Unbounded window - constraints inadequate 
C       3  Empty window - constraints inconsistent 
C       4  No accepted points within window 
C       5  Attempt to insert duplicate point 
C       6  Heap overflow 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Hard error for nonpositive dimensions. 
C 
        IF(JCNS.LE.0.OR.NPTS.LE.0.OR.LTOP.LE.0) STOP 701 
C 
C    Set up the window and check its validity. 
C 
        CALL WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Initialise other variables. 
C 
    1   CALL CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    LPTR,LBASE,IGARB) 
C 
C    Add the points, checking the validity of each, and 
C    occasionally resetting the starting-point of the 
C    nearest-neighbour walk to keep it near to the centroid 
C    of the accepted points. 
C 
        NRESET = 2 
        DO 2  N = 1,NPTS 
        XPT = PT(1,N) 
        YPT = PT(2,N) 
        CALL ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,0.0,0.0, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
        IF(IFLAG.EQ.7) STOP 702 
        IF(IFLAG.EQ.8) STOP 703 
        IF(IFLAG.NE.5.AND.IFLAG.NE.6) GOTO 3 
        RETURN 
    3   IF(IFLAG.NE.4) GOTO 4 
        NSAVE = -NADDR(NPTS) 
        NADDR(NPTS) = -N 
        NFREE = -NADDR(N) 
        NADDR(N) = -NSAVE 
        GOTO 2 
    4   IF(NPTSIN.LT.NRESET) GOTO 2 
        NRESET = 2*NRESET 
        CALL MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP, 
     1    IFLAG,XBAR,YBAR,NINDEX) 
        IF(IFLAG.EQ.4) STOP 704 
        IF(IFLAG.EQ.8) STOP 705 
        NSTART = NINDEX 
    2   CONTINUE 
C 
C    Final garbage collection and return. 
C 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1981 
C    SUBROUTINE WINLD 4.5 DATED 28 SEPTEMBER 1981 
C 
        SUBROUTINE WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
C 
C    Constructs the description of the window from the 
C    given constraints. 
C 
C    The subroutine checks the validity of the window.  If 
C    it is valid, then a data item is constructed as follows 
C    at the bottom of the heap L: L(1) is unset; L(2) is the 
C    number of effective constraints; succeeding locations in 
C    L, as many as are needed, hold the list of effective 
C    constraints in clockwise cyclic order round the boundary 
C    of the window, each pointer back to a constraint being 
C    flagged with a negative sign.  LBASE is set to the next 
C    free location in L above this data item.  Entries in JADDR 
C    are set to positive values for effective constraints and 
C    to zero for redundant ones.  On successful exit from 
C    the subroutine, IFLAG returns the value zero.  Four 
C    possible nonzero values may be returned, each denoting 
C    a catastrophic error condition, as follows. 
C 
C       1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C       2  Unbounded window - constraints inadequate 
C       3  Empty window - constraints inconsistent 
C       6  Heap overflow 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),L(LTOP) 
        LOGICAL NEWCLK,NEWACK 
C 
C    Check for dud values of LTOP and JCNS. 
C 
        IF(LTOP.LT.5) GOTO 66 
        IF(JCNS.LT.2) GOTO 62 
C 
C    Initialise JADDR to zero and check for dud constraints. 
C 
        DO 2 J = 1,JCNS 
        JADDR(J) = 0 
        IF(CN(1,J).EQ.0.0.AND.CN(2,J).EQ.0.0) GOTO 61 
    2   CONTINUE 
C 
C    Find the first effective constraint. 
C 
        JCOUNT = 0 
        LBASE = 3 
        DO 3 J = 1,JCNS 
        A = CN(1,J) 
        B = CN(2,J) 
        C = CN(3,J) 
        NEWCLK = .FALSE. 
        NEWACK = .FALSE. 
        DO 4 J1 = 1,JCNS 
        IF(JADDR(J1).NE.0.OR.J.EQ.J1) GOTO 4 
        CALL UU(A,B,C,CN(1,J1),CN(2,J1),CN(3,J1),U,IWIPE) 
        IF(IWIPE) 5,6,7 
C 
C    Check parallel constraints. 
C 
    6   AMDA = (A*CN(1,J1)+B*CN(2,J1))/(A*A+B*B) 
        IF(U) 27,27,28 
   27   IF(AMDA) 63,61,31 
   28   IF(AMDA) 4,61,32 
C 
C    Constraint J1 is made redundant by constraint J. 
C 
   32   JADDR(J1) = -1 
        GOTO 4 
C 
C    The wipe is anti-clockwise. 
C 
    5   IF(NEWACK) GOTO 9 
        UACK = U 
        NEWACK = .TRUE. 
        GOTO 4 
    9   IF(U.GT.UACK) UACK = U 
        GOTO 4 
C 
C    The wipe is clockwise. 
C 
    7   IF(NEWCLK) GOTO 10 
        UCLK = U 
        JNEXT = J1 
        NEWCLK = .TRUE. 
        GOTO 4 
   10   IF(U.GE.UCLK) GOTO 4 
        UCLK = U 
        JNEXT = J1 
    4   CONTINUE 
C 
C    Check for unbounded window. 
C 
        IF(.NOT.(NEWCLK.AND.NEWACK)) GOTO 62 
C 
C    Check if the wipes overlap.  If they do the constraint is 
C    ineffective, as nowhere along its length does it form part 
C    of the window boundary. 
C 
        IF(UCLK.GT.UACK) GOTO 26 
C 
C    The constraint is ineffective. 
C 
   31   JADDR(J) = -1 
        JCOUNT = JCOUNT+1 
    3   CONTINUE 
C 
C    If this point is reached no constraint is effective: the window 
C    is empty. 
C 
C 
C    The first effective constraint has been found.  It is J. 
C 
   26   JFIRST = J 
C 
C    Consider the next constraint in a clockwise direction. 
C 
   11   L(LBASE) = -J 
        JADDR(J) = 1 
        LBASE = LBASE+1 
        IF(LBASE.GT.LTOP) GOTO 66 
        JCOUNT = JCOUNT+1 
C 
C    Check for cycling round and round. 
C 
        IF(JCOUNT.GE.JCNS) STOP 706 
C 
C    Compute the anticlockwise wipe of constraint J on JNEXT. 
C 
        CALL UU(CN(1,JNEXT),CN(2,JNEXT),CN(3,JNEXT),A,B,C,U,IWIPE) 
        IF(IWIPE) 13,777,777 
  777   STOP 707 
   13   UACK = U 
C 
C    Find the next clockwise constraint from JNEXT. 
C 
        A = CN(1,JNEXT) 
        B = CN(2,JNEXT) 
        C = CN(3,JNEXT) 
        NEWCLK = .FALSE. 
        DO 14 J1 = 1,JCNS 
        IF(J1.EQ.JFIRST) GOTO 12 
        IF(JADDR(J1).NE.0.OR.J1.EQ.JNEXT) GOTO 14 
   12   CALL UU(A,B,C,CN(1,J1),CN(2,J1),CN(3,J1),U,IWIPE) 
        IF(IWIPE) 15,16,17 
C 
C    Check parallel constraints. 
C 
   16   AMDA = (A*CN(1,J1)+B*CN(2,J1))/(A*A+B*B) 
        IF(U) 127,127,128 
  127   IF(AMDA) 63,61,710 
  128   IF(AMDA) 14,61,14 
C 
C    Anti-clockwise wipe - check for inconsistency. 
C 
   15   IF(J1.EQ.JFIRST) GOTO 14 
        IF(U.LE.UACK) GOTO 14 
  710   STOP 710 
C 
C    Clockwise wipe. 
C 
   17   IF(NEWCLK) GOTO 20 
        JNXT2 = J1 
        UCLK = U 
        NEWCLK = .TRUE. 
        GOTO 14 
   20   IF(U.GE.UCLK) GOTO 14 
        UCLK = U 
        JNXT2 = J1 
   14   CONTINUE 
C 
C    Check for unbounded window. 
C 
        IF(.NOT.NEWCLK) GOTO 62 
C 
C    Check if the window has been closed. 
C 
        IF(JNXT2.EQ.JFIRST) GOTO 22 
C 
C    It has not.  Go back and add JNEXT to the list. 
        J = JNEXT 
        JNEXT = JNXT2 
        GOTO 11 
C 
C    It has.  Add JNEXT as the last constraint. 
C 
   22   L(LBASE) = -JNEXT 
        JADDR(JNEXT) = 1 
        LBASE = LBASE+1 
        IF(LBASE.GT.LTOP) GOTO 66 
        JCNSIN = LBASE-3 
        L(2) = JCNSIN 
C 
C    Entries in JADDR that are still 0 must be ineffective. 
C    Find any initial ineffective constraints and set them 
C    to zero as well. 
C 
        DO 23 J = 1,JCNS 
        IF(JADDR(J).GT.0) GOTO 60 
   23   JADDR(J) = 0 
        GOTO 63 
C 
C    Successful return. 
C 
   60   IFLAG = 0 
        RETURN 
C 
C    Error returns. 
C 
   61   IFLAG = 1 
        RETURN 
   62   IFLAG = 2 
        RETURN 
   63   IFLAG = 3 
        RETURN 
   66   IFLAG = 6 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CLEAR 4.1 DATED 17 AUGUST 1977 
C 
        SUBROUTINE CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    LPTR,LBASE,IGARB) 
C 
C    Initialises or reinitialises that part of the data 
C    structure not set up by subroutine WINLD. 
C 
C    The subroutine constructs a chain through the 
C    locations in NADDR with NFREE pointing to the start 
C    of the chain.  NSTART, NPTSIN, and IGARB are 
C    initialised to zero, and LPTR to LBASE, which has 
C    already been set up.  The effect of a call to 
C    subroutine CLEAR is to wipe the window clear of 
C    points ready for the construction of a completely 
C    new tessellation within the same window. 
C 
        DIMENSION NADDR(NPTS) 
        DO 1  N = 1,NPTS 
    1   NADDR(N) = -(N+1) 
        NFREE = 1 
        NSTART = 0 
        NPTSIN = 0 
        LPTR = LBASE 
        IGARB = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ADDPT 4.5 DATED 13 JANUARY 1981 
C 
        SUBROUTINE ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Builds the point (XPT,YPT) into the tessellation and 
C    returns its index as NINDEX. 
C 
C    The subroutine attempts to build the new point whose 
C    coordinates are (XPT,YPT) into the tessellation.  A 
C    check is made that the new point is adequately inside 
C    the window.  Then the current data structure is 
C    interrogated to determine whether there are already 
C    any accepted points; if not, all that needs to be done 
C    is some appropriate initialisation.  If points are 
C    already accepted, the (a) nearest neighbour among them 
C    to the new point is found, and a check is made that 
C    the new point is not too close to it.  A free index is 
C    reserved, and NINDEX set to its value; the corresponding 
C    entries in PT are set to the coordinates of the new 
C    point.  Finally, the contiguity list for the new point 
C    is constructed, and those of its neighbours are 
C    modified, so as to modify the tessellation to include 
C    it.  On successful completion of this procedure, IFLAG 
C    returns the value zero.  Nonzero values indicate 
C    failure to insert the new point, as follows. 
C 
C       4  Attempt to insert point outside window 
C       5  Attempt to insert duplicate point 
C       6  Heap overflow *** catastrophic error 
C       7  Address vector already full 
C       8  NSTART not the index of an accepted point 
C 
C    "Outside" and "duplicate" are each interpreted in terms 
C    of the tolerance values EPSCN,EPSPT in the argument 
C    list.  If there are no previously accepted points, 
C    return with IFLAG set to 8 cannot of course occur; 
C    in that case NSTART is set up to the index of the 
C    newly-accepted first point, and otherwise it is left 
C    unchanged. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Call subroutine SEEPT to check that (XPT,YPT) is far 
C    enough inside the window; return if not, with NINDEX 
C    giving the index of an unsatisfied constraint. 
C 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    If any points are already accepted, jump forward 
C    to the main part of the routine. 
C 
    1   IF(NPTSIN.GT.0) GOTO 2 
C 
C    This is the first point.  Check that NFREE defines a 
C    valid location for it, set NINDEX to this value, and 
C    copy the coordinates of the new point to the corresponding 
C    entries in PT.  Reset NPTSIN to 1. 
C 
        IF(NFREE.LE.NPTS) GOTO 3 
        IFLAG = 7 
        RETURN 
    3   NINDEX = NFREE 
        NFREE = -NADDR(NFREE) 
        PT(1,NINDEX) = XPT 
        PT(2,NINDEX) = YPT 
        NPTSIN = 1 
C 
C    Check that there is enough space available in the heap 
C    for initialisation to be possible. 
C 
        IF(LPTR+6*L(2)+1.LE.LTOP) GOTO 4 
        IFLAG = 6 
        RETURN 
C 
C    LPTR is the first free location in the heap.  Build up 
C    contiguity lists from that location; first, that for 
C    the new point NINDEX, 
C 
    4   NADDR(NINDEX) = LPTR 
        L(LPTR) = NINDEX 
        L(LPTR+1) = L(2) 
        LPTR = LPTR+2 
        LHI = 2+L(2) 
        DO 5  LL = 3,LHI 
        LC = LHI+3-LL 
        L(LPTR) = L(LC) 
    5   LPTR = LPTR+1 
C 
C    next, that for the first effective constraint, 
C 
        J = -L(3) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LHI) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(4) 
        LPTR = LPTR+5 
C 
C    next, that for the last effective constraint, 
C 
        J = -L(LHI) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LHI-1) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(3) 
        LPTR = LPTR+5 
C 
C    and finally, for the remaining effective constraints; 
C    since the window is bounded there are at least three 
C    effective constraints. 
C 
        LHI = LHI-1 
        DO 6  LL = 4,LHI 
        J = -L(LL) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LL-1) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(LL+1) 
    6   LPTR = LPTR+5 
C 
C    Set NSTART to NINDEX, IFLAG to zero, and return. 
C 
        NSTART = NINDEX 
        IFLAG = 0 
        RETURN 
C 
C    The main part of the routine deals with the insertion 
C    of a new point when there are already some accepted 
C    points.  Call subroutine LOCPT to find the (a) nearest 
C    neighbour of the new point (XPT,YPT).  If it is too 
C    close, return with the index of the neighbour in NINDEX. 
C 
    2   CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 7 
        RETURN 
C 
C    Check that NFREE defines a valid location for the new 
C    point, set NINDEX to this value having saved the index 
C    of the nearest neighbour as NEAR, and copy the coordinates 
C    of the new point to the corresponding entries in PT. 
C    Increment NPTSIN. 
C 
    7   IF(NFREE.LE.NPTS) GOTO 8 
        IFLAG = 7 
        RETURN 
    8   NEAR = NINDEX 
        NINDEX = NFREE 
        NFREE = -NADDR(NFREE) 
        PT(1,NINDEX) = XPT 
        PT(2,NINDEX) = YPT 
        NPTSIN = NPTSIN+1 
C 
C    Initialise NNEW, the first item to be placed on the 
C    contiguity list for NINDEX; it is the neighbour of 
C    NINDEX clockwise from NEAR.  Flag degeneracy if it occurs. 
C 
        LLO = NADDR(NEAR)+2 
        LHI = LLO-1+L(LLO-1) 
        XNEAR = PT(1,NEAR) 
        YNEAR = PT(2,NEAR) 
        NNEW = 0 
        DO 9  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XNEAR,YNEAR,N,T,IWIPE) 
        IF(IWIPE.LE.0) GOTO 9 
        IF(NNEW.EQ.0) GOTO 10 
        IF(TNEW-T) 9,51,10 
   51   IDGNEW = 1 
        IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 9 
   10   IDGNEW = 0 
        NNEW = N 
        TNEW = T 
    9   CONTINUE 
C 
C    Initialise NCURR to NEAR and LEFF to LTOP.  The 
C    contiguity list for NINDEX is constructed in a 
C    temporary position working down from LTOP. 
C 
        NCURR = NEAR 
        LEFF = LTOP 
C 
C    Enter the main loop.  First update NOLD and NCURR, 
C    and pass back the degeneracy flag. 
C 
   11   NOLD = NCURR 
        NCURR = NNEW 
        IDGCUR = IDGNEW 
C 
C    Is NCURR a point or a constraint? 
C 
        IF(NCURR) 12,701,13 
  701   STOP 712 
C 
C    Constraint.  (XCURR,YCURR) is the reflexion of NINDEX 
C    in it. 
C 
   12   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XPT+BJ*YPT+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XPT-AJ*DJ 
        YCURR = YPT-BJ*DJ 
C 
C    Pick up beginning and size of contiguity list, and 
C    garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = JADDR(JCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        IFLAG = 6 
        RETURN 
C 
C    Point, with coordinates (XCURR,YCURR). 
C 
   13   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
C 
C    Pick up beginning and size of contiguity list, and 
C    garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = NADDR(NCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        IFLAG = 6 
        RETURN 
C 
C    Find the end of the old contiguity list for NCURR, and flag 
C    this list as defunct.  Save the value of LPTR as the address 
C    of the new contiguity list for NCURR, flag the new 
C    list as active, and place NINDEX on it. 
C 
   14   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LCURR = LPTR 
        L(LPTR) = NCURR 
        L(LPTR+2) = NINDEX 
        LPTR = LPTR+3 
C 
C    Contiguity is a symmetric relation: find NOLD in 
C    the old contiguity list for NCURR. 
C 
        DO 15  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 16 
   15   CONTINUE 
        STOP 713 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for NINDEX after NCURR; it is the neighbour of 
C    NINDEX clockwise from NCURR.  Flag degeneracy if it 
C    occurs, and backspace LOLD if degeneracy occurred 
C    last time. 
C 
   16   NNEW = 0 
        IDGNEW = 0 
        LL = LOLD 
        IF(IDGCUR.EQ.0) GOTO 17 
        LOLD = LOLD-1 
        IF(LOLD.LT.LLO) LOLD = LHI 
   17   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.GT.0) GOTO 18 
        IF(NNEW) 19,17,19 
   18   IF(NNEW.EQ.0) GOTO 20 
        IF(TNEW-T) 19,52,20 
   52   IDGNEW = 1 
        GOTO 21 
   20   NNEW = N 
        TNEW = T 
        GOTO 17 
C 
C    Complete the new contiguity list for NCURR, and set the 
C    appropriate pointer to it from either JADDR or NADDR. 
C 
   19   L(LPTR) = NNEW 
        LPTR = LPTR+1 
   21   L(LPTR) = N 
        LPTR = LPTR+1 
        IF(LL.EQ.LOLD) GOTO 22 
        LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        GOTO 21 
   22   L(LCURR+1) = LPTR-LCURR-2 
        IF(NCURR) 23,703,24 
  703   STOP 714 
   23   JADDR(JCURR) = LCURR 
        GOTO 25 
   24   NADDR(NCURR) = LCURR 
C 
C    Add NCURR to the contiguity list for NINDEX. 
C 
   25   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NEAR to the contiguity list for 
C    NINDEX, then the main loop is completed. 
C 
        IF(NCURR.NE.NEAR) GOTO 11 
C 
C    Move the contiguity list for NINDEX down to the base of 
C    the free area, and set up pointers and size. 
C 
        NADDR(NINDEX) = LPTR 
        L(LPTR) = NINDEX 
        L(LPTR+1) = LTOP-LEFF 
        LPTR = LPTR+2 
        LEFF = LEFF+1 
        DO 26  LL = LEFF,LTOP 
        L(LPTR) = L(LL) 
   26   LPTR = LPTR+1 
C 
C    The new point has been built into the tessellation 
C    successfully.  Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SEEPT 4.2 DATED 21 OCTOBER 1980 
C 
        SUBROUTINE SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    JINDEX,TSTVAL) 
C 
C    Tests whether (XPT,YPT) is inside the window. 
C 
C    The subroutine evaluates the perpendicular distance 
C    from (XPT,YPT) to each effective constraint, flagged 
C    with a negative sign on the inwards side.  TSTVAL returns 
C    the largest such value, JINDEX the constraint where it 
C    occurs - arbitrary choice if equality.  IFLAG returns 
C    the value zero unless any such value equals or exceeds 
C    zero, in which case it returns the value 4.  If EPSCN 
C    is positive, then -EPSCN is used as the cutoff 
C    value instead of zero. 
C 
        DIMENSION CN(3,JCNS),JADDR(JCNS) 
C 
C    Examine the constraints to find TSTVAL and JINDEX. 
C 
        JINDEX = 0 
        DO 1  J = 1,JCNS 
        IF(JADDR(J).LE.0) GOTO 1 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        T = (AJ*XPT+BJ*YPT+CN(3,J))/SQRT(AJ*AJ+BJ*BJ) 
        IF(JINDEX.EQ.0) GOTO 2 
        IF(T.LE.TSTVAL) GOTO 1 
    2   TSTVAL = T 
        JINDEX = -J 
    1   CONTINUE 
C 
C    Find the cutoff value. 
C 
        CUTOFF = 0.0 
        IF(EPSCN.GT.0.0) CUTOFF = -EPSCN 
C 
C    Report if (XPT,YPT) does not lie in the window. 
C 
        IF(TSTVAL.LT.CUTOFF) GOTO 3 
        IFLAG = 4 
        RETURN 
    3   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCPT 4.2 DATED 18 JULY 1979 
C 
        SUBROUTINE LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
C 
C    Finds the nearest neighbour of (XPT,YPT). 
C 
C    The subroutine finds the index of the accepted point 
C    closest to the trial point (XPT,YPT), which does 
C    not need to be inside the window.  NINDEX returns the 
C    index of the closest point, and TSTVAL its squared 
C    distance from the trial point - arbitrary choice if 
C    equality.  IFLAG returns the value zero unless TSTVAL 
C    is zero, in which case it returns the value 5.  If 
C    EPSPT is positive then 4.0*EPSPT*EPSPT is used as the 
C    cutoff value instead of zero.  If the starting-point 
C    for the walk by which the nearest neighbour is found 
C    is not the index of an accepted point, then IFLAG 
C    returns the value 8 and no other useful information 
C    is provided. 
C 
C2      INTEGER*2 L 
        DIMENSION PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Check that NSTART is the index of an accepted point. 
C 
        IF(NSTART.LE.0.OR.NSTART.GT.NPTS) GOTO 5 
        IF(NADDR(NSTART).GT.0) GOTO 1 
    5   IFLAG = 8 
        RETURN 
C 
C    Initialise NINDEX to NSTART and then replace successively 
C    by any nearer accepted point until the (a) nearest is 
C    found. 
C 
    1   NINDEX = NSTART 
        XDIFF = XPT-PT(1,NINDEX) 
        YDIFF = YPT-PT(2,NINDEX) 
        TSTVAL = XDIFF*XDIFF+YDIFF*YDIFF 
    2   LLO = NADDR(NINDEX)+2 
        LHI = LLO-1+L(LLO-1) 
        DO 3  LL = LLO,LHI 
        NEXT = L(LL) 
        IF(NEXT.LE.0) GOTO 3 
        XDIFF = XPT-PT(1,NEXT) 
        YDIFF = YPT-PT(2,NEXT) 
        T = XDIFF*XDIFF+YDIFF*YDIFF 
        IF(T.GE.TSTVAL) GOTO 3 
        NINDEX = NEXT 
        TSTVAL = T 
        GOTO 2 
    3   CONTINUE 
C 
C    Find the cutoff value. 
C 
        CUTOFF = 0.0 
        IF(EPSPT.GT.0.0) CUTOFF = 4.0*EPSPT*EPSPT 
C 
C    Test the trial point for closeness to its nearest 
C    neighbour. 
C 
        IF(TSTVAL.GT.CUTOFF) GOTO 4 
        IFLAG = 5 
        RETURN 
    4   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE MIDPT 4.2 DATED 3 JULY 1979 
C 
        SUBROUTINE MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP, 
     1    IFLAG,XBAR,YBAR,NINDEX) 
C 
C    Finds the centroid of the accepted points and the index 
C    of the nearest accepted point to it. 
C 
C    The subroutine returns the coordinates of the centroid 
C    of the accepted points as (XBAR,YBAR), and the index of 
C    the (a) nearest accepted point to it as NINDEX - arbitrary 
C    choice if equality.  IFLAG returns the value zero on 
C    successful completion of the calculation.  If NSTART is 
C    not the index of an accepted point then (XBAR,YBAR) is 
C    returned correctly, but IFLAG returns the value 8 and 
C    no useful information is returned in NINDEX.  If there 
C    are no accepted points then no useful information is 
C    returned; IFLAG returns the value 4. 
C 
C2      INTEGER*2 L 
        DIMENSION PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Check that there are some accepted points. 
C 
        IF(NPTSIN.GT.0) GOTO 1 
        IFLAG = 4 
        RETURN 
C 
C    Find the centroid of the accepted points. 
C 
    1   XBAR = 0.0 
        YBAR = 0.0 
        DO 2  N = 1,NPTS 
        IF(NADDR(N).LE.0) GOTO 2 
        XBAR = XBAR+PT(1,N) 
        YBAR = YBAR+PT(2,N) 
    2   CONTINUE 
        XBAR = XBAR/FLOAT(NPTSIN) 
        YBAR = YBAR/FLOAT(NPTSIN) 
C 
C    Call subroutine LOCPT to find the index of the nearest 
C    accepted point. 
C 
        CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,0.0, 
     1    IFLAG,XBAR,YBAR,NINDEX,TSTVAL) 
C 
C    Suppress the possible value 5 for IFLAG returned 
C    by this subroutine - it does not matter if the 
C    centroid coincides with an accepted point. 
C 
        IF(IFLAG.EQ.5) IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TT 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TT(CN,JCNS,PT,NPTS,X0,Y0,X1,Y1,N,T,IWIPE) 
C 
C    Calculates intercepts on the perpendicular bisector of 
C    the line joining (X0,Y0) and (X1,Y1). 
C 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
        IWIPE = 1 
        IF(N) 2,701,1 
  701   STOP 715 
    1   XN = PT(1,N) 
        YN = PT(2,N) 
        D = (XN-X0)*(Y1-Y0)-(YN-Y0)*(X1-X0) 
        T = (XN-X0)*(XN-X1)+(YN-Y0)*(YN-Y1) 
    6   IF(D) 3,4,5 
    3   IWIPE = -1 
    5   T = T/D 
        RETURN 
    4   IWIPE = 0 
        RETURN 
    2   J = -N 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        D = AJ*(Y1-Y0)-BJ*(X1-X0) 
        T = -AJ*(X0+X1)-BJ*(Y0+Y1)-2.0*CN(3,J) 
        GOTO 6 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE UU 4.1 DATED 22 JULY 1976 
C 
        SUBROUTINE UU(A0,B0,C0,A1,B1,C1,U,IWIPE) 
C 
C    Calculates intercepts on a constraint. 
C 
        IWIPE = 1 
        D = B0*A1-A0*B1 
        U = -C1+C0*(A0*A1+B0*B1)/(A0*A0+B0*B0) 
        IF(D) 1,2,3 
    1   IWIPE = -1 
    3   U = U/D 
        RETURN 
    2   IWIPE = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GARBAJ 4.1 DATED 17 AUGUST 1977 
C 
        SUBROUTINE GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR, 
     1    LPROT) 
C 
C    Does a garbage collection on the heap.  Free space is 
C    returned to the top end, and LPTR, which points to the 
C    base of the free area, is reset.  The area below LPROT 
C    is protected. 
C 
C2      INTEGER*2 L 
        DIMENSION JADDR(JCNS),NADDR(NPTS),L(LTOP) 
C 
C    Return if there is nothing to garbage. 
C 
        IF(LPTR.EQ.LPROT) RETURN 
C 
C    The old value of LPTR is saved so that completion of the 
C    scan can be detected.  LPTR is reset to LPROT, and LSCAN 
C    is initialised to the same value. 
C 
        LSAVE = LPTR 
        LPTR = LPROT 
        LSCAN = LPROT 
C 
C    On entry to the scanning loop, LSCAN points to the first 
C    cell of the current data item.  If the entry in this cell 
C    is zero, the data item is defunct, and is skipped. 
C    Otherwise the data item is active, being the contiguity 
C    list either of an effective constraint (negative entry) 
C    or of an accepted point (positive entry).  In either of 
C    these latter cases the appropriate pointer into the 
C    heap is updated. 
C 
    1   N = L(LSCAN) 
        LLO = LSCAN 
        LSCAN = LSCAN+L(LSCAN+1)+2 
        IF(N) 2,3,4 
    2   J = -N 
        JADDR(J) = LPTR 
        GOTO 5 
    4   NADDR(N) = LPTR 
C 
C    Copy down the data item. 
C 
    5   LHI = LSCAN-1 
        DO 6  LL = LLO,LHI 
        L(LPTR) = L(LL) 
    6   LPTR = LPTR+1 
C 
C    Loop if there is more to do. 
C 
    3   IF(LSCAN.LT.LSAVE) GOTO 1 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SUBPT 4.5 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SUBPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,IFLAG,IGARB, 
     2    XPT,YPT,NINDEX) 
C 
C    Removes from the tessellation the point whose index 
C    is NINDEX and returns its coordinates as (XPT,YPT). 
C 
C    The subroutine attempts to remove from the tessellation 
C    the point whose index is NINDEX.  A check is made that 
C    this index refers to a currently accepted point.  Then 
C    the current data structure is interrogated to determine 
C    whether it is the only one; if so, all that is needed 
C    to delete it is a simple reinitialisation procedure 
C    similar to that carried out by subroutine CLEAR, but 
C    preserving the value of IGARB and the pointer chain in 
C    NADDR.  If not, a substantial calculation is needed. 
C    What is done is to carry out a scratchpad calculation 
C    of the tessellation of the whole window determined by 
C    those points contiguous to NINDEX.  This calculation 
C    is done by logic identical to that employed to 
C    construct a tessellation by repeated calls to 
C    subroutine ADDPT; the only reason why ADDPT cannot be 
C    used to do it is the need to administer a slightly more 
C    complicated data structure to maintain the scratchpad 
C    calculation distinct from the main data structure. 
C    The objects whose data items include scratchpad 
C    sections at the completion of the scratchpad 
C    calculation are (a) all effective constraints, and (b) 
C    all points contiguous to the point being removed. 
C    The form of such an extended data item is as follows: 
C    (1) back-reference; (2) total length of following part 
C    of data item; (3) main contiguity list; (4) marker zero; 
C    (5) scratchpad contiguity list.  Because the back- 
C    reference and length are maintained as in the normal 
C    form of the data base, garbage collection is not 
C    affected.  Once the scratchpad calculation is 
C    complete, a merging operation is carried out whereby 
C    the main and scratchpad lists for objects contiguous 
C    to the point being deleted are merged to form new 
C    lists describing the pattern of contiguities after 
C    removal of the point being deleted.  Then the 
C    scratchpad lists for effective constraints not 
C    contiguous to the point being deleted are killed; the 
C    data item for the deleted point is killed; NINDEX 
C    is chained into the beginning of the free space in NADDR; 
C    and finally the count of accepted points is decremented by 
C    one.  On successful completion of this procedure, 
C    return takes place with IFLAG set to zero; nonzero 
C    values indicate failure to remove the specified point, 
C    as follows. 
C 
C       6  Heap overflow *** catastrophic error 
C       9  NINDEX not the index of an accepted point 
C 
C    The coordinates of the removed point are returned as	 
C    (XPT,YPT).  If the removal leaves no accepted points, 
C    NSTART is reset to zero.  If the point being removed 
C    happens to coincide with NSTART, then NSTART is reset 
C    to a neighbouring point.  Otherwise it is not checked. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Check that NINDEX is the index of an accepted point, 
C    and return with IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 100 
        IF(NADDR(NINDEX).GT.0) GOTO 1 
  100   IFLAG = 9 
        RETURN 
C 
C    Copy the coordinates of NINDEX to (XPT,YPT); these 
C    values are not needed again. 
C 
    1   XPT = PT(1,NINDEX) 
        YPT = PT(2,NINDEX) 
C 
C    If NINDEX is the only accepted point, all that is 
C    needed to remove it is to reinitialise NSTART, NPTSIN, 
C    and LPTR, and to chain NINDEX back into NADDR. 
C 
        IF(NPTSIN.GT.1) GOTO 2 
        NSTART = 0 
        NADDR(NINDEX) = -NFREE 
        NFREE = NINDEX 
        NPTSIN = 0 
        LPTR = LBASE 
        RETURN 
C 
C    If NINDEX is not the only accepted point, a lengthy 
C    calculation is needed to remove it.  First pick up 
C    the length of its contiguity list.  Also initialise 
C    to zero the number of points accepted in the 
C    scratchpad calculation. 
C 
    2   LLOX = NADDR(NINDEX)+2 
        LSIZEX = L(LLOX-1) 
        NINX = 0 
C 
C    The outer loop of the scratchpad calculation is a 
C    scan through the objects contiguous to NINDEX. 
C 
        DO 3  LKX = 1,LSIZEX 
        LLX = NADDR(NINDEX)+1+LKX 
        NMOD = L(LLX) 
C 
C    Nothing to do unless NMOD is a point. 
C 
        IF(NMOD) 3,701,4 
  701   STOP 716 
C 
C    If any points are already accepted in the scratchpad 
C    calculation, jump forward to the main part of the 
C    loop. 
C 
    4   IF(NINX.GT.0) GOTO 5 
C 
C    This section initialises the scratchpad on detection 
C    of the first point contiguous to NINDEX.  Pick up 
C    beginning and size of its contiguity list, and garbage 
C    if necessary, checking that enough space can be recovered 
C    for its extended list. 
C 
        LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+L(2)+2.LE.LTOP) GOTO 6 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+L(2)+2.LE.LTOP) GOTO 6 
        IFLAG = 6 
        RETURN 
C 
C    Find top of main list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and kill 
C    old list.  Copy old main list across, and add marker. 
C 
    6   LHI = LLO-1+LSIZE 
        L(LPTR) = NMOD 
        L(LPTR+1) = LSIZE+1+L(2) 
        NADDR(NMOD) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 7  LL = LLO,LHI 
        L(LPTR) = L(LL) 
    7   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    The scratchpad part of the extended list is the 
C    boundary list in anticlockwise order. 
C 
        LHIC = 2+L(2) 
        DO 8  LL = 3,LHIC 
        LC = LHIC+3-LL 
        L(LPTR) = L(LC) 
    8   LPTR = LPTR+1 
C 
C    Now form the extended lists for the effective 
C    constraints; the scratchpad list in each case consists 
C    of the two contiguous effective constraints separated 
C    by the point NMOD. 
C 
        DO 9  LC = 3,LHIC 
        J = -L(LC) 
C 
C    Pick up beginning and size of main list, and garbage if 
C    necessary, checking that enough space can be recovered. 
C 
        LLO = JADDR(J)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LTOP) GOTO 10 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(J)+2 
        IF(LPTR+LSIZE+5.LE.LTOP) GOTO 10 
        IFLAG = 6 
        RETURN 
C 
C    Find top of main list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and kill 
C    old list.  Copy old main list across, and add marker. 
C 
   10   LHI = LLO-1+LSIZE 
        L(LPTR) = -J 
        L(LPTR+1) = LSIZE+4 
        JADDR(J) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 11  LL = LLO,LHI 
        L(LPTR) = L(LL) 
   11   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    Add scratchpad part of extended list, dealing 
C    properly with cases LC.EQ.3 and LC.EQ.LHIC. 
C 
        LC1 = LC-1 
        IF(LC1.LT.3) LC1 = LHIC 
        L(LPTR) = L(LC1) 
        L(LPTR+1) = NMOD 
        LC1 = LC+1 
        IF(LC1.GT.LHIC) LC1 = 3 
        L(LPTR+2) = L(LC1) 
    9   LPTR = LPTR+3 
C 
C    Initialise NEAR to NMOD, reset NINX to 1, and if 
C    necessary reset NSTART; the scratchpad calculation for 
C    the first point contiguous to NINDEX is then complete. 
C 
        NEAR = NMOD 
        NINX = 1 
        IF(NSTART.EQ.NINDEX) NSTART = NMOD 
        GOTO 3 
C 
C    Points subsequent to the first are added into the 
C    scratchpad data-base as in subroutine ADDPT, but apart 
C    from the need to copy items of the main data-base at each 
C    step, the calculations are a little simpler.  This is 
C    because various acceptability checks have already been 
C    made and because there is no need for sophistication 
C    over locating the nearest neighbour.  First extract the 
C    coordinates of the new point. 
C 
    5   XMOD = PT(1,NMOD) 
        YMOD = PT(2,NMOD) 
C 
C    Find its nearest neighbour among the points already in 
C    the scratchpad calculation. 
C 
        XDIFF = XMOD-PT(1,NEAR) 
        YDIFF = YMOD-PT(2,NEAR) 
        D = XDIFF*XDIFF+YDIFF*YDIFF 
        LLOX = NADDR(NINDEX)+2 
        LKX1 = LKX-1 
        DO 12  LKXX = 1,LKX1 
        LLXX = LLOX-1+LKXX 
        NTRY = L(LLXX) 
        IF(NTRY) 12,702,13 
  702   STOP 717 
   13   XDIFF = XMOD-PT(1,NTRY) 
        YDIFF = YMOD-PT(2,NTRY) 
        DTRY = XDIFF*XDIFF+YDIFF*YDIFF 
        IF(DTRY.GE.D) GOTO 12 
        D = DTRY 
        NEAR = NTRY 
   12   CONTINUE 
C 
C    Initialise NNEW, the first item to be placed on the 
C    scratchpad contiguity list for NMOD; it is the neighbour 
C    of NMOD clockwise from NEAR.  Flag degeneracy if it 
C    occurs. 
C 
C    First pick up the beginning and end of the extended 
C    contiguity list for NEAR, and its coordinates. 
C 
        LLO = NADDR(NEAR)+2 
        LHI = LLO-1+L(LLO-1) 
        XNEAR = PT(1,NEAR) 
        YNEAR = PT(2,NEAR) 
C 
C    Initialise NNEW to zero, and locate the marker to give 
C    the beginning of the scratchpad contiguity list. 
C 
        NNEW = 0 
        DO 14  LL = LLO,LHI 
        N = L(LL) 
        IF(N)  14,15,14 
   14   CONTINUE 
        STOP 720 
   15   LLO = LL+1 
C 
C    Now find NNEW as in subroutine ADDPT. 
C 
        DO 16  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XMOD,YMOD,XNEAR,YNEAR,N,T, 
     1    IWIPE) 
        IF(IWIPE.LE.0) GOTO 16 
        IF(NNEW.EQ.0) GOTO 17 
        IF(TNEW-T) 16,18,17 
   18   IDGNEW = 1 
        IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 16 
   17   IDGNEW = 0 
        NNEW = N 
        TNEW = T 
   16   CONTINUE 
C 
C    Initialise NCURR to NEAR and LEFF to LTOP.  The 
C    scratchpad contiguity list for NMOD is constructed 
C    in a temporary position working down from LTOP. 
C 
        NCURR = NEAR 
        LEFF = LTOP 
C 
C    Enter the main inner loop of the scratchpad 
C    calculation, carried out to incorporate NMOD into 
C    the scratchpad tessellation.  First update NOLD and 
C    NCURR, and pass back the degeneracy flag. 
C 
   19   NOLD = NCURR 
        NCURR = NNEW 
        IDGCUR = IDGNEW 
C 
C    Is NCURR a constraint or a point? 
C 
        IF(NCURR) 20,704,21 
  704   STOP 721 
C 
C    Constraint.  (XCURR,YCURR) is the reflexion of NMOD 
C    in it. 
C 
   20   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XMOD+BJ*YMOD+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XMOD-AJ*DJ 
        YCURR = YMOD-BJ*DJ 
C 
C    Pick up beginning and size of extended contiguity list,	 
C    and garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = JADDR(JCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        IFLAG = 6 
        RETURN 
C 
C    Point, with coordinates (XCURR,YCURR). 
C 
   21   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
C 
C    Pick up beginning and size of extended contiguity list, 
C    and garbage if necessary, checking that enough space 
C    can be recovered. 
C 
        LLO = NADDR(NCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        IFLAG = 6 
        RETURN 
C 
C    Find the end of the old contiguity list for NCURR, and 
C    kill that list.  Save the value of LPTR as the address 
C    of the new contiguity list for NCURR.  Flag the new 
C    list as active.  Copy the old main list across, and the 
C    marker, and start the new scratchpad list with NMOD. 
C 
   22   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LCURR = LPTR 
        L(LPTR) = NCURR 
        LPTR = LPTR+2 
        DO 23  LL = LLO,LHI 
        N = L(LL) 
        L(LPTR) = N 
        LPTR = LPTR+1 
        IF(N.EQ.0) GOTO 24 
   23   CONTINUE 
        STOP 722 
   24   LLO = LL+1 
        L(LPTR) = NMOD 
        LPTR = LPTR+1 
C 
C    Contiguity is a symmetric relation: find NOLD in the 
C    scratchpad part of the old extended contiguity list 
C    for NCURR.  Resetting LLO has already located the start 
C    of the scratchpad. 
C 
        DO 25  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 26 
   25   CONTINUE 
        STOP 723 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for NMOD after NCURR; it is the neighbour of 
C    NMOD clockwise from NCURR.  Flag degeneracy if it occurs, 
C    and backspace LOLD if degeneracy occurred last time. 
C 
   26   NNEW = 0 
        IDGNEW = 0 
        LL = LOLD 
        IF(IDGCUR.EQ.0) GOTO 27 
        LOLD = LOLD-1 
        IF(LOLD.LT.LLO) LOLD = LHI 
   27   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XMOD,YMOD,XCURR,YCURR,N,T, 
     1    IWIPE) 
        IF(IWIPE.GT.0) GOTO 28 
        IF(NNEW) 29,27,29 
   28   IF(NNEW.EQ.0) GOTO 30 
        IF(TNEW-T) 29,31,30 
   31   IDGNEW = 1 
        GOTO 32 
   30   NNEW = N 
        TNEW = T 
        GOTO 27 
C 
C    Complete the scratchpad part of the new contiguity list 
C    for NCURR, and set the appropriate pointer to it from 
C    either JADDR or NADDR. 
C 
   29   L(LPTR) = NNEW 
        LPTR = LPTR+1 
   32   L(LPTR) = N 
        LPTR = LPTR +1 
        IF(LL.EQ.LOLD) GOTO 33 
        LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        GOTO 32 
   33   L(LCURR+1) = LPTR-LCURR-2 
        IF(NCURR) 34,707,35 
  707   STOP 724 
   34   JADDR(JCURR) = LCURR 
        GOTO 36 
   35   NADDR(NCURR) = LCURR 
C 
C    Add NCURR to the temporary scratchpad contiguity list 
C    for NMOD. 
C 
   36   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NEAR to the temporary scratchpad 
C    contiguity list for NMOD, then the main inner loop is 
C    completed. 
C 
        IF(NCURR.NE.NEAR) GOTO 19 
C 
C    Pick up the beginning and size of the old contiguity 
C    list for NMOD, and garbage if necessary, checking 
C    that enough space can be recovered for the extended 
C    list. 
C 
        LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+2.LE.LEFF) GOTO 37 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+2.LE.LEFF) GOTO 37 
        IFLAG = 6 
        RETURN 
C 
C    Find top of old list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and 
C    kill old list.  Copy old list across and add marker. 
C 
   37   LHI = LLO-1+LSIZE 
        L(LPTR) = NMOD 
        L(LPTR+1) = LSIZE+1+LTOP-LEFF 
        NADDR(NMOD) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 38  LL = LLO,LHI 
        L(LPTR) = L(LL) 
   38   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    Copy the scratchpad list down from its temporary 
C    position at the top of the heap. 
C 
        LEFF = LEFF+1 
        DO 39  LL = LEFF,LTOP 
        L(LPTR) = L(LL) 
   39   LPTR = LPTR+1 
C 
C    Increment the count of points accepted in the 
C    scratchpad tessellation; this completes the outer 
C    main loop of the scratchpad calculation. 
C 
        NINX = NINX+1 
    3   CONTINUE 
C 
C    The scratchpad calculation is complete; the information 
C    it provides now has to be combined with that in the 
C    main data base. 
C 
C    The first stage is to merge main and scratchpad 
C    lists for objects contiguous to NINDEX.  In order 
C    to deal properly with degeneracy, a record has to 
C    be maintained of the objects preceding and 
C    following the current one; the scan through the 
C    contiguity list is designed to facilitate this. 
C 
        LLX = NADDR(NINDEX)+1+LSIZEX 
        NMOD = L(LLX-1) 
        NNEXT = L(LLX) 
        DO 40  LKX = 1,LSIZEX 
        LLX = NADDR(NINDEX)+1+LKX 
        NLAST = NMOD 
        NMOD = NNEXT 
        NNEXT = L(LLX) 
C 
C    Is NMOD a constraint or a point? 
C 
        IF(NMOD) 41,710,42 
  710   STOP 725 
C 
C    Constraint.  Pick up beginning and size of contiguity 
C    list, and garbage if necessary, checking that enough 
C    space can be recovered for merged list. 
C 
   41   JMOD = -NMOD 
        LLO = JADDR(JMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JMOD)+2 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        IFLAG = 6 
        RETURN 
C 
C    Point.  Pick up beginning and size of contiguity 
C    list, and garbage if necessary, checking that enough 
C    space can be recovered for merged list. 
C 
   42   LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        IFLAG = 6 
        RETURN 
C 
C    Find top of old list.  Save value of LPTR as 
C    address of new list, and set back pointer.  Kill old 
C    list. 
C 
   43   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LMOD = LPTR 
        L(LPTR) = NMOD 
        LPTR = LPTR+2 
C 
C    Find NINDEX in the main part of the extended list, 
C    copying across as far as that point, but not 
C    including it.  Then space one past NINDEX. 
C 
        DO 44  LL = LLO,LHI 
        N = L(LL) 
        IF(N.EQ.NINDEX) GOTO 45 
        L(LPTR) = N 
   44   LPTR = LPTR+1 
        STOP 726 
   45   LL = LL+1 
C 
C    Find the marker, and space on one to the beginning 
C    of the scratchpad list. 
C 
        DO 46  LMARK = LL,LHI 
        IF(L(LMARK).EQ.0) GOTO 47 
   46   CONTINUE 
        STOP 727 
   47   LMARK = LMARK+1 
C 
C    Find NNEXT in the scratchpad list. 
C 
        DO 48  LL1 = LMARK,LHI 
        IF(L(LL1).EQ.NNEXT) GOTO 49 
   48   CONTINUE 
        STOP 730 
C 
C    Copy the scratchpad list into the merged list, starting 
C    after (or, in the case of degeneracy, at) NNEXT and 
C    stopping at NLAST.  A special case arises if NINDEX 
C    occurs at position LLO. 
C 
   49   IF(LL.NE.(LLO+1).AND.NNEXT.EQ.L(LPTR-1)) LPTR = LPTR-1 
   50   N = L(LL1) 
        L(LPTR) = N 
        LPTR = LPTR+1 
        LL1 = LL1+1 
        IF(LL1.GT.LHI) LL1 = LMARK 
        IF(N.NE.NLAST) GOTO 50 
C 
C    Copy the rest of the main list into the merged list, 
C    starting after (or, in the case of degeneracy, at) 
C    NLAST.  Deal with special case. 
C 
        LHI = LMARK-2 
        IF(LL.GT.LHI) LL = LLO 
        IF(L(LL).EQ.NLAST) LPTR = LPTR-1 
        IF(LL.EQ.LLO) GOTO 51 
        DO 52  LL1 = LL,LHI 
        L(LPTR) = L(LL1) 
   52   LPTR = LPTR+1 
        IF(L(LPTR-1).EQ.L(LMOD+2)) LPTR = LPTR-1 
C 
C    Set up the length of the merged list and the address 
C    pointer to it. 
C 
   51   L(LMOD+1) = LPTR-LMOD-2 
        IF(NMOD) 53,714,54 
  714   STOP 731 
   53   JADDR(JMOD) = LMOD 
        GOTO 40 
   54   NADDR(NMOD) = LMOD 
   40   CONTINUE 
C 
C    Merging is now complete.  The next step is to kill 
C    the scratchpad lists for constraints not contiguous 
C    to NINDEX.  This is done, exceptionally, in situ, by 
C    making the marker and scratchpad section look like 
C    a dead data item, and resetting the length to that for 
C    the main section.  Nothing need be done to those 
C    constraints contiguous to NINDEX, as they have 
C    already been dealt with in the merging operation. 
C 
        DO 55  J = 1,JCNS 
        LLO = JADDR(J)+2 
        LHI = LLO-1+L(LLO-1) 
        DO 56  LMARK = LLO,LHI 
        IF(L(LMARK).EQ.0) GOTO 57 
   56   CONTINUE 
        GOTO 55 
   57   L(LMARK+1) = LHI-LMARK-1 
        L(LLO-1) = LMARK-LLO 
   55   CONTINUE 
C 
C    Copy the contiguity list for NINDEX, garbaging if 
C    necessary and checking that enough space can be 
C    recovered for this list. 
C 
        LLOX = NADDR(NINDEX)+2 
        LSIZEX = L(LLOX-1) 
        IF(LPTR+LSIZEX+1.LE.LTOP) GOTO 58 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLOX = NADDR(NINDEX)+2 
        IF(LPTR+LSIZEX+1.LE.LTOP) GOTO 58 
        IFLAG = 6 
        RETURN 
   58   LHIX = LLOX-1+LSIZEX 
        LL = LPTR+2 
        DO 59  LLX = LLOX,LHIX 
        L(LL) = L(LLX) 
   59   LL = LL+1 
        L(LPTR) = NINDEX 
        L(LPTR+1) = LSIZEX 
C 
C    Kill the original version of the contiguity list 
C    for NINDEX, chain NINDEX back into the beginning 
C    of the free space in NADDR, decrement NPTSIN, and 
C    return with IFLAG set to zero. 
C 
        L(LLOX-2) = 0 
        NADDR(NINDEX) = -NFREE 
        NFREE = NINDEX 
        NPTSIN = NPTSIN-1 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYPT 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Finds the contiguities which the trial point (XPT,YPT) 
C    would have if it were added to the tessellation, but 
C    does not actualy add it. 
C 
C    The subroutine constructs, at the exit value of LPTR, 
C    a data item (which, being above LPTR, is unprotected) 
C    whose back pointer is unset and which gives the 
C    contiguity list which the trial point (XPT,YPT) 
C    would have if it were added to the tessellation. 
C    Existing contiguity lists are not modified, although 
C    data items may be relocated as a result of garbage 
C    collection.  As with ADDPT, a check is made that the 
C    trial point is adequately within the window.  Then 
C    the data structure is interrogated to determine 
C    whether there are any accepted points; if not, the 
C    contiguity list for the trial point is obtained simply 
C    by reversing the boundary list.  If there are some 
C    accepted points the (a) nearest neighbour to the trial 
C    point is found and a check is made that it is not too 
C    close to the trial point.  The contiguity list for the 
C    trial point is constructed in a temporary position at 
C    the top of the heap - it is not necessary to take 
C    account of degeneracies, since other contiguity lists 
C    are not modified.  On succesful completion of the 
C    construction, the list is copied down to its final 
C    position at LPTR, and return takes place with IFLAG 
C    set to zero.  Nonzero values of IFLAG indicate failure 
C    to construct the contiguity list for the trial point, 
C    as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are each interpreted in terms 
C    of the tolerance values EPSCN and EPSPT in the argument 
C    list.  If there are no accepted points, return with 
C    IFLAG set to 8 cannot, of course, occur. NINDEX is 
C    the index of the (a) nearest neighbour to the trial 
C    point if there are any accepted points. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Call subroutine SEEPT to check that the trial point 
C    is far enough inside the window; return if not, with 
C    NINDEX giving the index (flagged with a minus sign) 
C    of an unsatisfied constraint. 
C 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    If any points are accepted, jump forward to the main 
C    part of the routine 
C 
    1   IF(NPTSIN.GT.0) GOTO 2 
C 
C    There are no accepted points.  Check that there is 
C    enough space available in the heap for the contiguity 
C    list of the trial point. 
C 
        IF(LPTR+L(2)+1.LE.LTOP) GOTO 3 
        IFLAG = 6 
        RETURN 
C 
C    Construct the contiguity list for the trial point, and 
C    its length, as part of a data item at LPTR, leaving 
C    the back pointer of that data item unset. 
C 
    3   L(LPTR+1) = L(2) 
        LTRIAL = LPTR+1+L(2) 
        LHI = 2+L(2) 
        DO 4  LL = 3,LHI 
        L(LTRIAL) = L(LL) 
    4   LTRIAL = LTRIAL-1 
C 
C    Return with IFLAG set to zero. 
C 
        IFLAG = 0 
        RETURN 
C 
C    The main part of the routine deals with the case where 
C    there are some accepted points.  Call subroutine LOCPT 
C    to find the (a) nearest neighbour of the trial point. 
C    If it is too close, return.  The index of the nearest 
C    neighbour is set in NINDEX. 
C 
    2   CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 5 
        RETURN 
C 
C    Initialise NNEW, the first item to be placed on the 
C    contiguity list for the trial point; it is the 
C    neighbour of the trial point clockwise from its 
C    nearest neighbour. 
C 
    5   NCURR = NINDEX 
        LLO = NADDR(NCURR)+2 
        LHI = LLO-1+L(LLO-1) 
        XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
        NNEW = 0 
        DO 6  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.LE.0) GOTO 6 
        IF(NNEW.EQ.0) GOTO 7 
        IF(TNEW-T) 6,51,7 
   51   IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 6 
    7   NNEW = N 
        TNEW = T 
    6   CONTINUE 
C 
C    Initialise LEFF to LTOP; the contiguity list for 
C    the trial point is constructed in a temporary position 
C    working down from LTOP. 
C 
        LEFF = LTOP 
C 
C    Enter the main loop.  First update NOLD and NCURR. 
C 
    8   NOLD = NCURR 
        NCURR = NNEW 
C 
C    Check that there is enough space to add another point 
C    to the temporary list for the trial point, and garbage 
C    if necessary, checking that enough space can be 
C    recovered. 
C 
        IF(LPTR+2.LE.LEFF) GOTO 9 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        IF(LPTR+2.LE.LEFF) GOTO 9 
        IFLAG = 6 
        RETURN 
C 
C     Is NCURR a constraint or a point? 
C 
    9   IF(NCURR) 10,701,11 
  701   STOP 732 
C 
C    Constraint: (XCURR,YCURR) is the reflexion of the 
C    trial point in it.  Pick up beginning of contiguity 
C    list. 
C 
   10   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XPT+BJ*YPT+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XPT-AJ*DJ 
        YCURR = YPT-BJ*DJ 
        LLO = JADDR(JCURR)+2 
        GOTO 12 
C 
C    Point, with coordinates (XCURR,YCURR): pick up 
C    beginning of contiguity list. 
C 
   11   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
        LLO = NADDR(NCURR)+2 
C 
C    Find the end of the contiguity list for NCURR and then, 
C    using symmetry of the contiguity relation, find NOLD 
C    in it. 
C 
   12   LHI = LLO-1+L(LLO-1) 
        DO 13  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 14 
   13   CONTINUE 
        STOP 733 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for the trial point after NCURR; it is the 
C    neighbour of the trial point clockwise from NCURR. 
C 
   14   NNEW = 0 
        LL = LOLD 
   15   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.GT.0) GOTO 16 
        IF(NNEW) 17,15,17 
   16   IF(NNEW.EQ.0) GOTO 18 
        IF(TNEW.LE.T) GO TO 17 
   18   NNEW = N 
        TNEW = T 
        GOTO 15 
C 
C    Add NCURR to the contiguity list for the trial point. 
C 
   17   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NINDEX to the contiguity list 
C    for the trial point, then the main loop is completed. 
C 
        IF(NCURR.NE.NINDEX) GOTO 8 
C 
C    Move the contiguity list for the trial point down to 
C    LPTR, as part of a data item at this location, and 
C    set up its size. 
C 
        L(LPTR+1) = LTOP-LEFF 
        LTRIAL = LPTR+2 
        LEFF = LEFF+1 
        DO 19  LL = LEFF,LTOP 
        L(LTRIAL) = L(LL) 
   19   LTRIAL = LTRIAL+1 
C 
C    The construction of the contiguity list for the trial 
C    point is complete.  Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILE4M 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TILE4M(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB) 
C 
C    A master routine for constructing a tessellation within an 
C    arbitrary window.  This routine gives access to the full TILE4 
C    database, and in new applications should be used in preference 
C    to subroutine TILE, which is supplied for compatibility with the 
C    routine of that name in TILE3.  Subroutine TILE4M can be followed 
C    by calls to ADDPT,SUBPT to manipulate the database, but this is 
C    not possible after a call to TILE.  Note that it is legitimate to 
C    call TILE4M with a temporary small value of NPTS in such cases - 
C    see below for details. 
C 
C    The user call to TILE4M must supply all the arguments.  The 
C    values of JCNS,NPTS,LTOP,EPSCN,EPSPT may be given as constants in 
C    the call, but all other arguments must be array or simple variable 
C    names, as appropriate.  The arguments are as follows. 
C 
C      CN      REAL array of dimension (3,JCNS).  User supplies 
C              constraints according to the convention 
C 
C                      CN(1,J)*x + CN(2,J)*y + CN(3,J) .LT. 0.0 
C 
C              for J = 1,...,JCNS. 
C 
C      JCNS    INTEGER variable in which user supplies number of 
C              constraints and also second dimension of CN and 
C              dimension of JADDR. 
C 
C      JADDR   INTEGER array of dimension JCNS supplied by user as 
C              uninitialised space for part of the TILE4 database. 
C              Contains heap addresses of constraint data items 
C              on return. 
C 
C      PT      REAL array of dimension (2,NPTS).  User supplies 
C              point coordinates according to the convention 
C 
C                      (PT(1,N),PT(2,N))  is  (x,y) 
C 
C              for N = 1,...,NPTS. 
C 
C      NPTS    INTEGER variable in which user supplies number of 
C              points and also second dimension of PT and 
C              dimension of NADDR.  In applications where a call to 
C              TILE4M is to be followed by subsequent manipulation 
C              of the database, it may be desirable to set up PT(2, ) 
C              and NADDR( ) to a larger true dimension in the 
C              calling program than the dimension actually passed 
C              to TILE4M.  This is quite legitimate provided that 
C              the call to TILE4M is followed by a call to 
C              subroutine EXTEND to initialise the remainder of 
C              NADDR before calls to ADDPT,SUBPT are made with 
C              the larger value of NPTS in the calling program. 
C 
C      NADDR   INTEGER array of dimension NPTS supplied by user 
C              as uninitialised space for part of TILE4 database. 
C              Contains heap addresses of point data items on 
C              return. 
C 
C      NFREE   INTEGER variable not initialised by user.  Returns a 
C              pointer to the first free location in PT,NADDR, or 
C              returns NPTS+1 if there are none. 
C 
C      NSTART  INTEGER variable holding index of starting-point for 
C              nearest-neighbour walk.  User initialisation controls 
C              updating as follows.  If NSTART is positive or zero on 
C              entry, then it is updated occasionally to lie near the 
C              centroid of the accepted points.  If it is negative on 
C              entry, then each nearest-neighbour walk commences from 
C              the previous accepted point.  Use 0,-1 as initial values 
C              by convention.  Return value is the index of the next 
C              starting-point which would have been used. 
C              The 0 option should be used when the order of the points 
C              in the array PT bears little relation to their actual 
C              position in the plane.  The -1 option should be used 
C              when the points in nearby locations in the array PT 
C              are likely to have nearby positions in the plane. 
C 
C      NPTSIN  INTEGER variable not initialised by user.  Returns the 
C              number of currently accepted points. 
C 
C      L       INTEGER (or INTEGER*2) array of dimension LTOP. 
C              Supplied by user as uninitialised space for part of the 
C              TILE4 database.  Holds a heap of data items on return. 
C 
C      LTOP    INTEGER variable in which user supplies dimension 
C              of L.  A value of 9*NPTS, or 1000 if that is greater, 
C              is recommended, although a larger value may be 
C              needed if the database is to be manipulated 
C              subsequently by ADDPT,SUBPT. 
C 
C      LPTR    INTEGER variable not initialised by user.  Returns a 
C              pointer to the first free location in L. 
C 
C      LBASE   INTEGER variable not initialised by user.  Returns a 
C              pointer to the base of the working space in L. 
C 
C      EPSCN   REAL variable in which user supplies tolerance 
C              value for "inside" window.  Negative sign is treated 
C              as zero. 
C 
C      EPSPT   REAL variable in which user supplies tolerance 
C              value for "identical" points.  Negative sign is treated 
C              as zero. 
C 
C      IFLAG   INTEGER variable not initialised by user.  Returns 
C              value zero on completely successful exit from TILE4M. 
C              Nonzero values indicate complete or partial failure 
C              of TILE4M as follows. 
C 
C              1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C                 for some value of J.  TILE4M fails. 
C 
C              2  Unbounded window - constraints inadequate. 
C                 TILE4M fails. 
C 
C              3  Empty window - constraints inconsistent. 
C                 TILE4M fails. 
C 
C              4  Not all points "inside" window.  TILE4M succeeds 
C                 by rejecting such points. 
C 
C              5  Some points "identical".  TILE4M succeeds by 
C                 rejecting points which duplicate an already 
C                 accepted point.  If both errors 4 and 5 occur, 
C                 error 4 is reported. 
C 
C              6  Heap overflow - LTOP too small.  TILE4M fails. 
C 
C      IGARB   INTEGER variable not initialised by user.  Returns 
C              number of garbage collections necessarily carried out 
C              in TILE4M, that is, excluding the final "cosmetic" 
C              garbage collection carried out immediately before 
C              control is returned to the calling program. 
C 
C    On return from TILE4M, the TILE4 database is held as explained 
C    above.  Some of the values returned are of no direct interest 
C    to the user, and are returned simply so that they can be passed 
C    on to other routines.  The user resets items in the database 
C    at his own risk, and can easily produce chaos by attempting to 
C    do so other than by calling routines supplied as part of the 
C    TILE4 package.  Users intending to access the database directly, 
C    rather than via the supplied interrogation routines, will need 
C    to understand its structure, which is as follows.  Each 
C    effective constraint, and each accepted point, has a contiguity 
C    list of neighbouring points and constraints held as part of a 
C    data item in the heap L.  The data item for constraint J is 
C    at LL = JADDR(J), and L(LL) is set to -J.  The data item for 
C    point N is at LL = NADDR(N), and L(LL) is set to N.  Ineffective 
C    constraints have JADDR(J) set to zero.  Rejected points have 
C    NADDR(N) set to a negative value.  These negative values when 
C    stripped of the sign form a chain entered at NFREE and terminated 
C    by (-)(NPTS+1).  A data item in L consists of the back pointer 
C    already described, then a location holding the length of the 
C    succeeding list, then the contiguity list itself, in as many 
C    locations (always at least 3) as are required, with points 
C    referenced by their index and constraints by their index with 
C    a minus sign flagging them.  The order of objects within 
C    the lists is anticlockwise cyclic.  Storage in L is administered 
C    as a heap, that is to say, data items are not modified 
C    in situ but are replaced by new data items at the base of the 
C    free space.  Abandoned data items are flagged as dead by 
C    having their back-pointers set to zero, and subroutine GARBAJ 
C    returns the space occupied by such items to the top 
C    of the heap as free space when required.  At the very bottom 
C    of the heap is a special data item protected against 
C    removal by GARBAJ.  This has the following form.  L(1) 
C    is unset.  L(2) is the length of the succeeding list. 
C    Succeeding entries in L, as many as are needed, hold a 
C    list of effective constraints (each as usual flagged 
C    with a minus sign) in clockwise cyclic order round the 
C    window boundary.  LBASE is the first location above this data 
C    item. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Hard error for nonpositive dimensions. 
C 
        IF(JCNS.LE.0.OR.NPTS.LE.0.OR.LTOP.LE.0) STOP 734 
C 
C    Set up the window and check its validity. 
C 
        CALL WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Save information from NSTART, then complete initialisation. 
C 
    1   IRESET = 0 
        IF(NSTART.LT.0) IRESET = -1 
        CALL CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN,LPTR,LBASE,IGARB) 
        IFSAVE = 0 
        NRESET = 2 
C 
C    Scan the points one by one, attempting to insert each. 
C 
        DO 2  NPT = 1,NPTS 
        XPT = PT(1,NPT) 
        YPT = PT(2,NPT) 
        CALL ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
        IF(IFLAG.EQ.7) STOP 735 
        IF(IFLAG.EQ.8) STOP 736 
        IF(IFLAG.NE.6) GOTO 3 
        RETURN 
C 
C    Nothing disastrous has happened.  But something mildly nasty 
C    may have done. 
C 
    3   IF(IFLAG.EQ.0) GOTO 4 
        IF(IFLAG.EQ.4) IFSAVE = 4 
        IF(IFLAG.EQ.5.AND.IFSAVE.EQ.0) IFSAVE = 5 
        NSAVE = -NADDR(NPTS) 
        NADDR(NPTS) = -NPT 
        NFREE = -NADDR(NPT) 
        NADDR(NPT) = -NSAVE 
        GOTO 2 
C 
C    The point has been accepted.  Deal with NSTART. 
C 
    4   IF(IRESET.EQ.-1) GOTO 5 
        IF(NPTSIN.LT.NRESET) GOTO 2 
        NRESET = 2*NRESET 
        CALL MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP,IFLAG, 
     1    XBAR,YBAR,NINDEX) 
        IF(IFLAG.EQ.4) STOP 737 
        IF(IFLAG.EQ.8) STOP 740 
    5   NSTART = NINDEX 
    2   CONTINUE 
C 
C    Final garbage collection, set IFLAG, and return. 
C 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IFLAG = IFSAVE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE EXTEND 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE EXTEND(NADDR,NPTS,NSUB) 
C 
C    Extends the free space chain in NADDR from NSUB to NPTS. 
C 
        DIMENSION NADDR(NPTS) 
        IF(NSUB.GT.NPTS) STOP 741 
        IF(NSUB.EQ.NPTS) RETURN 
        NS1 = NSUB+1 
        DO 1  N = NS1,NPTS 
    1   NADDR(N) = -(N+1) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    INTERPOLATION ROUTINES LAST MODIFIED 13 MARCH 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RECT1G 4.1 DATED 13 MARCH 1981 
C 
        SUBROUTINE RECT1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,SBAREA, 
     2    DELSBA,PTOFF,KNBMAX,VAL,GRAD,Z,ZX,ZY,MM,M,N, 
     3    XMIN,XMAX,YMIN,YMAX,ZLO,ZHI) 
C 
C    Evaluates the C1 natural neighbour interpolant and its gradient 
C    at all points on a rectangular grid. 
C 
C    User provides the TILE4 database; working arrays SBAREA(KNBMAX), 
C    DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), and KNBMAX to dimension them; 
C    arrays VAL(NPTS),GRAD(2,NPTS) holding data values and gradients at 
C    data sites, with GRAD usually set by a call to GRLD; 
C    arrays Z(MM,N),ZX(MM,N),ZY(MM,N) in which the routine 
C    returns value and components of gradient at the grid points, with 
C    MM as the true first dimension of these arrays, and M,N as the 
C    numbers of grid points in the X and Y directions respectively; 
C    XMIN,XMAX,YMIN,YMAX to specify the position of the grid; and ZLO, 
C    ZHI in which the routine returns the minimum and maximum values of 
C    the interpolant on the grid.  The position and spacing of the grid 
C    is defined as follows.  The rectangle ÝXMIN,XMAX¨*ÝYMIN,YMAX¨ is 
C    divided into M*N small rectangles each (XMAX-XMIN)/FLOAT(M) by 
C    (YMAX-YMIN)/FLOAT(N); the interpolant is evaluated at the centres 
C    of these.  Thus grid point (MSCAN,NSCAN) is at 
C 
C       (XMIN+(FLOAT(MSCAN)-0.5)*(XMAX-XMIN)/FLOAT(M), 
C          YMIN+(FLOAT(NSCAN)-0.5)*(YMAX-YMIN)/FLOAT(N)) 
C 
C    The grid need not lie strictly within the window; at grid points 
C    outside the window, ZX and ZY are set conventionally to zero, and 
C    Z is set to 2.0*ZLO-ZHI to allow such points to be detected without 
C    the need for an array of flags.  IFLAG is zero on successful 
C    return.  Nonzero values indicate error or nonstandard returns as 
C    follows. 
C 
C       4  Some grid point outside window 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    Note that heap overflow is not here a catastrophic error.  A return 
C    value of 4 for IFLAG is not an error in this context - it simply 
C    alerts the user to the fact that his grid is partly outside the 
C    window.  If grid points ever duplicate data sites, a value of 5 for 
C    IFLAG is generated internally by NNBR1G, but this again is not an 
C    error in this context; value and gradient are set correctly, and 
C    the flag is not passed back to the user. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS), 
     2    GRAD(2,NPTS),Z(MM,N),ZX(MM,N),ZY(MM,N) 
C 
C    Initialise values. 
C 
        ZLO = 1.0 
        ZHI = 0.0 
        IFLAG = 0 
        NS2 = NSTART 
        DELTAX = (XMAX-XMIN)/FLOAT(M) 
        DELTAY = (YMAX-YMIN)/FLOAT(N) 
        X0 = XMIN-0.5*DELTAX 
        Y0 = YMIN-0.5*DELTAY 
        XSCAN = X0 
C 
C     Scan the rectangular grid. 
C 
        DO 1  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        IF(NS2.EQ.0) NS2 = NSTART 
        NS1 = NS2 
        NS2 = 0 
        DO 1  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL NNBR1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NS1,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFL,IGARB,XSCAN,YSCAN,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL,GRAD, 
     3    Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
C 
C    Check IFL and take appropriate action. 
C 
        IF(IFL.EQ.0.OR.IFL.EQ.5) GOTO 2 
        IFLAG = IFL 
        IF(IFL.EQ.4) GOTO 1 
        RETURN 
C 
C    We are dealing with a grid point inside the window. 
C 
    2   Z(MSCAN,NSCAN) = Z1 
        ZX(MSCAN,NSCAN) = ZX1 
        ZY(MSCAN,NSCAN) = ZY1 
        IF(ZLO.LE.ZHI) GOTO 3 
        ZLO = Z1 
        ZHI = Z1 
        GOTO 4 
    3   IF(ZLO.GT.Z1) ZLO = Z1 
        IF(ZHI.LT.Z1) ZHI = Z1 
    4   IF(NS2.EQ.0) NS2 = NINDEX 
        NS1 = NINDEX 
    1   CONTINUE 
C 
C    Deal with the flagging of points outside the window if any. 
C 
        IF(IFLAG.EQ.0) RETURN 
        ZVAL = 2.0*ZLO-ZHI 
        XSCAN = X0 
        DO 5  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        DO 5  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFL,XSCAN,YSCAN,JINDEX,TSTVAL) 
        IF(IFL.EQ.0) GOTO 5 
        Z(MSCAN,NSCAN) = ZVAL 
        ZX(MSCAN,NSCAN) = 0.0 
        ZY(MSCAN,NSCAN) = 0.0 
    5   CONTINUE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RECTHG 4.1 DATED 13 MARCH 1981 
C 
        SUBROUTINE RECTHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,SBAREA, 
     2    DELSBA,PTOFF,KNBMAX,VAL,Z,ZX,ZY,MM,M,N, 
     3    XMIN,XMAX,YMIN,YMAX,ZLO,ZHI) 
C 
C    Evaluates the C1 natural neighbour interpolant and its gradient 
C    at all points on a rectangular grid, with the gradients 
C    at the data sites all forced to be zero and GRAD 
C    accordingly not required.  Other comments are as for 
C    RECT1G. 
C 
C    User provides the TILE4 database; working arrays SBAREA(KNBMAX), 
C    DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), and KNBMAX to dimension them; 
C    array VAL(NPTS) holding values at data sites; 
C    arrays Z(MM,N),ZX(MM,N),ZY(MM,N) in which the routine 
C    returns value and components of gradient at the grid points, with 
C    MM as the true first dimension of these arrays, and M,N as the 
C    numbers of grid points in the X and Y directions respectively; 
C    XMIN,XMAX,YMIN,YMAX to specify the position of the grid; and ZLO, 
C    ZHI in which the routine returns the minimum and maximum values of 
C    the interpolant on the grid.  The position and spacing of the grid 
C    is defined as follows.  The rectangle ÝXMIN,XMAX¨*ÝYMIN,YMAX¨ is 
C    divided into M*N small rectangles each (XMAX-XMIN)/FLOAT(M) by 
C    (YMAX-YMIN)/FLOAT(N); the interpolant is evaluated at the centres 
C    of these.  Thus grid point (MSCAN,NSCAN) is at 
C 
C       (XMIN+(FLOAT(MSCAN)-0.5)*(XMAX-XMIN)/FLOAT(M), 
C          YMIN+(FLOAT(NSCAN)-0.5)*(YMAX-YMIN)/FLOAT(N)) 
C 
C    The grid need not lie strictly within the window; at grid points 
C    outside the window, ZX and ZY are set conventionally to zero, and 
C    Z is set to 2.0*ZLO-ZHI to allow such points to be detected without 
C    the need for an array of flags.  IFLAG is zero on successful 
C    return.  Nonzero values indicate error or nonstandard returns as 
C    follows. 
C 
C       4  Some grid point outside window 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    Note that heap overflow is not here a catastrophic error.  A return 
C    value of 4 for IFLAG is not an error in this context - it simply 
C    alerts the user to the fact that his grid is partly outside the 
C    window.  If grid points ever duplicate data sites, a value of 5 for 
C    IFLAG is generated internally by NNBRHG, but this again is not an 
C    error in this context; value and gradient are set correctly, and 
C    the flag is not passed back to the user. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS), 
     2    Z(MM,N),ZX(MM,N),ZY(MM,N) 
C 
C    Initialise values. 
C 
        ZLO = 1.0 
        ZHI = 0.0 
        IFLAG = 0 
        NS2 = NSTART 
        DELTAX = (XMAX-XMIN)/FLOAT(M) 
        DELTAY = (YMAX-YMIN)/FLOAT(N) 
        X0 = XMIN-0.5*DELTAX 
        Y0 = YMIN-0.5*DELTAY 
        XSCAN = X0 
C 
C     Scan the rectangular grid. 
C 
        DO 1  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        IF(NS2.EQ.0) NS2 = NSTART 
        NS1 = NS2 
        NS2 = 0 
        DO 1  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL NNBRHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NS1,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFL,IGARB,XSCAN,YSCAN,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL, 
     3    Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
C 
C    Check IFL and take appropriate action. 
C 
        IF(IFL.EQ.0.OR.IFL.EQ.5) GOTO 2 
        IFLAG = IFL 
        IF(IFL.EQ.4) GOTO 1 
        RETURN 
C 
C    We are dealing with a grid point inside the window. 
C 
    2   Z(MSCAN,NSCAN) = Z1 
        ZX(MSCAN,NSCAN) = ZX1 
        ZY(MSCAN,NSCAN) = ZY1 
        IF(ZLO.LE.ZHI) GOTO 3 
        ZLO = Z1 
        ZHI = Z1 
        GOTO 4 
    3   IF(ZLO.GT.Z1) ZLO = Z1 
        IF(ZHI.LT.Z1) ZHI = Z1 
    4   IF(NS2.EQ.0) NS2 = NINDEX 
        NS1 = NINDEX 
    1   CONTINUE 
C 
C    Deal with the flagging of points outside the window if any. 
C 
        IF(IFLAG.EQ.0) RETURN 
        ZVAL = 2.0*ZLO-ZHI 
        XSCAN = X0 
        DO 5  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        DO 5  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFL,XSCAN,YSCAN,JINDEX,TSTVAL) 
        IF(IFL.EQ.0) GOTO 5 
        Z(MSCAN,NSCAN) = ZVAL 
        ZX(MSCAN,NSCAN) = 0.0 
        ZY(MSCAN,NSCAN) = 0.0 
    5   CONTINUE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR1G 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE NNBR1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,GRAD,Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
C 
C    Calculates the C0 and C1 natural neighbour interpolants and 
C    their gradients at the point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolants are wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) as working space, 
C    and KNBMAX to dimension them, arrays VAL(NPTS),GRAD(2,NPTS) 
C    holding the data values and gradients at the data sites.  GRAD 
C    will usually have been loaded by a call to GRLD.  Z0 returns the 
C    C0 natural neighbour interpolant, (ZX0,ZY0) its gradient; Z1 
C    returns the C1 natural neighbour interpolant, (ZX1,ZY1) its 
C    gradient.  IFLAG is zero on successful return.  Nonzero values 
C    indicate error or nonstandard returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    "Outside" and "duplicate" are interpreted in terms of the 
C    tolerance values EPSCN,EPSPT.  If duplication occurs, NINDEX 
C    returns the index of the duplicated point.  Z0 returns the 
C    value at NINDEX, and (ZX0,ZY0) conventionally return (0.0,0.0). 
C    (ZX1,ZY1) return the gradient at NINDEX, and Z1 returns the 
C    value at (XPT,YPT) on the plane through Z0 at NINDEX with 
C    gradient (ZX1,ZY1).  Note that heap overflow is not here a 
C    catastrophic error.  ICASE returns the edge effect indicator 
C    from the call to TRYSBG, that is, 5 at the boundary and 0 
C    elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), 
     2    VAL(NPTS),GRAD(2,NPTS) 
C 
C    Call TRYSBG to calculate subareas and their gradients. 
C 
        CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   UPT = XPT-PT(1,NINDEX) 
        VPT = YPT-PT(2,NINDEX) 
        Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        ZX1 = GRAD(1,NINDEX) 
        ZY1 = GRAD(2,NINDEX) 
        Z1 = Z0+ZX1*UPT+ZY1*VPT 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        S1X = 0.0 
        S1Y = 0.0 
        S2X = 0.0 
        S2Y = 0.0 
        SM1X = 0.0 
        SM1Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
        TM1X = 0.0 
        TM1Y = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        GX = 0.0 
        GY = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
        GXX = 0.0 
        GXY = 0.0 
        GYX = 0.0 
        GYY = 0.0 
        EDGEXX = 0.0 
        EDGEXY = 0.0 
        EDGEYX = 0.0 
        EDGEYY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up and calculate values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT*UPT+VPT*VPT 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
        GRX = GRAD(1,N) 
        GRY = GRAD(2,N) 
        ZETA = ZPT+GRX*UPT+GRY*VPT 
C 
C    Accumulate main values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZETA/DPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        S1X = S1X+SX*DPT+S*UPT/DPT 
        S1Y = S1Y+SY*DPT+S*VPT/DPT 
        S2X = S2X+SX*DSQ 
        S2Y = S2Y+SY*DSQ 
        SM1X = SM1X+(SX-S*UPT/DSQ)/DPT 
        SM1Y = SM1Y+(SY-S*VPT/DSQ)/DPT 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
        TM1X = TM1X+(SX*ZETA+S*(GRX-ZETA*UPT/DSQ))/DPT 
        TM1Y = TM1Y+(SY*ZETA+S*(GRY-ZETA*VPT/DSQ))/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        GX = GX+S*GRX 
        GY = GY+S*GRY 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
        GXX = GXX+SX*GRX 
        GXY = GXY+SY*GRX 
        GYX = GYX+SX*GRY 
        GYY = GYY+SY*GRY 
        EDGEXX = EDGEXX+SX*UPT+S 
        EDGEXY = EDGEXY+SY*UPT 
        EDGEYX = EDGEYX+SX*VPT 
        EDGEYY = EDGEYY+SY*VPT+S 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        IF(ICASE.EQ.0) GOTO 5 
        S2X = S2X+2.0*EDGEX 
        S2Y = S2Y+2.0*EDGEY 
        CORR = (GX*EDGEX+GY*EDGEY)/S0 
        T0 = T0+CORR 
        T0X = T0X+(GXX*EDGEX+GX*EDGEXX+GYX*EDGEY+GY*EDGEYX-CORR*S0X)/S0 
        T0Y = T0Y+(GXY*EDGEX+GX*EDGEXY+GYY*EDGEY+GY*EDGEYY-CORR*S0Y)/S0 
    5   Z1DEN = S1*S0+S2*SM1 
        Z1 = (S1*T0+S2*TM1)/Z1DEN 
        ZX1 = S1X*T0+S1*T0X+S2X*TM1+S2*TM1X 
        ZX1 = (ZX1-Z1*(S1X*S0+S1*S0X+S2X*SM1+S2*SM1X))/Z1DEN 
        ZY1 = S1Y*T0+S1*T0Y+S2Y*TM1+S2*TM1Y 
        ZY1 = (ZY1-Z1*(S1Y*S0+S1*S0Y+S2Y*SM1+S2*SM1Y))/Z1DEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 DATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBRHG 4.1 DATED 10 MARCH 1981 
C 
        SUBROUTINE NNBRHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,Z0,ZX0,ZY0,ZH,ZXH,ZYH) 
C 
C    Calculates the C0 natural neighbour interpolant and its gradient, 
C    and also the C1 natural neighbour interpolant and its gradient with 
C    all data site gradients forced to zero, at the point (XPT,YPT). 
C    This routine thus produces the same effect as a call to NNBR1G 
C    with zeroes entered into GRAD, but the array GRAD does not need 
C    to be passed to it and the calculation is more efficient.  The 
C    C1 interpolant value is returned as ZH, its gradient as (ZXH,ZYH). 
C    Other details are as for NNBR1G. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 
C 
C    Call TRYSBG to calculate subareas and their gradients. 
C 
       CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        ZH = Z0 
        ZXH = 0.0 
        ZYH = 0.0 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        S1X = 0.0 
        S1Y = 0.0 
        S2X = 0.0 
        S2Y = 0.0 
        SM1X = 0.0 
        SM1Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
        TM1X = 0.0 
        TM1Y = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        EDGEX = 0.0 
        EDGEY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT**2+VPT**2 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZPT/DPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        S1X = S1X+SX*DPT+S*UPT/DPT 
        S1Y = S1Y+SY*DPT+S*VPT/DPT 
        S2X = S2X+SX*DSQ 
        S2Y = S2Y+SY*DSQ 
        SM1X = SM1X+(SX-S*UPT/DSQ)/DPT 
        SM1Y = SM1Y+(SY-S*VPT/DSQ)/DPT 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
        TM1X = TM1X+(SX-S*UPT/DSQ)*ZPT/DPT 
        TM1Y = TM1Y+(SY-S*VPT/DSQ)*ZPT/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        IF(ICASE.EQ.0) GOTO 5 
        S2X = S2X+2.0*EDGEX 
        S2Y = S2Y+2.0*EDGEY 
    5   ZHDEN = S1*S0+S2*SM1 
        ZH = (S1*T0+S2*TM1)/ZHDEN 
        ZXH = S1X*T0+S1*T0X+S2X*TM1+S2*TM1X 
        ZXH = (ZXH-ZH*(S1X*S0+S1*S0X+S2X*SM1+S2*SM1X))/ZHDEN 
        ZYH = S1Y*T0+S1*T0Y+S2Y*TM1+S2*TM1Y 
        ZYH = (ZYH-ZH*(S1Y*S0+S1*S0Y+S2Y*SM1+S2*SM1Y))/ZHDEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR1 4.2 DATED 13 MARCH 1981 
C 
        SUBROUTINE NNBR1(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,GRAD,Z0,Z1) 
C 
C    Calculates the C0 and C1 natural neighbour interpolants at the 
C    point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolants are wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX) as working space, and KNBMAX to 
C    dimension them, arrays VAL(NPTS),GRAD(2,NPTS) holding the data 
C    values and gradients at the data sites.  GRAD will usually have 
C    been loaded by a call to GRLD.  Z0 returns the C0 natural 
C    neighbour interpolant, Z1 the C1 natural neighbour interpolant. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    or nonstandard returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    "Outside" and "duplicate" are interpreted in terms of the 
C    tolerance values EPSCN,EPSPT.  If duplication occurs, NINDEX 
C    returns the index of the duplicated point, and Z0,Z1 return 
C    correctly.  Note that heap overflow is not here a catastrophic 
C    error.  ICASE returns the edge effect indicator value from the 
C    call to TRYSBA, that is, 5 at the boundary and 0 elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS),GRAD(2,NPTS) 
C 
C    Call TRYSBA to calculate subareas. 
C 
        CALL TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is zero or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   UPT = XPT-PT(1,NINDEX) 
        VPT = YPT-PT(2,NINDEX) 
        Z0 = VAL(NINDEX) 
        Z1 = Z0+GRAD(1,NINDEX)*UPT+GRAD(2,NINDEX)*VPT 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        GX = 0.0 
        GY = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up and calculate values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT**2+VPT**2 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        ZPT = VAL(N) 
        GRX = GRAD(1,N) 
        GRY = GRAD(2,N) 
        ZETA = ZPT+GRX*UPT+GRY*VPT 
C 
C    Accumulate main values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZETA/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        GX = GX+S*GRX 
        GY = GY+S*GRY 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        IF(ICASE.EQ.0) GOTO 5 
        CORR = (GX*EDGEX+GY*EDGEY)/S0 
        T0 = T0+CORR 
    5   Z1DEN = S1*S0+S2*SM1 
        Z1 = (S1*T0+S2*TM1)/Z1DEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR0G 4.3 DATED 9 MARCH 1981 
C 
        SUBROUTINE NNBR0G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL, 
     3    Z0,ZX0,ZY0) 
C 
C    Calculates the C0 natural neighbour interpolant and its gradient 
C    at the point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolant is wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) as working space, 
C    and KNBMAX to dimension them, and an array VAL(NPTS) holding the 
C    data values at the data sites.  Z0 returns the C0 natural 
C    neighbour interpolant, (ZX0,ZY0) its gradient.  IFLAG is zero on 
C    successful return.  Nonzero values indicate error or nonstandard 
C    returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If duplication occurs, NINDEX returns the 
C    index of the duplicated point.  Z0 returns the value at NINDEX. 
C    (ZX0,ZY0) conventionally return (0.0,0.0).  Note that heap overflow 
C    is not here a catastrophic error.  ICASE returns the edge effect 
C    indicator flag from the call to TRYSBG, that is, 5 at the boundary 
C    and 0 elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 
C 
C    Call TRYSBG to calculate subareas and gradients. 
C 
        CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators.  There is no need to worry about the edge case, 
C    because we do not try to compensate in the C0 natural neighbour 
C    interpolant. 
C 
    1   S0 = 0.0 
        T0 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
        LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 3 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        T0 = T0+S*ZPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
    3   CONTINUE 
C 
C    Calculate values for return. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR0 4.2 DATED 12 MARCH 1981 
C 
        SUBROUTINE NNBR0(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE,VAL,Z0
     3   ,VSHAPE) 
C 
C    Calculates the C0 natural neighbour interpolant at the point 
C    (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolant is wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX) as working space, and KNBMAX to 
C    dimension them, and an array VAL(NPTS) holding the data values at 
C    the data sites.  Z0 returns the C0 natural neighbour interpolant. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    or nonstandard returns as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If duplication occurs, NINDEX returns the 
C    index of the duplicated point, and Z0 returns the value at NINDEX. 
C    Note that heap overflow is not here a catastrophic error.  ICASE 
C    returns the edge effect indicator flag from the call to TRYSBA, 
C    that is, 5 at the boundary and 0 elsewhere. 
C 

!/
!This subroutine NNINT was modified by Frank Tsai for the main program Lithology.90
!/

C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 

	DIMENSION GGG(NPTS), VSHAPE(NPTS)
	REAL WW
	GGG=0.D0
!
!	DIMENSION GG(KNBMAX),PP(KNBMAX)

C 
C    Call TRYSBA to calculate subareas. 
C 
        CALL TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX)
		GGG(NINDEX)=1.D0 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators.  There is no need to worry about the edge case, 
C    because we do not try to compensate in the C0 natural neighbour 
C    interpolant. 
C 
    1   S0 = 0.0 
        T0 = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
        LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 		
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 3
	  IF (N.NE.NINDEX) THEN
        S = SBAREA(K)*VSHAPE(N)
	  GGG(N)=S
        ZPT = VAL(N) 		
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        T0 = T0+S*ZPT
		ELSE
		S1=SBAREA(K)		
		V1=VAL(N)
		ENDIF
    3   CONTINUE 

	WW=S0+S1
	Z0=1.D0/WW*(T0+(WW-S0)*V1)
	
	DO I=1,NPTS
	IF (I.EQ.NINDEX) THEN
	GGG(NINDEX)=1.D0-(S0/WW)
	ELSE
	GGG(I)=GGG(I)/WW
	ENDIF
	ENDDO

C 
C    Calculate value for return. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GRLD 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE GRLD(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VAL, 
     1    SBAREA,PTOFF,VALOFF,KNBMAX,GRAD) 
C 
C    Loads the gradient at each data site into GRAD. 
C 
C    User provides the appropriate part of the TILE4 database, and 
C    the values at the accepted points (data sites) in VAL. 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VALOFF(KNBMAX) are needed as 
C    working space. 
C 
C    If KNBMAX is too small, error return with IFLAG set to 10 occurs. 
C    Otherwise IFLAG is zero on return.  For accepted point NPT, the 
C    gradient is loaded into (GRAD(1,NPT),GRAD(2,NPT)). 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP),VAL(NPTS), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VALOFF(KNBMAX),GRAD(2,NPTS) 
C 
C    Scan the points. 
C 
        DO 1  NPT = 1,NPTS 
C 
C    Calculate subareas. 
C 
        CALL ACCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NPT,AREA,SBAREA, 
     1    PTOFF,KNBMAX,KNB) 
C 
C    Check IFLAG.  Return if 10, loop if 9. 
C 
        IF(IFLAG.EQ.0) GOTO 2 
        IF(IFLAG.EQ.9) GOTO 1 
        RETURN 
C 
C    Calculate gradient. 
C 
    2   CALL CURLYD(NPTS,NADDR,L,LTOP,IFLAG,VAL,NPT,SBAREA,PTOFF,VALOFF, 
     1    KNB,BETAX,BETAY,GAMMA,ICASE) 
C 
C    Check IFLAG to make sure nothing has gone wrong. 
C 
        IF(IFLAG.EQ.0) GOTO 3 
        STOP 766 
C 
C    Save gradient. 
C 
    3   GRAD(1,NPT) = BETAX 
        GRAD(2,NPT) = BETAY 
    1   CONTINUE 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CURLYD 4.4 DATED 13 JANUARY 1981 
C 
        SUBROUTINE CURLYD(NPTS,NADDR,L,LTOP,IFLAG,VAL,NINDEX, 
     1    SBAREA,PTOFF,VALOFF,KNB,BETAX,BETAY,GAMMA,ICASE) 
C 
C    Applies the discrete differentiator CURLY-D at NINDEX. 
C 
C    NPTS,NADDR,L,LTOP are from the TILE4 database.  VAL contains values 
C    at accepted points.  NINDEX is the index of the point (data site) 
C    at which CURLY-D is to be applied.  KNB is the number of its 
C    neighbours, and for K = 1,...,KNB the subtile area for the Kth 
C    neighbour is held in SBAREA(K) and the offset to it as a real or 
C    virtual point in (PTOFF(1,K),PTOFF(2,K)).  KNB,SBAREA,PTOFF can be 
C    loaded by a call to ACCSBA.  VALOFF(K) returns the real or virtual 
C    value offset at the Kth neighbour, (BETAX,BETAY) the gradient and 
C    GAMMA the spherical coefficient at NINDEX.  IFLAG returns the value 
C    zero unless NINDEX is not the index of an accepted point, in which 
C    case it returns the value 9.  ICASE returns a value indicating the 
C    extent of the edge effects, according to the following conventions. 
C 
C       0  All neighbours of NINDEX are points.  Equations for BETAX, 
C          BETAY,GAMMA assumed well-conditioned.  Problems could arise 
C          only if NINDEX is too close to another point.  This should 
C          have been trapped earlier. 
C 
C       1  Three or more neighbours of NINDEX, but not all, are points 
C          and equations are well-conditioned.  This is the normal edge 
C          case.  The calculation is similar to case 0 but with extra 
C          terms to compensate for the edge effects. 
C 
C       2  Two neighbours of NINDEX are points.  More precisely, the 
C          equations for BETAX,BETAY,GAMMA are ill-conditioned or 
C          indeterminate, but on putting GAMMA equal to zero and 
C          discarding an equation, well-conditioned equations for 
C          BETAX,BETAY result. 
C 
C       3  One neighbour of NINDEX is a point.  More precisely, the 
C          equations for BETAX,BETAY,GAMMA are ill-conditioned or 
C          indeterminate, and on putting GAMMA equal to zero and 
C          discarding an equation, the resultant equations for 
C          BETAX,BETAY are still ill-conditioned or indeterminate. 
C          A direction for beta is constructed conventionally, and an 
C          appropriate magnitude is then determined. 
C 
C       4  No neighbours of NINDEX are points.  This occurs only if 
C          NINDEX is the sole accepted point.  The case is included 
C          for completeness, but will not normally arise. 
C 
C2      INTEGER*2 L 
        DIMENSION NADDR(NPTS),L(LTOP),VAL(NPTS),SBAREA(KNB), 
     1    PTOFF(2,KNB),VALOFF(KNB) 
C 
C    Test value for conditioning. 
C 
        DATA EPS1,EPS2/1.0E-6,1.0E-6/ 
C 
C    Check that NINDEX is the index of an accepted point, and if so pick 
C    up LLO as the base of its contiguity list.  If not, return with 
C    IFLAG set to 9. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    ZINDEX as the value at NINDEX.  Initialise accumulator variables. 
C 
    2   LHI = LLO-1+L(LLO-1) 
        ZINDEX = VAL(NINDEX) 
        AREA = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
        SUMSQ = 0.0 
        Q = 0.0 
        HXX = 0.0 
        HXY = 0.0 
        HYY = 0.0 
        PX = 0.0 
        PY = 0.0 
C 
C    Scan the neighbours, counting how many are points and how many are 
C    constraints.  Load VALOFF for points, and accumulate all necessary 
C    values. 
C 
        NBR = 0 
        JBR = 0 
        K = 0 
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N) 4,701,5 
  701   STOP 767 
    4   JBR = JBR+1 
        GOTO 3 
    5   NBR = NBR+1 
        Z = VAL(N)-ZINDEX 
        VALOFF(K) = Z 
        S = SBAREA(K) 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        AREA = AREA+S 
        EDGEX = EDGEX+S*U 
        EDGEY = EDGEY+S*V 
        SUMSQ = SUMSQ+S*(U*U+V*V) 
        Q = Q+S*Z 
        S = S/(U*U+V*V) 
        HXX = HXX+S*U*U 
        HXY = HXY+S*U*V 
        HYY = HYY+S*V*V 
        PX = PX+S*Z*U 
        PY = PY+S*Z*V 
    3   CONTINUE 
        IF(K.NE.KNB) STOP 770 
C 
C    Split off edge cases, and deal with case 0. 
C 
        IF(JBR.GT.0) GOTO 6 
        ICASE = 0 
        DET = HXX*HYY-HXY*HXY 
        BETAX = (HYY*PX-HXY*PY)/DET 
        BETAY = (HXX*PY-HXY*PX)/DET 
        GAMMA = Q/SUMSQ 
        IFLAG = 0 
        RETURN 
C 
C    Split off bad edge cases, and try for case 1. 
C 
    6   IF(NBR.LE.2) GOTO 7 
        ICASE = 1 
        HXXMOD = HXX-EDGEX*EDGEX/SUMSQ 
        HXYMOD = HXY-EDGEX*EDGEY/SUMSQ 
        HYYMOD = HYY-EDGEY*EDGEY/SUMSQ 
        DET = HXXMOD*HYYMOD-HXYMOD*HXYMOD 
        IF(DET/AREA.LT.EPS1) GOTO 7 
        PXMOD = PX-Q*EDGEX/SUMSQ 
        PYMOD = PY-Q*EDGEY/SUMSQ 
        BETAX = (HYYMOD*PXMOD-HXYMOD*PYMOD)/DET 
        BETAY = (HXXMOD*PYMOD-HXYMOD*PXMOD)/DET 
        GAMMA = (Q-EDGEX*BETAX-EDGEY*BETAY)/SUMSQ 
        GOTO 10 
C 
C    Try for case 2. 
C 
    7   IF(NBR.LE.1) GOTO 8 
        ICASE = 2 
        DET = HXX*HYY-HXY*HXY 
        IF(DET/AREA.LT.EPS2) GOTO 8 
        BETAX = (HYY*PX-HXY*PY)/DET 
        BETAY = (HXX*PY-HXY*PX)/DET 
        GAMMA = 0.0 
        GOTO 10 
C 
C    Case 3. 
C 
    8   IF(NBR.EQ.0) GOTO 9 
        ICASE = 3 
        BETAX = PX/AREA 
        BETAY = PY/AREA 
        GAMMA = 0.0 
        GOTO 10 
C 
C    Case 4. 
C 
    9   ICASE = 4 
        BETAX = 0.0 
        BETAY = 0.0 
        GAMMA = 0.0 
C 
C    Calculate and load virtual value offsets. 
C 
   10   LL = LLO-1 
        DO 11  K = 1,KNB 
        LL = LL+1 
        IF(L(LL).GE.0) GOTO 11 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        VALOFF(K) = BETAX*U+BETAY*V+GAMMA*(U*U+V*V) 
   11   CONTINUE 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYSBG 4.1 DATED 6 MARCH 1981 
C 
        SUBROUTINE TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    For a trial point (XPT,YPT) calculates the contiguities this point 
C    would have if inserted into the tessellation, the offsets to its 
C    neighbours, the area of its tile, and the areas and gradients 
C    of the areas of its subtiles. 
C 
C    The user supplies the TILE4 data structure, the coordinates 
C    (XPT,YPT) of the trial point, unset variables NINDEX,AREA,KNB, 
C    ICASE, and arrays SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX). 
C    NINDEX returns the index of the accepted point nearest to the trial 
C    point.  AREA returns the area of the tile of the trial point. 
C    KNB returns the number of neighbours.  ICASE returns 0 if all these 
C    are points, and 5 if some are constraints.  The data item 
C    containing the contiguity list is returned at LPTR in the heap L, 
C    and neighbours are referenced in the order in which they occur in 
C    this contiguity list.  SBAREA(K) for K = 1,...,KNB returns the area 
C    of the subtile for the Kth neighbour.  (DELSBA(1,K),DELSBA(2,K)) 
C    for K = 1,...,KNB return the components of the gradient of the 
C    area of the subtile for the Kth neighbour.  (PTOFF(1,K),PTOFF(2,K)) 
C    for K = 1,...,KNB return the offsets to the Kth neighbour as a real 
C    or virtual point.  IFLAG is zero on successful return.  Nonzero 
C    values indicate error returns as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If error return with IFLAG set to 5 occurs, 
C    NINDEX correctly returns the index of the duplicated point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Call TRYPT to find the contiguity list for the trial point. 
C 
        CALL TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Check IFLAG and return if it is nonzero. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Call LOCSBG for the contiguity list at LPTR constructed for the 
C    trial point (XPT,YPT), to work out the subareas and gradients. 
C 
    1   CALL LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LPTR,XPT,YPT, 
     1    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBG, and these values can 
C    be passed straight back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ACCSBG 4.1 DATED 6 MARCH 1981 
C 
        SUBROUTINE ACCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB) 
C 
C    Calculates the coordinate offsets of the neighbours of NINDEX, the 
C    areas of their subtiles of the tile of NINDEX, the gradients of 
C    these subtile areas, and the total area of the tile of NINDEX. 
C 
C    Note: this routine is provided as a tool for future use.  It is not 
C    called as part of the standard interpolation procedure in TILE 4. 
C 
C    The user supplies arrays SBAREA(KNBMAX),DELSBA(2,KNBMAX), 
C    PTOFF(2,KNBMAX).  KNB returns the number of neighbours of NINDEX. 
C    If this would exceed KNBMAX, error return with IFLAG set to 10 
C    occurs.  If NINDEX is not the index of an accepted point, error 
C    return with IFLAG set to 9 occurs.  Otherwise IFLAG is zero on 
C    successful return.  AREA returns the total area of the tile of 
C    NINDEX.  (PTOFF(1,K),PTOFF(2,K)) for K = 1,...,KNB return the 
C    coordinate offsets of the neighbours of NINDEX as real or virtual 
C    points, in the order in which they are encountered in its 
C    contiguity list.  SBAREA(K) for K = 1,...,KNB returns the area of 
C    the subtile associated with the Kth neighbour, (DELSBA(1,K), 
C    DELSBA(2,K)) the components of the gradient of its area. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check that NINDEX is the index of an accepted point.  If not, 
C    return with IFLAG set to 9.  If so, LOCN is the address of its 
C    data item in the heap. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN.GT.0) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    2   XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Call subroutine LOCSBG to calculate the subareas and their 
C    gradients. 
C 
        CALL LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XINDEX,YINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBG, and these values can 
C    just be passed back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCSBG 4.2 DATED 6 MARCH 1981 
C 
        SUBROUTINE LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XPT,YPT,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Calculates coordinate offsets to neighbours, tile area, subtile 
C    areas, and gradients of subtile areas for a point (XPT,YPT) whose 
C    data item in the heap is held at LOCN.  The data item may be either 
C    that of an accepted point, or that of a trial point prepared at 
C    location LPTR by a call to TRYPT.  The user should not normally 
C    call subroutine LOCSBG directly, but rather via a call to ACCSBG 
C    to find the subareas and gradients for an accepted point, or via 
C    a call to TRYSBG to find those for a trial point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    No check is made on LOCN.  If it is not in fact the address of a 
C    data item, chaos will result.  Assuming this does not happen, pick 
C    up LLO as the base of the contiguity list, and LHI as the top of 
C    it.  Check that the list is not too long. 
C 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        IF(LHI-LLO+1.LE.KNBMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Find the coordinate offsets of all the neighbours.  Initialise 
C    DELSBA to zero. 
C 
    1   KNB = 0 
        DO 2  LL = LLO,LHI 
        KNB = KNB+1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XPT,YPT,L(LL),U,V) 
        PTOFF(1,KNB) = U 
        PTOFF(2,KNB) = V 
        DELSBA(1,KNB) = 0.0 
    2   DELSBA(2,KNB) = 0.0 
C 
C    Initialise AREA to zero, and UCURR,VCURR,DCURR to the components 
C    and squared length of the last offset.  Initial and ICASE. 
C 
        AREA = 0.0 
        UCURR = PTOFF(1,KNB) 
        VCURR = PTOFF(2,KNB) 
        DCURR = UCURR**2+VCURR**2 
        LL = LLO-1 
        ICASE = 0 
C 
C    The main loop scans through the neighbours.  It is a DO loop on 
C    KCURR, with LL scanning the original contiguity list in parallel. 
C 
        DO 3  KCURR = 1,KNB 
        LL = LL+1 
C 
C    Save old values and pick up new ones.  Set ICASE to 5 if KCURRth 
C    neighbour is virtual. 
C 
        UPREV = UCURR 
        VPREV = VCURR 
        DPREV = DCURR 
        UCURR = PTOFF(1,KCURR) 
        VCURR = PTOFF(2,KCURR) 
        DCURR = UCURR**2+VCURR**2 
        N = L(LL) 
        IF(N.LE.0) ICASE = 5 
C 
C    The subarea for KCURR is calculated by breaking the subtile into 
C    triangles and accumulating their areas in SBA, which is initialised 
C    to zero.  All these triangles have as a common vertex the vertex of 
C    the tile clockwise from KCURR.  Find the offset of this hinge point 
C    as (UHINGE,VHINGE). 
C 
        SBA = 0.0 
        C = 0.5/(UPREV*VCURR-UCURR*VPREV) 
        UHINGE = (DPREV*VCURR-DCURR*VPREV)*C 
        VHINGE = (DCURR*UPREV-DPREV*UCURR)*C 
C 
C    The first such triangle has the KCURRth face of the tile itself as 
C    a side.  One end of this side is the hinge point.  Find the other 
C    end, as (UVTX,VVTX). 
C 
        KEXTR = KCURR+1 
        IF(KEXTR.GT.KNB) KEXTR = 1 
        UEXTR = PTOFF(1,KEXTR) 
        VEXTR = PTOFF(2,KEXTR) 
        DEXTR = UEXTR**2+VEXTR**2 
        C = 0.5/(UCURR*VEXTR-UEXTR*VCURR) 
        UVTX = (DCURR*VEXTR-DEXTR*VCURR)*C 
        VVTX = (DEXTR*UCURR-DCURR*UEXTR)*C 
C 
C    If the KCURRth neighbour is real, calculate the main term in the 
C    subtile gradient - the only one if all neighbours are real. 
C 
        IF(N.LE.0) GOTO 4 
        FACTOR = 0.5*SQRT(((UVTX-UHINGE)**2+(VVTX-VHINGE)**2)/DCURR) 
        DELSBA(1,KCURR) = DELSBA(1,KCURR)+FACTOR*(UVTX+UHINGE) 
        DELSBA(2,KCURR) = DELSBA(2,KCURR)+FACTOR*(VVTX+VHINGE) 
C 
C    Now enter a test-controlled inner loop to find the remaining 
C    vertices of the subtile in anticlockwise order.  These vertices are 
C    internal to the tile.  Each new vertex, its predecessor, and the 
C    hinge point together form a triangle, and these triangles give the 
C    desired decomposition of the subtile.  We begin by saving old 
C    values. 
C 
    4   KSCAN = KEXTR 
        USCAN = UEXTR 
        VSCAN = VEXTR 
        UOLD = UVTX 
        VOLD = VVTX 
C 
C    Attempt to find the next vertex.  The candidates for defining the 
C    next vertex are all the neighbours following KSCAN and 
C    preceding KCURR.  They are scanned in an inner test-controlled loop 
C    to find which one gives the tightest intercept on the perpendicular 
C    bisector of (UCURR,VCURR) and (USCAN,VSCAN).  The winner, or the 
C    last such in the case of ties, becomes KEXTR.  If there are no 
C    such candidates, KEXTR is not advanced from KSCAN, the attempt to 
C    find a new vertex fails, and the subtile has been exhausted.  In 
C    that case set (UVTX,VVTX) to (UHINGE,VHINGE) for the benefit of 
C    the gradient calculation if KCURR is virtual. 
C 
        K = KSCAN 
    5   K = K+1 
        IF(K.GT.KNB) K = 1 
        IF(K.EQ.KCURR) GOTO 6 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        D = (U-USCAN)*(VCURR-VSCAN)-(UCURR-USCAN)*(V-VSCAN) 
        IF(D.LE.0.0) GOTO 5 
        T = ((U-USCAN)*(U-UCURR)+(V-VSCAN)*(V-VCURR))/D 
        IF(KEXTR.EQ.KSCAN) GOTO 7 
        IF(TEXTR-T) 5,8,7 
    7   TEXTR = T 
    8   KEXTR = K 
        UEXTR = U 
        VEXTR = V 
        GOTO 5 
    6   IF(KEXTR.NE.KSCAN) GOTO 9 
        UVTX = UHINGE 
        VVTX = VHINGE 
        GOTO 10 
C 
C    The new vertex has been identified.  Find its offset as (UVTX,VVTX) 
C    and hence find the area of the triangle and add it to SBA. 
C 
    9   U2 = USCAN-UCURR 
        V2 = VSCAN-VCURR 
        U3 = UEXTR-UCURR 
        V3 = VEXTR-VCURR 
        S2 = U2*(USCAN+UCURR)+V2*(VSCAN+VCURR) 
        S3 = U3*(UEXTR+UCURR)+V3*(VEXTR+VCURR) 
        C = 0.5/(U2*V3-U3*V2) 
        UVTX = (S2*V3-S3*V2)*C 
        VVTX = (S3*U2-S2*U3)*C 
        SBA = SBA+0.5*((UOLD-UHINGE)*(VVTX-VHINGE)-(UVTX-UHINGE)* 
     1    (VOLD-VHINGE)) 
C 
C    If the KCURRth neighbour is virtual then an infinitesimal 
C    displacement of the point itself causes an infinitesimal 
C    displacement of the KCURRth neighbour, related to that of the 
C    point by reflexion in the boundary.  The boundary is left fixed, 
C    and rather than calculate the usual term in the gradient for 
C    the corresponding subtile and then cancel it, we have simply 
C    omitted it.  But two systems of compensating terms are needed, 
C    which cancel one another in pairs but are allocated to different 
C    sums and so have an effect.  Those allocated to the gradient of the 
C    subtile area of the KCURRth neighbour are, like that subtile area 
C    itself, omitted from the summations involved in interpolant 
C    calculations and are thus in practice lost, but the balancing terms 
C    paired with them are allocated to the gradients of the areas of 
C    those subtiles contiguous to that of KCURR.  Where such subtiles 
C    correspond to real points, the terms have an effect in the 
C    interpolant calculations. 
C 
C    Check if the KCURRth neighbour is virtual. 
C 
   10   IF(N.GT.0) GOTO 11 
C 
C    It is.  Find the gradient contribution arising from the 
C    KCURR - KSCAN common boundary w.r.t. displacement of KCURR. 
C 
        FACTOR = 0.5*SQRT(((UVTX-UOLD)**2+(VVTX-VOLD)**2)/ 
     1    ((USCAN-UCURR)**2+(VSCAN-VCURR)**2)) 
        GX = FACTOR*(UVTX+UOLD-2.0*UCURR) 
        GY = FACTOR*(VVTX+VOLD-2.0*VCURR) 
C 
C    Reflect it in the generating point - KCURR boundary to obtain the 
C    gradient contribution w.r.t. displacement of the generating 
C    point. 
C 
        FACTOR = 2.0*(GX*UCURR+GY*VCURR)/DCURR 
        GX = GX-FACTOR*UCURR 
        GY = GY-FACTOR*VCURR 
C 
C    Add it into the gradient for KCURR to keep the books straight 
C    even though we do not use this value, and subtract it from the 
C    gradient for KSCAN, which is an effective change unless KSCAN 
C    is itself virtual. 
C 
        DELSBA(1,KCURR) = DELSBA(1,KCURR)+GX 
        DELSBA(2,KCURR) = DELSBA(2,KCURR)+GY 
        DELSBA(1,KSCAN) = DELSBA(1,KSCAN)-GX 
        DELSBA(2,KSCAN) = DELSBA(2,KSCAN)-GY 
   11   IF(KEXTR.NE.KSCAN) GOTO 4 
C 
C    Store the subarea and add it to the area. 
C 
        SBAREA(KCURR) = SBA 
    3   AREA = AREA+SBA 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYSBA 4.4 DATED 6 MARCH 1981 
C 
        SUBROUTINE TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    For a trial point (XPT,YPT) calculates the contiguities this point 
C    would have if inserted into the tessellation, the offsets to its 
C    neighbours, the area of its tile, and the areas of its subtiles. 
C 
C    The user supplies the TILE4 data structure, the coordinates 
C    (XPT,YPT) of the trial point, unset variables NINDEX,AREA,KNB, 
C    ICASE, and arrays SBAREA(KNBMAX),PTOFF(2,KNBMAX).  NINDEX returns 
C    the index of the accepted point nearest to the trial point.  AREA 
C    returns the area of the tile of the trial point.  KNB returns the 
C    number of neighbours.  ICASE returns 0 if all these are points, and 
C    5 if some are constraints.  The data item containing the contiguity 
C    list is returned at LPTR in the heap L, and neighbours are 
C    referenced in the order in which they occur in this contiguity 
C    list.  SBAREA(K) for K = 1,...,KNB returns the area of the subtile 
C    for the Kth neighbour.  (PTOFF(1,K),PTOFF(2,K)) for K = 1,...,KNB 
C    return the offsets to the Kth neighbour as a real or virtual point. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap  overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If error return with IFLAG set to 5 occurs, 
C    NINDEX correctly returns the index of the duplicated point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Call TRYPT to find the contiguity list for the trial point. 
C 
        CALL TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Check IFLAG and return if it is nonzero. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Call LOCSBA for the contiguity list at LPTR constructed for the 
C    trial point (XPT,YPT), to work out subareas. 
C 
    1   CALL LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LPTR,XPT,YPT, 
     1    AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBA, and these values can 
C    be passed straight back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ACCSBA 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE ACCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    AREA,SBAREA,PTOFF,KNBMAX,KNB) 
C 
C    Calculates the coordinate offsets of the neighbours of NINDEX, 
C    the areas of their subtiles of the tile of NINDEX, and the area 
C    of the tile of NINDEX. 
C 
C    The user supplies arrays SBAREA(KNBMAX),PTOFF(2,KNBMAX).  KNB 
C    returns the number of neighbours of NINDEX.  If this would exceed 
C    KNBMAX, error return with IFLAG set to 10 occurs.  If NINDEX is 
C    not the index of an accepted point, error return with IFLAG set to 
C    9 occurs.  Otherwise IFLAG is zero on successful return.  AREA 
C    returns the total area of the tile of NINDEX.  (PTOFF(1,K), 
C    PTOFF(2,K)) for K = 1,...,KNB return the coordinate offsets of 
C    the neighbours of NINDEX as actual or virtual points, in the order 
C    in which they are encountered in its contiguity list.  SBAREA(K) 
C    for K = 1,...,KNB returns the area of the subtile associated with 
C    the Kth neighbour. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check that NINDEX is the index of an accepted point.  If not, 
C    return with IFLAG set to 9.  If so, LOCN is the address of its 
C    data item in the heap. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN.GT.0) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    2   XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Call subroutine LOCSBA to calculate the subareas. 
C 
        CALL LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XINDEX,YINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBA, and these values 
C    can just be passed back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCSBA 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XPT,YPT,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Calculates coordinate offsets to neighbours, tile area, and 
C    subtile areas for a point (XPT,YPT) whose data item in the heap 
C    is held at LOCN.  The data item may be either that of an 
C    accepted point, or that of a trial point prepared at location 
C    LPTR by a call to TRYPT.  The user should not normally call 
C    subroutine LOCSBA directly, but rather via a call to ACCSBA to 
C    find the subareas for an accepted point, or via a call to TRYSBA 
C    to find those for a trial point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    No check is made on LOCN.  If it is not in fact the address of a 
C    data item, chaos will result.  Assuming this does not happen, 
C    pick up LLO as the base of the contiguity list, and LHI as 
C    the top of it.  Check that the list is not too long. 
C 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        IF(LHI-LLO+1.LE.KNBMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Find the coordinate offsets of all the neighbours. 
C 
    1   KNB = 0 
        DO 4  LL = LLO,LHI 
        KNB = KNB+1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XPT,YPT,L(LL),U,V) 
        PTOFF(1,KNB) = U 
    4   PTOFF(2,KNB) = V 
C 
C    Initialise AREA to zero, and UCURR,VCURR,DCURR to the components 
C    and squared length of the last offset. 
C 
        AREA = 0.0 
        UCURR = PTOFF(1,KNB) 
        VCURR = PTOFF(2,KNB) 
        DCURR = UCURR*UCURR+VCURR*VCURR 
        LL = LLO-1 
        ICASE = 0 
C 
C    The main loop scans through the neighbours.  Save old values and 
C    pick up new ones. 
C 
        DO 5  KCURR = 1,KNB 
        LL = LL+1 
        UPREV = UCURR 
        VPREV = VCURR 
        DPREV = DCURR 
        UCURR = PTOFF(1,KCURR) 
        VCURR = PTOFF(2,KCURR) 
        DCURR = UCURR*UCURR+VCURR*VCURR 
        IF(L(LL).LE.0) ICASE = 5 
C 
C    The subarea for KCURR is calculated by breaking the subtile into 
C    triangles and accumulating their areas in SBA, which is initialised 
C    to zero.  All these triangles have as a common vertex the vertex of 
C    the tile clockwise from KCURR.  Find the offset of this hinge point 
C    as (UHINGE,VHINGE). 
C 
        SBA = 0.0 
        C = 0.5/(UPREV*VCURR-UCURR*VPREV) 
        UHINGE = (DPREV*VCURR-DCURR*VPREV)*C 
        VHINGE = (DCURR*UPREV-DPREV*UCURR)*C 
C 
C    The first such triangle has the Kth face of the tile itself as a 
C    side.  One end of this side is the hinge point.  Find the other 
C    end, as (UVTX,VVTX). 
C 
        KEXTR = KCURR+1 
        IF(KEXTR.GT.KNB) KEXTR = 1 
        UEXTR = PTOFF(1,KEXTR) 
        VEXTR = PTOFF(2,KEXTR) 
        DEXTR = UEXTR*UEXTR+VEXTR*VEXTR 
        C = 0.5/(UCURR*VEXTR-UEXTR*VCURR) 
        UVTX = (DCURR*VEXTR-DEXTR*VCURR)*C 
        VVTX = (DEXTR*UCURR-DCURR*UEXTR)*C 
C 
C    Now enter a test-controlled inner loop to find the remaining 
C    vertices of the subtile in anticlockwise order.  These vertices 
C    are internal to the tile.  Each new vertex, its predecessor, and 
C    the hinge point together form a triangle, and these triangles give 
C    the desired decomposition of the subtile.  We begin by saving 
C    old values. 
C 
    6   KSCAN = KEXTR 
        USCAN = UEXTR 
        VSCAN = VEXTR 
        UOLD = UVTX 
        VOLD = VVTX 
C 
C    Attempt to find the next vertex.  The candidates for defining the 
C    next vertex are the neighbours following KSCAN and 
C    preceding KCURR.  They are scanned in an inner test-controlled loop 
C    to find which one gives the tightest intercept on the perpendicular 
C    bisector of (UCURR,VCURR) and (USCAN,VSCAN).  The winner, or the 
C    last such in the case of ties, becomes KEXTR.  If there are no 
C    such candidates, KEXTR is not advanced from KSCAN, the attempt to 
C    find a new vertex fails, and the subtile has been exhausted. 
C 
        K = KSCAN 
    7   K = K+1 
        IF(K.GT.KNB) K = 1 
        IF(K.EQ.KCURR) GOTO 8 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        D = (U-USCAN)*(VCURR-VSCAN)-(UCURR-USCAN)*(V-VSCAN) 
        IF(D.LE.0.0) GOTO 7 
        T = ((U-USCAN)*(U-UCURR)+(V-VSCAN)*(V-VCURR))/D 
        IF(KEXTR.EQ.KSCAN) GOTO 9 
        IF(TEXTR-T) 7,10,9 
    9   TEXTR = T 
   10   KEXTR = K 
        UEXTR = U 
        VEXTR = V 
        GOTO 7 
    8   IF(KEXTR.EQ.KSCAN) GOTO 11 
C 
C    The new vertex has been identified.  Find its offset as (UVTX,VVTX) 
C    and hence find the area of the triangle and add it to SBA. 
C 
        U2 = USCAN-UCURR 
        V2 = VSCAN-VCURR 
        U3 = UEXTR-UCURR 
        V3 = VEXTR-VCURR 
        S2 = U2*(USCAN+UCURR)+V2*(VSCAN+VCURR) 
        S3 = U3*(UEXTR+UCURR)+V3*(VEXTR+VCURR) 
        C = 0.5/(U2*V3-U3*V2) 
        UVTX = (S2*V3-S3*V2)*C 
        VVTX = (S3*U2-S2*U3)*C 
        SBA = SBA+0.5*((UOLD-UHINGE)*(VVTX-VHINGE)-(UVTX-UHINGE)* 
     1    (VOLD-VHINGE)) 
        GOTO 6 
C 
C    Store the subarea and add it to the area. 
C 
   11   SBAREA(KCURR) = SBA 
    5   AREA = AREA+SBA 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    GRAPHICS ROUTINES LAST MODIFIED 13 JANUARY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE WINPLT 4.1 DATED 25 JULY 1979 
C 
        SUBROUTINE WINPLT(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,IWIN) 
C 
C    Initialises the plot, and plots the window in logical pen colour 
C    IWIN.  A negative value of IWIN suppresses plotting, but 
C    still initialises the plot. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP),VTX(2,KVTMAX) 
        COMMON/PLTBLK/DUMMY 
C 
C    Call WINVTX to calculate the window vertices and rectangular 
C    frame. 
C 
        CALL WINVTX(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,KVTX, 
     1    XMIN,XMAX,YMIN,YMAX) 
C 
C    Check IFLAG and return if it is nonzero. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Call PLTON to initialise the plot. 
C 
    1   CALL PLTON(XMIN,XMAX,YMIN,YMAX) 
C 
C    Return if IWIN is negative, otherwise set the pen colour 
C    for plotting the window. 
C 
        IF(IWIN.LT.0) RETURN 
        CALL PLTPEN(IWIN) 
C 
C    Plot the window. 
C 
        CALL PLTMOV(VTX(1,KVTX),VTX(2,KVTX)) 
        DO 2  K = 1,KVTX 
    2   CALL PLTLIN(VTX(1,K),VTX(2,K)) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE POLPLT 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE POLPLT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG, 
     1    VTX,KVTMAX,NCODE,NTAG,IDIFF,ISAME,IPTS,IMARK) 
C 
C    Plots the tessellation and the points.  NCODE is an array of 
C    point codes which must be positive for accepted points.  Entries 
C    in this array are modified and restored by this routine.  The 
C    boundary between tiles whose points have distinct codes is plotted 
C    in logical pen colour IDIFF, and logical pen colour ISAME is used 
C    where the codes are the same.  NTAG is an array of identifiers for 
C    the points.  Where appropriate NTAG and NCODE may be the same 
C    array, but neither should coincide with NADDR.  Logical pen colour 
C    IPTS is used for plotting the points.  In all cases a negative 
C    logical pen colour suppresses plotting efficiently.  IMARK controls 
C    the style of point plotting, in a manner dependent on how the 
C    interface routine PLTPT is implemented.  The usual convention 
C    is as follows. 
C 
C       1  no position mark, integer annotation 
C       2  no position mark, character annotation 
C       3  position mark, no annotation 
C       4  position mark, integer annotation 
C       5  position mark, character annotation 
C 
C    If IDIFF and ISAME are both zero, then points alone are plotted in 
C    an efficient manner. 
C 
C    REMINDER: Entries in NCODE must be POSITIVE.  Zero is NOT positive. 
C 
C    KVTMAX must be at least as large as the number of vertices of any 
C    one tile.  Error return with IFLAG set to 10 takes place if this 
C    does not hold.  In this case the plot is abandoned tidily, and 
C    NCODE is restored.  IFLAG is zero on normal return. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX),NCODE(NPTS),NTAG(NPTS) 
        COMMON/PLTBLK/DUMMY 
C 
C    Anything to do?  If not, return immediately. 
C 
        IF(IDIFF.GE.0.OR.ISAME.GE.0.OR.IPTS.GE.0) GOTO 1 
        IFLAG = 0 
        RETURN 
C 
C    Initialise the pen colour as unset. 
C 
    1   IPEN = -1 
C 
C    Scan the points to make sure nothing is missed, but drop without 
C    further ado indices which do not reference accepted points, also 
C    points already dealt with. 
C 
        DO 2  NPT = 1,NPTS 
        IF(NADDR(NPT).LE.0) GOTO 2 
        IF(NCODE(NPT).LE.0) GOTO 2 
C 
C    Start a chain of points. 
C 
        NNEW = NPT 
C 
C    Pursue the chain. 
C 
    3   NCURR = NNEW 
        MC = NCODE(NCURR) 
C 
C    Points only?  If so skip the vertex calculations. 
C 
        IF(IDIFF.LT.0.AND.ISAME.LT.0) GOTO 4 
C 
C    Some plotting likely on the tile boundary, so find it. 
C 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NCURR, 
     1    VTX,KVTMAX,KVTX,.FALSE.) 
C 
C    Check IFLAG.  Tidy up and go home if it is 10.  It ought not 
C    to be 9, as we know already that NCURR is an accepted point. 
C 
        IF(IFLAG.EQ.0) GOTO 4 
        IF(IFLAG.EQ.9) STOP 756 
        CALL PLTOFF 
        DO 5  N = 1,NPTS 
        MC = NCODE(N) 
        IF(NADDR(N).GT.0.AND.MC.LT.0) NCODE(N) = -MC 
    5   CONTINUE 
        RETURN 
C 
C    Go round the tile boundary plotting in the appropriate way anything 
C    that needs plotting. 
C 
    4   LLO = NADDR(NCURR)+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        KPEN = 0 
        DO 6  LL = LLO,LHI 
        K = K+1 
        K1 = K+1 
        IF(K1.GT.KVTX) K1 = 1 
        N = L(LL) 
        IF(N.LE.0) GOTO 6 
        M = NCODE(N) 
        IF(M.LE.0) GOTO 6 
        NNEW = N 
        IF(M.EQ.MC) GOTO 7 
        IF(IDIFF.LT.0) GOTO 6 
        IF(IPEN.EQ.IDIFF) GOTO 8 
        IPEN = IDIFF 
        CALL PLTPEN(IPEN) 
        GOTO 8 
    7   IF(ISAME.LT.0) GOTO 6 
        IF(IPEN.EQ.ISAME) GOTO 8 
        IPEN = ISAME 
        CALL PLTPEN(IPEN) 
    8   IF(KPEN.LT.K) CALL PLTMOV(VTX(1,K),VTX(2,K)) 
        CALL PLTLIN(VTX(1,K1),VTX(2,K1)) 
        KPEN = K1 
    6   CONTINUE 
C 
C    Plot the point if wanted. 
C 
        IF(IPTS.LT.0) GOTO 9 
        IF(IPEN.EQ.IPTS) GOTO 10 
        IPEN = IPTS 
        CALL PLTPEN(IPEN) 
   10   CALL PLTPT(PT(1,NCURR),PT(2,NCURR),NTAG(NCURR),IMARK) 
C 
C    Flag the code with a minus sign to show that the point is done. 
C 
    9   NCODE(NCURR) = -MC 
C 
C    Pursue the chain if it is still alive, otherwise find the next 
C    unplotted point. 
C 
        IF(NNEW.NE.NCURR) GOTO 3 
    2   CONTINUE 
C 
C    Restore the code values. 
C 
        DO 11  N = 1,NPTS 
        MC = NCODE(N) 
        IF(NADDR(N).GT.0.AND.MC.LT.0) NCODE(N) = -MC 
   11   CONTINUE 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRIPLT 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TRIPLT(JCNS,PT,NPTS,NADDR,L,LTOP,NCODE,NTAG,IDIFF, 
     1    ISAME,IPTS,IMARK) 
C 
C    Plots the Delaunay triangulation and the points.  NCODE is an array 
C    of point codes which must be positive for accepted points.  Entries 
C    in this array, and entries in the heap L, are modified and restored 
C    by this routine.  Edges in the triangulation linking points with 
C    distinct codes are plotted in logical pen colour IDIFF, and logical 
C    pen colour ISAME is used where the codes are the same.  NTAG is an 
C    array of identifiers for the points.  Where appropriate NTAG and 
C    NCODE may be the same array, but neither should coincide with 
C    NADDR.  Logical pen colour IPTS is used for plotting the points. 
C    In all cases a negative logical pen colour suppresses plotting. 
C    IMARK controls the style of point plotting, in a manner dependent 
C    on how the interface routine PLTPT has been implemented.  The usual 
C    convention is as follows. 
C 
C       1  no position mark, integer annotation 
C       2  no position mark, character annotation 
C       3  position mark, no annotation 
C       4  position mark, integer annotation 
C       5  position mark, character annotation 
C 
C    If IDIFF and ISAME are both negative, points alone will be plotted, 
C    but not as efficiently as if POLPLT is used. 
C 
C    REMINDER: Entries in NCODE must be POSITIVE.  Zero is NOT positive. 
C 
C    No transfer array is used by this routine, and there are no error 
C    returns. 
C 
C2      INTEGER*2 L 
        DIMENSION PT(2,NPTS),NADDR(NPTS),L(LTOP),NCODE(NPTS),NTAG(NPTS) 
        COMMON/PLTBLK/DUMMY 
C 
C    Anything to do?  If not, return immediately. 
C 
        IF(IDIFF.LT.0.AND.ISAME.LT.0.AND.IPTS.LT.0) RETURN 
C 
C    Initialise the pen colour as unset. 
C 
        IPEN = -1 
C 
C    Scan the points to make sure nothing is missed, but drop without 
C    further ado indices which do not reference accepted points, also 
C    points already dealt with. 
C 
        DO 1  NPT = 1,NPTS 
        IF(NADDR(NPT).LE.0) GOTO 1 
    2   IF(NCODE(NPT).LE.0) GOTO 1 
C 
C    Start or restart a chain of points.  Move to the first point of the 
C    chain. 
C 
        NNEW = NPT 
        CALL PLTMOV(PT(1,NNEW),PT(2,NNEW)) 
C 
C    Pursue the chain. 
C 
    3   NCURR = NNEW 
        MC = NCODE(NCURR) 
C 
C    Have we visited this point before?  If not, plot it if wanted and 
C    flag it as visited. 
C 
        LLO = NADDR(NCURR)+2 
        IF(L(LLO-2).EQ.0) GOTO 4 
        L(LLO-2) = 0 
        IF(IPTS.LT.0) GOTO 4 
        IF(IPEN.EQ.IPTS) GOTO 5 
        IPEN = IPTS 
        CALL PLTPEN(IPEN) 
    5   CALL PLTPT(PT(1,NCURR),PT(2,NCURR),NTAG(NCURR),IMARK) 
C 
C    Now try to find an edge leading away from NCURR which needs to be 
C    plotted.  Failing that, a neighbour of NCURR which is not already 
C    completely dealt with will enable us to keep the chain going. 
C 
    4   LHI = LLO-1+L(LLO-1) 
        DO 6  LL = LLO,LHI 
        N = L(LL) 
C 
C    If N is positive, we have a candidate for plotting. 
C 
        IF(N.GT.0) GOTO 7 
C 
C    If N is negative but at least -JCNS, we have a constraint, which 
C    gets us nowhere, so try again. 
C 
        IF(N.GE.-JCNS) GOTO 6 
C 
C    We are looking at an edge which has been dealt with.  Unscramble it 
C    to find out where it went.  It might provide a way of extending the 
C    chain if the point it leads to is not completely dealt with. 
C 
        N = -(N+JCNS) 
        IF(NCODE(N).GT.0) NNEW = N 
        GOTO 6 
C 
C    Now look at a candidate for plotting.  Even if we do not want to 
C    plot it, it provides a way of extending the chain.  Either way, it 
C    must be flagged as dealt with from both ends. 
C 
    7   NNEW = N 
        L(LL) = -(N+JCNS) 
        LLONEW = NADDR(NNEW)+2 
        LHINEW = LLONEW-1+L(LLONEW-1) 
        DO 8  LLNEW = LLONEW,LHINEW 
        IF(L(LLNEW).EQ.NCURR) GOTO 9 
    8   CONTINUE 
        STOP 757 
    9   L(LLNEW) = -(NCURR+JCNS) 
C 
C    Does the edge need plotting?  Keep looking if not.  Otherwise 
C    plot it. 
C 
        IF(NCODE(NNEW).EQ.MC) GOTO 10 
        IF(IDIFF.LT.0) GOTO 6 
        IF(IPEN.EQ.IDIFF) GOTO 11 
        IPEN = IDIFF 
        CALL PLTPEN(IPEN) 
        GOTO 11 
   10   IF(ISAME.LT.0) GOTO 6 
        IF(IPEN.EQ.ISAME) GOTO 11 
        IPEN = ISAME 
        CALL PLTPEN(IPEN) 
   11   CALL PLTLIN(PT(1,NNEW),PT(2,NNEW)) 
C 
C    Flag NCURR if it is now completely dealt with, then move on. 
C 
        IF(LL.EQ.LHI) NCODE(NCURR) = -MC 
        GOTO 3 
    6   CONTINUE 
C 
C    If we have dropped through this DO loop, we have failed to find 
C    a plottable edge leading from NCURR, so it is completely dealt 
C    with. 
C 
        NCODE(NCURR) = -MC 
C 
C    But we might still be able to extend the chain by a pen-up move. 
C    If not, try to restart the chain from NPT. 
C 
        IF(NNEW.EQ.NCURR) GOTO 2 
        CALL PLTMOV(PT(1,NNEW),PT(2,NNEW)) 
        GOTO 3 
    1   CONTINUE 
C 
C    We have finished.  Restore NCODE and L, and return. 
C 
        DO 12  NPT = 1,NPTS 
        LLO = NADDR(NPT)+2 
        IF(LLO.LE.2) GOTO 12 
        IF(L(LLO-2).NE.0) STOP 760 
        L(LLO-2) = NPT 
        MC = NCODE(NPT) 
        IF(MC.GT.0) STOP 761 
        NCODE(NPT) = -MC 
        LHI = LLO-1+L(LLO-1) 
        DO 13  LL = LLO,LHI 
        N = L(LL) 
        IF(N.GE.0) STOP 762 
        IF(N.GE.-JCNS) GOTO 13 
        L(LL) = -(N+JCNS) 
   13   CONTINUE 
   12   CONTINUE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE DIRPLT 4.2 DATED 7 DECEMBER 1979 
C 
        SUBROUTINE DIRPLT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG, 
     1    VTX,KVTMAX,NCODE,NTAG,IWIN,IPDIFF,IPSAME,ITDIFF,ITSAME, 
     2    IPTS,IMARK) 
C 
C    A master plotting routine to plot the window, tessellation, and 
C    triangulation, and close the plot.  VTX is a transfer array of 
C    dimension (2,KVTMAX).  The second dimension should usually be 
C    set to 50.  If this is inadequate, error return with IFLAG set 
C    to 10 occurs, otherwise IFLAG returns zero.  NCODE is an array 
C    of point codes which must be POSITIVE for accepted points. 
C    Recall that zero is NOT positive.  For a pair of contiguous 
C    points with distinct codes, logical pen colour IPDIFF 
C    is used for the inter-tile boundary and logical pen colour 
C    ITDIFF for the edge in the triangulation.  Logical pen 
C    colours IPSAME,ITSAME are used if the codes are the same. 
C    The window is plotted in logical pen colour IWIN. 
C    Logical pen colour IPTS is used to mark the points.  A negative 
C    logical pen colour suppresses plotting efficiently in each case. 
C    NTAG is an array of tags for the points.  NTAG and NCODE may 
C    be the same array where appropriate, but neither should coincide 
C    with NADDR.  The style of point annotation is controlled by 
C    IMARK in a manner dependent on how the interface routine 
C    PLTPT has been implemented.  The usual convention is as follows. 
C 
C       1  no position mark, integer annotation 
C       2  no position mark, character annotation 
C       3  position mark, no annotation 
C       4  position mark, integer annotation 
C       5  position mark, character annotation 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX),NCODE(NPTS),NTAG(NPTS) 
        CALL WINPLT(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,IWIN) 
        IF(IFLAG.EQ.0) GOTO 1 
        GOTO 4 
    1   IPPTS = IPTS 
        ITPTS = -1 
        IF(ITDIFF.LT.0.AND.ITSAME.LT.0) GOTO 2 
        IPPTS = -1 
        ITPTS = IPTS 
    2   CALL POLPLT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTX,KVTMAX, 
     1    NCODE,NTAG,IPDIFF,IPSAME,IPPTS,IMARK) 
        IF(IFLAG.EQ.0) GOTO 3 
        GOTO 4 
    3   CALL TRIPLT(JCNS,PT,NPTS,NADDR,L,LTOP,NCODE,NTAG,ITDIFF, 
     1    ITSAME,ITPTS,IMARK) 
    4   CALL PLTOFF 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE DIRHAT 4.1 DATED 7 DECEMBER 1979 
C 
        SUBROUTINE DIRHAT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTX, 
     1    KVTMAX,NCODE,NTAG,IWIN,IDIFF,ISAME,IPTS,IMARK,AV,UV, 
     2    IV,JV,MCMAX) 
C 
C    A master plotting routine to plot the window and the tessellation 
C    and crosshatch the tessellation, and then close the plot.  VTX is 
C    a transfer array of dimension (2,KVTMAX).  The second dimension 
C    should usually be set to 50.  If this is inadequate, error return 
C    with IFLAG set to 10 occurs, otherwise IFLAG is zero on return. 
C    NCODE is an array of point codes, which must lie in the range 1 to 
C    MCMAX.  Logical pen colour IWIN is used for the window, IDIFF for 
C    an inter-tile boundary where the codes differ, ISAME where they 
C    are the same.  A tile with code N is crosshatched with hatch angle 
C    AV(N), scale UV(N), style IV(N), logical pen colour JV(N).  Any 
C    action may be suppressed efficiently by using a negative logical 
C    pen colour.  Crosshatching may also be suppressed by using a zero 
C    hatch style.  Points may be plotted, in logical pen colour IPTS 
C    and with annotation taken from NTAG and controlled by IMARK under 
C    the usual conventions.  NTAG may coincide with NCODE, but neither 
C    should coincide with NADDR. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX),NCODE(NPTS),NTAG(NPTS),AV(MCMAX),UV(MCMAX), 
     2    IV(MCMAX),JV(MCMAX) 
        CALL WINPLT(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,IWIN) 
        IF(IFLAG.EQ.0) GOTO 1 
        GOTO 2 
    1   CALL POLHAT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTX, 
     1    KVTMAX,NCODE,NTAG,IDIFF,ISAME,IPTS,IMARK,AV,UV,IV,JV,MCMAX) 
    2   CALL PLTOFF 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE POLHAT 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE POLHAT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTX, 
     1    KVTMAX,NCODE,NTAG,IDIFF,ISAME,IPTS,IMARK,AV,UV,IV,JV,MCMAX) 
C 
C    Plots the tessellation and the points, and crosshatches the tiles. 
C    NCODE is an array of point codes which must be positive for 
C    accepted points.  Entries in this array are modified and restored 
C    by this routine.  The boundary between tiles whose points have 
C    distinct codes is plotted in logical pen colour IDIFF, and logical 
C    pen colour ISAME is used where the codes are the same.  NTAG is 
C    an array of identifiers for the points.  Where appropriate NTAG and 
C    NCODE may be the same array, but neither should coincide with 
C    NADDR.  Logical pen colour IPTS is used for plotting the points. 
C    In all cases a negative logical pen colour suppresses plotting 
C    efficiently.  IMARK controls the style of point plotting, in a 
C    manner dependent on how the interface routine PLTPT is implemented. 
C    The usual convention is as follows. 
C 
C       1  no position mark, integer annotation 
C       2  no position mark, character annotation 
C       3  position mark, no annotation 
C       4  position mark, integer annotation 
C       5  position mark, character annotation 
C 
C    The values in NCODE control the type of crosshatching applied to 
C    each tile, by lookup in the arrays AV,UV,IV,JV.  Since these are 
C    of dimension MCMAX, the values in NCODE must lie in the range 
C    1 to MCMAX.  The Nth hatch style is specified by angle AV(N), 
C    scale UV(N), style IV(N), and logical pen colour JV(N). 
C    Hatching may be suppressed either by giving a zero hatch style or 
C    a negative logical pen colour. 
C 
C    The array VTX, of dimension (2,KVTMAX), is used as working space 
C    and must be large enough to hold the vertex description of each 
C    tile.  Usually KVTMAX is set to 50.  Error return with IFLAG set 
C    to 10 takes place if KVTMAX is too small, in which case the plot 
C    is abandoned tidily, and NCODE is restored before return.  IFLAG 
C    is zero on normal return. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX),NCODE(NPTS),NTAG(NPTS),AV(MCMAX),UV(MCMAX), 
     2    IV(MCMAX),JV(MCMAX) 
        COMMON/PLTBLK/DUMMY 
C 
C    Initialise the pen colour as unset. 
C 
        IPEN = -1 
C 
C    Scan the points to make sure nothing is missed, but drop without 
C    further ado indices which do not reference accepted points, also 
C    points already dealt with. 
C 
        DO 1  NPT = 1,NPTS 
        IF(NADDR(NPT).LE.0) GOTO 1 
        IF(NCODE(NPT).LE.0) GOTO 1 
C 
C    Start a chain of points. 
C 
        NNEW = NPT 
C 
C    Pursue the chain. 
C 
    2   NCURR = NNEW 
        MC = NCODE(NCURR) 
C 
C    Find the tile boundary by a call to TILVTX. 
C 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NCURR,VTX, 
     1    KVTMAX,KVTX,.FALSE.) 
C 
C    Check IFLAG.  Tidy up and go home if it is 10.  It ought not to be 
C    9, as we already know that NCURR is an accepted point. 
C 
        IF(IFLAG.EQ.0) GOTO 3 
        IF(IFLAG.EQ.9) STOP 763 
        CALL PLTOFF 
        DO 4  N = 1,NPTS 
        MC = NCODE(N) 
        IF(NADDR(N).GT.0.AND.MC.LT.0) NCODE(N) = -MC 
    4   CONTINUE 
        RETURN 
C 
C    Go round the tile boundary plotting in the appropriate way 
C    anything that needs to be plotted. 
C 
    3   LLO = NADDR(NCURR)+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        KPEN = 0 
        DO 5  LL = LLO,LHI 
        K = K+1 
        K1 = K+1 
        IF(K1.GT.KVTX) K1 = 1 
        N = L(LL) 
        IF(N.LE.0) GOTO 5 
        M = NCODE(N) 
        IF(M.LE.0) GOTO 5 
        NNEW = N 
        IF(M.EQ.MC) GOTO 6 
        IF(IDIFF.LT.0) GOTO 5 
        IF(IPEN.EQ.IDIFF) GOTO 7 
        IPEN = IDIFF 
        CALL PLTPEN(IPEN) 
        GOTO 7 
    6   IF(ISAME.LT.0) GOTO 5 
        IF(IPEN.EQ.ISAME) GOTO 7 
        IPEN = ISAME 
        CALL PLTPEN(IPEN) 
    7   IF(KPEN.LT.K) CALL PLTMOV(VTX(1,K),VTX(2,K)) 
        CALL PLTLIN(VTX(1,K1),VTX(2,K1)) 
        KPEN = K1 
    5   CONTINUE 
C 
C    Do the crosshatching. 
C 
        JPEN = JV(MC) 
        IF(JPEN.LT.0) GOTO 8 
        IF(IPEN.EQ.JPEN) GOTO 9 
        IPEN = JPEN 
        CALL PLTPEN(IPEN) 
    9   CALL XHATCH(VTX,KVTX,AV(MC),UV(MC),IV(MC)) 
C 
C    Plot the point if wanted. 
C 
    8   IF(IPTS.LT.0) GOTO 10 
        IF(IPEN.EQ.IPTS) GOTO 11 
        IPEN = IPTS 
        CALL PLTPEN(IPEN) 
   11   CALL PLTPT(PT(1,NCURR),PT(2,NCURR),NTAG(NCURR),IMARK) 
C 
C    Flag the code with a minus sign to show that the point is done. 
C 
   10   NCODE(NCURR) = -MC 
C 
C    Pursue the chain if it is still alive, otherwise find the next 
C    unplotted point. 
C 
        IF(NNEW.NE.NCURR) GOTO 2 
    1   CONTINUE 
C 
C    Restore the code values. 
C 
        DO 12  N = 1,NPTS 
        MC = NCODE(N) 
        IF(NADDR(N).GT.0.AND.MC.LT.0) NCODE(N) = -MC 
   12   CONTINUE 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE HATCH 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE HATCH(VTX,KVTX,ANG,Q0,Q1,P0,P1,P2,P3) 
C 
C    Hatches a polygon. 
C 
C    This subroutine assumes that the (X,Y) coordinates of 
C    the vertices of a convex polygon are held in VTX in 
C    anticlockwise order.  There are KVTX vertices, assumed 
C    to be at least 3.  Any polygon specified in this way 
C    can be hatched, but in particular this is the form of 
C    description produced from the tile database by 
C    subroutine TILVTX.  One call to the routine produces 
C    one set of hatch lines in the currently set pen colour. 
C    The remaining arguments of the subroutine specify the 
C    nature of the hatching.  ANG is the angle in radians 
C    (anticlockwise positive) between the (X,Y) axes and the 
C    (P,Q) axes.  Hatch lines are produced parallel to the 
C    P axis.  Q1 (assumed positive) is the spacing of hatch 
C    lines in the Q direction, and Q0 is the offset of some 
C    line in the raster (visible or not within the polygon) 
C    from the origin.  P3 is assumed nonnegative.  If it is 
C    zero, solid lines are produced, and P0,P1,P2 are 
C    ignored.  If it is positive, it is taken as the gap length 
C    in broken lines, with P2, assumed positive if P3 is, 
C    as the dash length.  On the line at Q, a dash starts 
C    at P0+Q*P1 and extends in the positive direction.  If 
C    any of the assumptions do not hold, no hatching is done. 
C 
        DIMENSION VTX(2,KVTX) 
        COMMON/PLTBLK/DUMMY 
C 
C    Check validity of KVTX and control values. 
C 
        IF(KVTX.LT.3) RETURN 
        IF(Q1.LE.0.0.OR.P3.LT.0.0) RETURN 
        IF(P3.GT.0.0.AND.P2.LE.0.0) RETURN 
C 
C    Set up P4 as P2+P3 if needed, and pick up the sine and 
C    cosine of the hatch angle. 
C 
        IF(P3.GT.0.0) P4 = P2+P3 
        S = SIN(ANG) 
        C = COS(ANG) 
C 
C    Find the maximum and minimum Q coordinates and identify 
C    vertices where they occur. 
C 
        KMAX = 1 
        KMIN = 1 
        QMAX = C*VTX(2,1)-S*VTX(1,1) 
        QMIN = QMAX 
        DO 1  K = 2,KVTX 
        Q = C*VTX(2,K)-S*VTX(1,K) 
        IF(Q.LE.QMAX) GOTO 2 
        KMAX = K 
        QMAX = Q 
        GOTO 1 
    2   IF(Q.GE.QMIN) GOTO 1 
        KMIN = K 
        QMIN = Q 
    1   CONTINUE 
C 
C    Locate the raster in the position prior to the first 
C    line of the hatch to be drawn, and set up a direction. 
C    To reduce pen-up movement, lines are drawn in 
C    alternate directions. 
C 
        Q = AINT((QMIN-Q0)/Q1)*Q1+Q0 
        IF(Q.GT.QMIN) Q = Q-Q1 
        IDIR = -1 
C 
C    Initialise the vertices immediately above and below 
C    the current line of the hatch in the positive and 
C    negative directions (right and left relative to 
C    the (P,Q) axes), and initialise appropriate P and 
C    Q coordinates. 
C 
        KRUP = KMIN 
        KRDN = KMIN 
        PRUP = C*VTX(1,KMIN)+S*VTX(2,KMIN) 
        QRUP = C*VTX(2,KMIN)-S*VTX(1,KMIN) 
        KLUP = KMIN 
        KLDN = KMIN 
        PLUP = PRUP 
        QLUP = QRUP 
C 
C    Scan through the lines of the hatch, alternately 
C    in either direction. 
C 
    3   Q = Q+Q1 
        IF(Q.GT.QMAX) RETURN 
        IDIR = -IDIR 
C 
C    Update the vertices above and below the current line, 
C    and their (P,Q) coordinates. 
C 
    4   IF(QRUP.GE.Q) GOTO 5 
        KRDN = KRUP 
        KRUP = KRUP+1 
        IF(KRUP.GT.KVTX) KRUP = 1 
        QRDN = QRUP 
        QRUP = C*VTX(2,KRUP)-S*VTX(1,KRUP) 
        PRDN = PRUP 
        PRUP = C*VTX(1,KRUP)+S*VTX(2,KRUP) 
        IF(KRDN.EQ.KMAX) STOP 764 
        GOTO 4 
    5   IF(QLUP.GE.Q) GOTO 6 
        KLDN = KLUP 
        KLUP = KLUP-1 
        IF(KLUP.LT.1) KLUP = KVTX 
        QLDN = QLUP 
        QLUP = C*VTX(2,KLUP)-S*VTX(1,KLUP) 
        PLDN = PLUP 
        PLUP = C*VTX(1,KLUP)+S*VTX(2,KLUP) 
        IF(KLDN.EQ.KMAX) STOP 765 
        GOTO 5 
C 
C    Find the P coordinate at each end of the 
C    current line of the hatch. 
C 
    6   PR = ((Q-QRDN)*PRUP+(QRUP-Q)*PRDN)/(QRUP-QRDN) 
        PL = ((Q-QLDN)*PLUP+(QLUP-Q)*PLDN)/(QLUP-QLDN) 
C 
C    Loop if nothing is to be drawn. 
C 
        IF(PR.LE.PL) GOTO 3 
C 
C    If a broken line is to be drawn, move on to the 
C    appropriate section of the main loop. 
C 
        IF(P3.GT.0.0) GOTO 7 
C 
C    Draw a solid line in the current direction. 
C 
        IF(IDIR.EQ.-1) GOTO 8 
C 
C    Left to right. 
C 
        CALL PLTMOV(C*PL-S*Q,S*PL+C*Q) 
        CALL PLTLIN(C*PR-S*Q,S*PR+C*Q) 
        GOTO 3 
C 
C    Right to left. 
C 
    8   CALL PLTMOV(C*PR-S*Q,S*PR+C*Q) 
        CALL PLTLIN(C*PL-S*Q,S*PL+C*Q) 
        GOTO 3 
C 
C    Draw a broken line in the current direction. 
C 
    7   POFF = P0+Q*P1 
        IF(IDIR.EQ.-1) GOTO 9 
C 
C    Left to right. 
C 
        P = AINT((PL-POFF)/P4)*P4+POFF 
        IF(P.GT.PL) P = P-P4 
   10   PA = AMAX1(P,PL) 
        PB = AMIN1(P+P2,PR) 
        IF(PA.GE.PB) GOTO 11 
        CALL PLTMOV(C*PA-S*Q,S*PA+C*Q) 
        CALL PLTLIN(C*PB-S*Q,S*PB+C*Q) 
   11   P = P+P4 
        IF(P.LT.PR) GOTO 10 
        GOTO 3 
C 
C    Right to left. 
C 
    9   P = AINT((PR-POFF)/P4)*P4+POFF 
        IF(P.GT.PR) P = P-P4 
   13   PA = AMAX1(P,PL) 
        PB = AMIN1(P+P2,PR) 
        IF(PA.GE.PB) GOTO 12 
        CALL PLTMOV(C*PB-S*Q,S*PB+C*Q) 
        CALL PLTLIN(C*PA-S*Q,S*PA+C*Q) 
   12   IF(P.LE.PL) GOTO 3 
        P = P-P4 
        GOTO 13 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE XHATCH 4.3 DATED 19 FEBRUARY 1981 
C 
        SUBROUTINE XHATCH(VTX,KVTX,A,U,I) 
C 
C    Cross-hatches a polygon in a variety of styles. 
C 
C    This subroutine assumes that the vertices of a convex 
C    polygon are held in VTX in anticlockwise order.  There are 
C    KVTX (GE.3) vertices.  A is angle of hatch to X axis, U is 
C    scale of hatch (U.LE.0.0 suppresses hatching).  I is style 
C    of hatch.  I.LE.0 or I too large suppresses hatching.  By 
C    convention use 0 for this.  Values 1 upwards give a variety 
C    of different styles. 
C 
        DIMENSION VTX(2,KVTX) 
        COMMON/PLTBLK/DUMMY 
        DATA P2,P3,P4/1.570796327,1.047197551,0.7853981634/ 
        DATA R2,R3/1.414213562,1.732050808/ 
        DATA Z,H,B,D,T,F/0.0,0.5,1.0,2.0,3.0,4.0/ 
        IF(U.LE.Z) RETURN 
        S = H*U 
        V = D*U 
        IF(I.GT.40) GOTO 41 
        IF(I.GT.20) GOTO 21 
        IF(I.GT.0) GOTO 1 
        RETURN 
C 
C    1   Simple hatch 
C 
    1   IF(I.GT.1) GOTO 2 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        RETURN 
C 
C    2   Square crosshatch 
C 
    2   IF(I.GT.2) GOTO 3 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,Z,Z) 
        RETURN 
C 
C    3   Rectangular crosshatch 
C 
    3   IF(I.GT.3) GOTO 4 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,Z,Z,Z) 
        RETURN 
C 
C    4   Long rectangular crosshatch 
C 
    4   IF(I.GT.4) GOTO 5 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,T*U,Z,Z,Z,Z) 
        RETURN 
C 
C    5   Very long rectangular crosshatch 
C 
    5   IF(I.GT.5) GOTO 6 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,F*U,Z,Z,Z,Z) 
        RETURN 
C 
C    6   Stepped rectangular crosshatch 
C 
    6   IF(I.GT.6) GOTO 7 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,B,U,U) 
        RETURN 
C 
C    7   Stepped long rectangular crosshatch 
C 
    7   IF(I.GT.7) GOTO 8 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,B,U,V) 
        RETURN 
C 
C    8   Stepped very long rectangular crosshatch 
C 
    8   IF(I.GT.8) GOTO 9 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,B,U,T*U) 
        RETURN 
C 
C    9   Bonded crosshatch 
C 
    9   IF(I.GT.9) GOTO 10 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,H,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,U,V,-H*U,H,U,U) 
        RETURN 
C 
C    10  Divided square crosshatch 
C 
   10   IF(I.GT.10) GOTO 11 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A,U,V,-U,B,V,V) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,U,V,-U,B,V,V) 
        RETURN 
C 
C    11  Herringbone crosshatch 
C 
   11   IF(I.GT.11) GOTO 12 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B,T*U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,-U,-B,T*U,U) 
        RETURN 
C 
C    12  Long herringbone crosshatch 
C 
   12   IF(I.GT.12) GOTO 13 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B,F*U,V) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,-U,-B,F*U,V) 
        RETURN 
C 
C    13  Square boxes 
C 
   13   IF(I.GT.13) GOTO 14 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        RETURN 
C 
C    14  Rectangular boxes 
C 
   14   IF(I.GT.14) GOTO 15 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,V,V) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,Z,U,U) 
        RETURN 
C 
C    15  Steps 
C 
   15   IF(I.GT.15) GOTO 16 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,-U,-B,U,U) 
        RETURN 
C 
C    16  Square waves in phase 
C 
   16   IF(I.GT.16) GOTO 17 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        RETURN 
C 
C    17  Square waves out of phase 
C 
   17   IF(I.GT.17) GOTO 18 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,H,U,U) 
        CALL HATCH(VTX,KVTX,A,U,V,S,H,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        RETURN 
C 
C    18  Square labels 
C 
   18   IF(I.GT.18) GOTO 19 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A,U,V,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        RETURN 
C 
C    19  Ladders 
C 
   19   IF(I.GT.19) GOTO 20 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        RETURN 
C 
C    20  Digit 0 
C 
   20   IF(I.GT.20) GOTO 21 
        CALL HATCH(VTX,KVTX,A,Z,T*U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A,V,T*U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,V,U) 
        RETURN 
C 
C    21  Digit 1 
C 
   21   IF(I.GT.21) GOTO 22 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,V,U) 
        RETURN 
C 
C    22  Digit 2 
C 
   22   IF(I.GT.22) GOTO 23 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,Z,U,V) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,U,Z,U,V) 
        RETURN 
C 
C    23  Digit 3 
C 
   23   IF(I.GT.23) GOTO 24 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,V,U) 
        RETURN 
C 
C    24  Digit 4 
C 
   24   IF(I.GT.24) GOTO 25 
        CALL HATCH(VTX,KVTX,A,U,T*U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,U,Z,U,V) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,V,U) 
        RETURN 
C 
C    25  Digit 5 
C 
   25   IF(I.GT.25) GOTO 26 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,U,Z,U,V) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,U,V) 
        RETURN 
C 
C    26  Digit 6 
C 
   26   IF(I.GT.26) GOTO 27 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,Z,Z,V,U) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,U,V) 
        RETURN 
C 
C    27  Digit 7 
C 
   27   IF(I.GT.27) GOTO 28 
        CALL HATCH(VTX,KVTX,A,V,T*U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,V,U) 
        RETURN 
C 
C    28  Digit 8 
C 
   28   IF(I.GT.28) GOTO 29 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,V,U) 
        RETURN 
C 
C    29  Digit 9 
C 
   29   IF(I.GT.29) GOTO 30 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,U,Z,U,V) 
        CALL HATCH(VTX,KVTX,A+P2,-U,V,Z,Z,V,U) 
        RETURN 
C 
C    30  Diamond crosshatch 
C 
   30   IF(I.GT.30) GOTO 31 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P3,Z,U,Z,Z,Z,Z) 
        RETURN 
C 
C    31  Triangular crosshatch 
C 
   31   IF(I.GT.31) GOTO 32 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P3,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+D*P3,Z,U,Z,Z,Z,Z) 
        RETURN 
C 
C    32  Triangle and hexagon crosshatch 
C 
   32   IF(I.GT.32) GOTO 33 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P3,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+D*P3,U,V,Z,Z,Z,Z) 
        RETURN 
C 
C    33  Honeycomb 
C 
   33   IF(I.GT.33) GOTO 34 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,R3,V/R3,F*U/R3) 
        CALL HATCH(VTX,KVTX,A+D*P3,Z,U,Z,R3,V/R3,F*U/R3) 
        CALL HATCH(VTX,KVTX,A+F*P3,Z,U,Z,R3,V/R3,F*U/R3) 
        RETURN 
C 
C    34  Diamond boxes 
C 
   34   IF(I.GT.34) GOTO 35 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+P3,Z,U,Z,-B/R3,V/R3,V/R3) 
        RETURN 
C 
C    35  Triangular boxes 
C 
   35   IF(I.GT.35) GOTO 36 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+P3,Z,V,Z,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+D*P3,U,V,Z,B/R3,V/R3,V/R3) 
        RETURN 
C 
C    36  Hexagonal boxes 
C 
   36   IF(I.GT.36) GOTO 37 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,Z,V/R3,F*U/R3) 
        CALL HATCH(VTX,KVTX,A+D*P3,Z,V,Z,Z,V/R3,F*U/R3) 
        CALL HATCH(VTX,KVTX,A+F*P3,U,V,R3*U,Z,V/R3,F*U/R3) 
        RETURN 
C 
C    37  Triangular labels 
C 
   37   IF(I.GT.37) GOTO 38 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P3,Z,V,Z,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+D*P3,U,V,Z,B/R3,V/R3,V/R3) 
        RETURN 
C 
C    38  Indentures 
C 
   38   IF(I.GT.38) GOTO 39 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+P3,Z,U,Z,B/R3,V/R3,V/R3) 
        RETURN 
C 
C    39  Six stars 
C 
   39   IF(I.GT.39) GOTO 40 
        CALL HATCH(VTX,KVTX,A,Z,V,-U/R3,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+P3,Z,V,-U/R3,B/R3,V/R3,V/R3) 
        CALL HATCH(VTX,KVTX,A+D*P3,Z,V,-U/R3,B/R3,V/R3,V/R3) 
        RETURN 
C 
C    40  Plus signs 
C 
   40   IF(I.GT.40) GOTO 41 
        CALL HATCH(VTX,KVTX,A,Z,V,-S,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,-S,Z,U,U) 
        RETURN 
C 
C    41  Multiplication signs 
C 
   41   IF(I.GT.41) GOTO 42 
        CALL HATCH(VTX,KVTX,A+P4,Z,R2*U,-S,B,U,R2*V-U) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,R2*U,-S,B,U,R2*V-U) 
        RETURN 
C 
C    42  Eight stars 
C 
   42   IF(I.GT.42) GOTO 43 
        CALL HATCH(VTX,KVTX,A,Z,V,-S,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P4,Z,R2*U,-S,B,U,R2*V-U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,V,-S,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,R2*U,-S,B,U,R2*V-U) 
        RETURN 
C 
C    43  Square crosshatch with one diagonal 
C 
   43   IF(I.GT.43) GOTO 44 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P4,Z,R2*S,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,Z,Z) 
        RETURN 
C 
C    44  Square crosshatch with two diagonals 
C 
   44   IF(I.GT.44) GOTO 45 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P4,Z,R2*S,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,R2*S,Z,Z,Z,Z) 
        RETURN 
C 
C    45  Square crosshatch with alternate diagonals 
C 
   45   IF(I.GT.45) GOTO 46 
        CALL HATCH(VTX,KVTX,A,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P4,Z,R2*U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,Z,Z) 
        CALL HATCH(VTX,KVTX,A+T*P4,R2*S,R2*U,Z,Z,Z,Z) 
        RETURN 
C 
C    46  Interwoven crosshatch 
C 
   46   IF(I.GT.46) GOTO 47 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,-S,B,U,U) 
        CALL HATCH(VTX,KVTX,A,Z,U,S,B,U,U) 
        RETURN 
C 
C    47  Basketwork 
C 
   47   IF(I.GT.47) GOTO 48 
        CALL HATCH(VTX,KVTX,A,Z,U,-S/D,B,U+S/D,S*T/D) 
        CALL HATCH(VTX,KVTX,A,S/D,U,S,B,U+S/D,S*T/D) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,B,U+S/D,S*T/D) 
        CALL HATCH(VTX,KVTX,A+P2,S/D,U,S*T/D,B,U+S/D,S*T/D) 
        RETURN 
C 
C    48  Conifers 
C 
   48   IF(I.GT.48) GOTO 49 
        CALL HATCH(VTX,KVTX,A+P2,Z,R2*U,R2*U,B,R2*U,R2*U) 
        CALL HATCH(VTX,KVTX,A+P4,Z,V,-S,Z,S,S*T) 
        CALL HATCH(VTX,KVTX,A+P4,-S/D,V,-T*S/D,Z,S,T*S) 
        CALL HATCH(VTX,KVTX,A+P4,-S,V,-U,Z,S,T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,V,-S,Z,S,T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,S/D,V,-T*S/D,Z,S,T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,S,V,-U,Z,S,T*S) 
        RETURN 
C 
C    49  Deciduous forest 
C 
   49   IF(I.GT.49) GOTO 50 
        CALL HATCH(VTX,KVTX,A+P4,Z,T*R2*U/D,Z,Z,R2*U/D,R2*U) 
        CALL HATCH(VTX,KVTX,A+P4,R2*U/D,T*R2*U/D,Z,Z,R2*U/D,R2*U) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,T*R2*U/D,Z,Z,R2*U/D,R2*U) 
        CALL HATCH(VTX,KVTX,A+T*P4,R2*U,T*R2*U/D,Z,Z,R2*U/D,R2*U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,T*S,-S,B,S,V+S) 
        RETURN 
C 
C    50  Mixed forest 
C 
   50   IF(I.GT.50) GOTO 51 
        CALL HATCH(VTX,KVTX,A+P2,Z,D*R2*U,R2*U,Z,R2*U,R2*U) 
        CALL HATCH(VTX,KVTX,A+P4,Z,V,-S,Z,S,S*T) 
        CALL HATCH(VTX,KVTX,A+P4,-S/D,V,-S,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+P4,-S,V,-S,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,Z,V,-S,Z,S,T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,S/D,V,-U,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,S,V,-T*U/D,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+P4,-S,V,V,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+T*P4,S,V,U,B,S,V+T*S) 
        CALL HATCH(VTX,KVTX,A+P2,R2*U,D*R2*U,R2*U/F,Z,R2*U/F,R2*(V-U/F)) 
        RETURN 
C 
C    51  Orchard 
C 
   51   IF(I.GT.51) GOTO 52 
        CALL HATCH(VTX,KVTX,A,Z,V,Z,H,U,U) 
        CALL HATCH(VTX,KVTX,A,U,V,-S,H,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,Z,U,Z,Z,U,U) 
        CALL HATCH(VTX,KVTX,A+P2,S,U,S,D,S,D*V-S) 
        RETURN 
C 
C    The crosshatch routine is user-extendable.  Begin here. 
C 
   52   RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    INTERROGATION ROUTINES LAST MODIFIED 3 APRIL 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE WINVTX 4.1 DATED 4 JULY 1979 
C 
        SUBROUTINE WINVTX(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,KVTX, 
     1    XMIN,XMAX,YMIN,YMAX) 
C 
C    Calculates the vertices of the window in clockwise cyclic order, 
C    the first vertex being that between the last and first constraints 
C    in the boundary list.  The coordinates of the Kth vertex are 
C    returned in (VTX(1,K),VTX(2,K)).  KVTX returns the number of 
C    vertices.  If this would exceed KVTMAX, which is the second 
C    dimension of the user-supplied array VTX, error return with IFLAG 
C    set to 10 occurs.  Otherwise IFLAG is zero on return.  XMIN,XMAX, 
C    YMIN,YMAX return the smallest rectangle containing the window. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP),VTX(2,KVTMAX) 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
        IF(L(2).LE.KVTMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Calculate the vertex coordinates by solving all the relevant 
C    simultaneous equations. 
C 
    1   LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        JCN = -L(LHI) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        KVTX = 0 
        DO 2  LL = LLO,LHI 
        AJOLD = AJ 
        BJOLD = BJ 
        CJOLD = CJ 
        JCN = -L(LL) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        KVTX = KVTX+1 
        D = AJOLD*BJ-AJ*BJOLD 
        VTX(1,KVTX) = (BJOLD*CJ-BJ*CJOLD)/D 
    2   VTX(2,KVTX) = (CJOLD*AJ-CJ*AJOLD)/D 
C 
C    Calculate the extreme coordinate values. 
C 
        XMIN = VTX(1,1) 
        XMAX = XMIN 
        YMIN = VTX(2,1) 
        YMAX = YMIN 
        DO 3  K = 2,KVTX 
        XMIN = AMIN1(XMIN,VTX(1,K)) 
        XMAX = AMAX1(XMAX,VTX(1,K)) 
        YMIN = AMIN1(YMIN,VTX(2,K)) 
    3   YMAX = AMAX1(YMAX,VTX(2,K)) 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILVTX 4.1 DATED 4 JULY 1979 
C 
        SUBROUTINE TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTX,KVTMAX,KVTX,RELAT) 
C 
C    Calculates the coordinates (if RELAT is .FALSE.) or coordinate 
C    offsets from NINDEX (if RELAT is .TRUE.) of the vertices of 
C    the tile of NINDEX.  The first vertex is that between the last 
C    and first neighbours of NINDEX, and the order is anticlockwise 
C    cyclic.  The coordinates, or coordinate offsets, of the Kth 
C    vertex are returned in (VTX(1,K),VTX(2,K)).  KVTX returns the 
C    number of vertices.  If this would exceed KVTMAX, which is the 
C    second dimension of the user-supplied array VTX, error return 
C    with IFLAG set to 10 occurs.  If NINDEX is not the index of an 
C    accepted point, error return with IFLAG set to 9 occurs. 
C    Otherwise IFLAG is zero on return. 
C 
C2      INTEGER*2 L 
        LOGICAL RELAT 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX) 
C 
C    Check that NINDEX is the index of an accepted point, and if so 
C    pick up LLO as the beginning of its contiguity list.  Return with 
C    IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
    2   IF(L(LLO-1).LE.KVTMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    3   LHI = LLO-1+L(LLO-1) 
        XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Find the offset of LHI, and the squared length of the offset. 
C 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LHI),U,V) 
        D = U*U+V*V 
C 
C    Loop through the neighbours of NINDEX, calculating vertex offsets 
C    as we go. 
C 
        KVTX = 0 
        DO 4  LL = LLO,LHI 
        UOLD = U 
        VOLD = V 
        DOLD = D 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LL),U,V) 
        D = U*U+V*V 
        KVTX = KVTX+1 
        C = 0.5/(UOLD*V-U*VOLD) 
        VTX(1,KVTX) = (DOLD*V-D*VOLD)*C 
    4   VTX(2,KVTX) = (D*UOLD-DOLD*U)*C 
C 
C    Modify to absolute coordinates if wanted. 
C 
        IF(RELAT) GOTO 5 
        DO 6  K = 1,KVTX 
        VTX(1,K) = VTX(1,K)+XINDEX 
    6   VTX(2,K) = VTX(2,K)+YINDEX 
C 
C    Set IFLAG to zero and return. 
C 
    5   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE OFFSET 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE OFFSET(CN,JCNS,PT,NPTS,XORIG,YORIG,NPT,UPT,VPT) 
C 
C    If NPT is positive, returns as (UPT,VPT) the offset of the 
C    position of NPT from (XORIG,YORIG).  If NPT is negative, 
C    returns as (UPT,VPT) the offset from (XORIG,YORIG) to its 
C    reflexion in constraint -NPT.  A zero value of NPT is trapped 
C    as a hard error, but otherwise no check on NPT is made. 
C 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
        IF(NPT) 1,701,2 
  701   STOP 742 
    1   JCN = -NPT 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        DJ = 2.0*(AJ*XORIG+BJ*YORIG+CN(3,JCN))/(AJ*AJ+BJ*BJ) 
        UPT = -AJ*DJ 
        VPT = -BJ*DJ 
        RETURN 
    2   UPT = PT(1,NPT)-XORIG 
        VPT = PT(2,NPT)-YORIG 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INFORM 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE INFORM(VTX,KVTX,P,AREA,X,Y,XX,XY,YY) 
C 
C    Returns size and shape parameters for a convex polygon. 
C 
C    (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX hold the vertices of a 
C    convex polygon in anticlockwise cyclic order.  P returns the 
C    perimeter, AREA the area.  For a uniform probability measure on 
C    the polygon, (X,Y) return the mean, (XX,XY,YY) the variances and 
C    covariance.  If the polygon is regarded as a thin uniform 
C    lamina of unit density, then AREA is its mass, (X,Y) its centroid 
C    (centre of mass), and (AREA*XX,AREA*XY,AREA*YY) its moments and 
C    product of inertia. 
C 
        DIMENSION VTX(2,KVTX) 
        IF(KVTX.LT.3) STOP 743 
C 
C    Calculate the perimeter. 
C 
        P = 0.0 
        U = VTX(1,KVTX) 
        V = VTX(2,KVTX) 
        DO 1  K = 1,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K) 
        V = VTX(2,K) 
        UDIFF = U-UOLD 
        VDIFF = V-VOLD 
    1   P = P+SQRT(UDIFF*UDIFF+VDIFF*VDIFF) 
C 
C    Calculate zeroth, first, and second moments about the first 
C    vertex by breaking the polygon into triangles with it as a 
C    common vertex. 
C 
        AREA = 0.0 
        X = 0.0 
        Y = 0.0 
        XX = 0.0 
        XY = 0.0 
        YY = 0.0 
        X1 = VTX(1,1) 
        Y1 = VTX(2,1) 
        U = VTX(1,2)-X1 
        V = VTX(2,2)-Y1 
        DO 2  K = 3,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K)-X1 
        V = VTX(2,K)-Y1 
C 
C    Find the area of the current triangle and add it to AREA. 
C 
        A = (UOLD*V-U*VOLD)/2.0 
        AREA = AREA+A 
C 
C    Accumulate mean relative to (X1,Y1) times area. 
C 
        X = X+A*(UOLD+U)/3.0 
        Y = Y+A*(VOLD+V)/3.0 
C 
C    Accumulate second moments about (X1,Y1) times area. 
C 
        XX = XX+A*(UOLD*UOLD+UOLD*U+U*U)/6.0 
        XY = XY+A*(UOLD*VOLD+(UOLD*V+U*VOLD)/2.0+U*V)/6.0 
    2   YY = YY+A*(VOLD*VOLD+VOLD*V+V*V)/6.0 
C 
C    Normalise, centre, and return. 
C 
        X = X/AREA 
        Y = Y/AREA 
        XX = XX/AREA-X*X 
        XY = XY/AREA-X*Y 
        YY = YY/AREA-Y*Y 
        X = X+X1 
        Y = Y+Y1 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GRAZE 4.3 DATED 3 APRIL 1981 
C 
        SUBROUTINE GRAZE(VTXOFF,KVTX,R,AREA,EATEN,ANGLE) 
C 
C    (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,...,KVTX give in anticlockwise 
C    order the coordinate offsets from a reference point of the vertices 
C    of a polygonal field, not necessarily convex or containing the 
C    reference point, although usually convex in applications. 
C    A goat is tethered to the reference point by a rope of length R. 
C    AREA returns the total area of the field.  EATEN returns the area 
C    of the grazed part of the field.  ANGLE returns the angle in 
C    radians over which the rope is an effective restraint. 
C 
        DIMENSION VTXOFF(2,KVTX) 
C 
C    Only the square of R is used.  Initialise AREA,EATEN,ANGLE. 
C 
        RSQ = R*R 
        AREA = 0.0 
        EATEN = 0.0 
        ANGLE = 0.0 
C 
C    Return immediately if the field is degenerate. 
C 
        IF(KVTX.LT.3) RETURN 
C 
C    Pick up the offset of the last vertex. 
C 
        U1 = VTXOFF(1,KVTX) 
        V1 = VTXOFF(2,KVTX) 
C 
C    Loop through the vertices, saving the offset of the previous vertex 
C    and picking up the offset of the current vertex. 
C 
C    Comments refer to the usual case where the field is starlike about 
C    the reference point, but all area and angle calculations are 
C    carried out with signed values which always accumulate correctly. 
C 
        DO 1  K = 1,KVTX 
        U0 = U1 
        V0 = V1 
        U1 = VTXOFF(1,K) 
        V1 = VTXOFF(2,K) 
C 
C    Find in terms of a barycentric coordinate T the points where the 
C    line through the previous and current vertices meets the circle of 
C    radius R centre the reference point.  First find the discriminant 
C    of the relevant quadratic.  If it is nonpositive a sector is eaten 
C    within the current triangle. 
C 
        X = U0*V1-U1*V0 
        U = U1-U0 
        V = V1-V0 
        D = U*U+V*V 
        H = RSQ*D-X*X 
        IF(H.LE.0.0) GOTO 2 
C 
C    Otherwise the roots are real and distinct.  Find them. 
C 
        H = SQRT(H) 
        P = U0*U+V0*V 
        TA = (-P-H)/D 
        TB = (-P+H)/D 
C 
C    The eaten portion of the current triangle is as follows. 
C 
C      0.0.LT.1.0.LE.TA.LT.TB  sector 
C      TA.LT.TB.LE.0.0.LT.1.0  sector 
C      0.0.LT.TA.LT.TB.LT.1.0  sector plus triangle plus sector 
C      TA.LE.0.0.LT.1.0.LE.TB  triangle 
C      TA.LE.0.0.LT.TB.LT.1.0  triangle plus sector 
C      0.0.LT.TA.LT.1.0.LE.TB  sector plus triangle 
C 
C    Test economically to find which case applies. 
C 
        IF(TB.LE.0.0) GOTO 2 
        IF(TB.GE.1.0) GOTO 3 
        IF(TA.LE.0.0) GOTO 4 
C 
C    Sector plus triangle plus sector. 
C 
        UA = U0+TA*U 
        VA = V0+TA*V 
        UB = U0+TB*U 
        VB = V0+TB*V 
        EATEN = EATEN+UA*VB-UB*VA 
        ANGLE = ANGLE+ATAN2(U0*VA-UA*V0,U0*UA+V0*VA) 
     1    +ATAN2(UB*V1-U1*VB,UB*U1+VB*V1) 
        GOTO 1 
C 
C    Triangle plus sector. 
C 
    4   UB = U0+TB*U 
        VB = V0+TB*V 
        EATEN = EATEN+U0*VB-UB*V0 
        ANGLE = ANGLE+ATAN2(UB*V1-U1*VB,UB*U1+VB*V1) 
        GOTO 1 
C 
C    More tests. 
C 
    3   IF(TA.LE.0.0) GOTO 5 
        IF(TA.GE.1.0) GOTO 2 
C 
C    Sector plus triangle. 
C 
        UA = U0+TA*U 
        VA = V0+TA*V 
        EATEN = EATEN+UA*V1-U1*VA 
        ANGLE = ANGLE+ATAN2(U0*VA-UA*V0,U0*UA+V0*VA) 
        GOTO 1 
C 
C    Triangle. 
C 
    5   EATEN = EATEN+X 
        GOTO 1 
C 
C    Sector. 
C 
    2   ANGLE = ANGLE+ATAN2(X,U0*U1+V0*V1) 
C 
C    Increment AREA. 
C 
    1   AREA = AREA+X 
C 
C    On completion of the loop, combine and normalise values and return. 
C 
        AREA = 0.5*AREA 
        EATEN = 0.5*(EATEN+RSQ*ANGLE) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDISC 4.2 DATED 19 NOVEMBER 1980 
C 
        SUBROUTINE INDISC(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTXOFF, 
     1    KVTMAX,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    For values of R in RVAL(M) for M = 1,...,MRVAL, calculates the 
C    area of the region of the window within R of accepted points, 
C    and returns it as EVAL(M), and also calculates the total angle 
C    subtended at accepted points by the curved boundary of this 
C    region, and returns it as AVAL(M).  ARETOT returns the total area 
C    of the window, and ANGTOT returns the number of accepted points 
C    times two pi.  IFLAG is zero on successful return.  The only error 
C    return is 
C 
C      10  Overflow in VTXOFF. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL),AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and rejecting if necessary. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SUBDIV 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SUBDIV(VTX,KVTX,PROP) 
C 
C    The KVTX vertices of a convex polygon are held in (either) cyclic 
C    order in (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX.  PROP(K) returns 
C    the proportion of the total area of the polygon contained within 
C    the polygon defined by the first K vertices.  PROP(1) and 
C    PROP(2) are always zero, and PROP(KVTX) is always unity. 
C 
        DIMENSION VTX(2,KVTX),PROP(KVTX) 
        IF(KVTX.LT.3) STOP 744 
        PROP(1) = 0.0 
        PROP(2) = 0.0 
        A = 0.0 
        X1 = VTX(1,1) 
        Y1 = VTX(2,1) 
        U = VTX(1,2)-X1 
        V = VTX(2,2)-Y1 
        DO 1  K = 3,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K)-X1 
        V = VTX(2,K)-Y1 
        A = A+UOLD*V-U*VOLD 
    1   PROP(K) = A 
        DO 2  K = 3,KVTX 
    2   PROP(K) = PROP(K)/A 
        PROP(KVTX) = 1.0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDICO 4.2 DATED 19 NOVEMBER 1980 
C 
        SUBROUTINE INDICO(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NCODE, 
     1    VTXOFF,KVTMAX,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT,NC) 
C 
C    Carries out the same calculation as subroutine INDISC, but 
C    accumulates contributions only from those points whose code 
C    values as looked up in NCODE are equal to NC.  ARETOT returns 
C    the total tile area associated with such points, and ANGTOT 
C    returns the number of such points times two pi. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP),NCODE(NPTS), 
     1    VTXOFF(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL),AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and rejecting if necessary. 
C 
        DO 2  NINDEX = 1,NPTS 
        IF(NCODE(NINDEX).NE.NC) GOTO 2 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE WINBD 4.1 DATED 27 SEPTEMBER 1979 
C 
        SUBROUTINE WINBD(CN,JCNS,L,LTOP,XMIN,XMAX,YMIN,YMAX) 
C 
C    XMIN,XMAX,YMIN,YMAX return the smallest rectangle in the coordinate 
C    directions which contains the window. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP) 
C 
C    Calculate the vertex coordinates, updating extreme values. 
C 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        JCN = -L(LHI) 
        AJOLD = CN(1,JCN) 
        BJOLD = CN(2,JCN) 
        CJOLD = CN(3,JCN) 
        JCN = -L(LLO) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        D = AJOLD*BJ-AJ*BJOLD 
        XMIN = (BJOLD*CJ-BJ*CJOLD)/D 
        YMIN = (CJOLD*AJ-CJ*AJOLD)/D 
        XMAX = XMIN 
        YMAX = YMIN 
        LLO = LLO+1 
        DO 1  LL = LLO,LHI 
        AJOLD = AJ 
        BJOLD = BJ 
        CJOLD = CJ 
        JCN = -L(LL) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        D = AJOLD*BJ-AJ*BJOLD 
        X = (BJOLD*CJ-BJ*CJOLD)/D 
        Y = (CJOLD*AJ-CJ*AJOLD)/D 
        XMIN = AMIN1(XMIN,X) 
        XMAX = AMAX1(XMAX,X) 
        YMIN = AMIN1(YMIN,Y) 
    1   YMAX = AMAX1(YMAX,Y) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILVTE 4.1 DATED 21 OCTOBER 1980 
C 
        SUBROUTINE TILVTE(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTX,KVTMAX,KVTX,RELAT,IEDGE) 
C 
C    Calculates the coordinates (if RELAT is .FALSE.) or coordinate 
C    offsets from NINDEX (if RELAT is .TRUE.) of the vertices of 
C    the tile of NINDEX.  The first vertex is that between the last 
C    and first neighbours of NINDEX, and the order is anticlockwise 
C    cyclic.  The coordinates, or coordinate offsets, of the Kth 
C    vertex are returned in (VTX(1,K),VTX(2,K)).  KVTX returns the 
C    number of vertices.  If this would exceed KVTMAX, which is the 
C    second dimension of the user-supplied array VTX, error return 
C    with IFLAG set to 10 occurs.  If NINDEX is not the index of an 
C    accepted point, error return with IFLAG set to 9 occurs. 
C    Otherwise IFLAG is zero on return.  Thus this routine does 
C    exactly the same as TILVTX, which it closely resembles.  However 
C    it returns one additional argument IEDGE.  This returns -1 if 
C    NINDEX is contiguous to some constraint, and 0 or 1 if not.  The 
C    value 1 is returned if the tile of NINDEX could be affected by 
C    points outside the current window if it were enlarged, and 
C    otherwise 0 is returned. 
C 
C2      INTEGER*2 L 
        LOGICAL RELAT 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX) 
C 
C    Check that NINDEX is the index of an accepted point, and if so 
C    pick up LLO as the beginning of its contiguity list.  Return with 
C    IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
    2   IF(L(LLO-1).LE.KVTMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    3   LHI = LLO-1+L(LLO-1) 
        XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Initialise IEDGE to zero. 
C 
        IEDGE = 0 
C 
C    Find the offset of LHI, and the squared length of the offset. 
C 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LHI),U,V) 
        D = U*U+V*V 
C 
C    Loop through the neighbours of NINDEX, calculating vertex offsets 
C    as we go, and resetting IEDGE to -1 if any neighbour is a 
C    constraint. 
C 
        KVTX = 0 
        DO 4  LL = LLO,LHI 
        UOLD = U 
        VOLD = V 
        DOLD = D 
        N = L(LL) 
        IF(N.LT.0) IEDGE = -1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,N,U,V) 
        D = U*U+V*V 
        KVTX = KVTX+1 
        C = 0.5/(UOLD*V-U*VOLD) 
        VTX(1,KVTX) = (DOLD*V-D*VOLD)*C 
    4   VTX(2,KVTX) = (D*UOLD-DOLD*U)*C 
C 
C    Unless IEDGE is already -1, scan vertices against constraints.  If 
C    any vertex is nearer to a constraint than to NINDEX, reset IEDGE 
C    to 1 and drop out of the scan. 
C 
        IF(IEDGE.EQ.-1) GOTO 7 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        DO 5  K = 1,KVTX 
        UVTX = VTX(1,K) 
        VVTX = VTX(2,K) 
        DVTX = SQRT(UVTX**2+VVTX**2) 
        XVTX = XINDEX+UVTX 
        YVTX = YINDEX+VVTX 
        DO 5  LL = LLO,LHI 
        J = -L(LL) 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        D = -(AJ*XVTX+BJ*YVTX+CN(3,J))/SQRT(AJ**2+BJ**2) 
        IF(D.LT.DVTX) GOTO 6 
    5   CONTINUE 
        GOTO 7 
    6   IEDGE = 1 
C 
C    Modify to absolute coordinates if wanted. 
C 
    7   IF(RELAT) GOTO 8 
        DO 9  K = 1,KVTX 
        VTX(1,K) = VTX(1,K)+XINDEX 
    9   VTX(2,K) = VTX(2,K)+YINDEX 
C 
C    Set IFLAG to zero and return. 
C 
    8   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE MEAD 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE MEAD(VTXOFF,KVTX,AREA,ECCIR,ABCEN) 
C 
C    Calculates tile area, eccircularity, and abcentricity as defined 
C    by R. Mead, "A relationship between individual plant spacing and 
C    yield", Annals of Botany (New Series) 30 (1966) pp. 301-309. 
C 
C    The statistics are calculated for a tile the coordinates offsets 
C    of whose KVTX vertices are held in anticlockwise cyclic order 
C    in (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,...,KVTX.  These values 
C    will usually have been set up by a call to TILVTX with RELAT 
C    set to .TRUE..  AREA, ECCIR, ABCEN return the area, 
C    eccircularity, and abcentricity respectively. 
C 
        DIMENSION VTXOFF(2,KVTX) 
        IF(KVTX.LT.3) STOP 745 
C 
C    Calculate the area of the tile and the offset of its centroid. 
C 
        UCURR = VTXOFF(1,KVTX) 
        VCURR = VTXOFF(2,KVTX) 
        AREA = 0.0 
        UCENTR = 0.0 
        VCENTR = 0.0 
        DO 1  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = VTXOFF(1,K) 
        VCURR = VTXOFF(2,K) 
        A = 0.5*(UOLD*VCURR-UCURR*VOLD) 
        AREA = AREA+A 
        UCENTR = UCENTR+A*(UOLD+UCURR)/3.0 
    1   VCENTR = VCENTR+A*(VOLD+VCURR)/3.0 
        UCENTR = UCENTR/AREA 
        VCENTR = VCENTR/AREA 
        OFFCEN = SQRT(UCENTR**2+VCENTR**2) 
C 
C    Calculate the weighted sum of centroid-to-vertex distances, using 
C    the exterior angles at the vertices as weights; they must in fact 
C    add up to two pi, but we obtain the normalising constant by 
C    accumulating them. 
C 
        WTSUM = 0.0 
        DISSUM = 0.0 
        UNEW = VTXOFF(1,KVTX) 
        VNEW = VTXOFF(2,KVTX) 
        UCURR = VTXOFF(1,KVTX-1) 
        VCURR = VTXOFF(2,KVTX-1) 
        DO 2  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = UNEW 
        VCURR = VNEW 
        UNEW = VTXOFF(1,K) 
        VNEW = VTXOFF(2,K) 
        DIST = SQRT((UCURR-UCENTR)**2+(VCURR-VCENTR)**2) 
        PN = UNEW-UCURR 
        QN = VNEW-VCURR 
        PO = UOLD-UCURR 
        QO = VOLD-VCURR 
        ANG = ATAN2(PN*QO-PO*QN,-PO*PN-QO*QN) 
        WTSUM = WTSUM+ANG 
    2   DISSUM = DISSUM+ANG*DIST 
        DISSUM = DISSUM/WTSUM 
C 
C    Calculate the statistics. 
C 
        ECCIR = DISSUM*SQRT(WTSUM*0.5/AREA) 
        ABCEN = OFFCEN/DISSUM 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NLIST 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE NLIST(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    NBR,PTOFF,KNBMAX,KNB) 
C 
C    Returns in NBR(K) the Kth neighbour of NINDEX and in (PTOFF(1,K), 
C    PTOFF(2,K)) the offset to it as a real or virtual point for 
C    K = 1,...,KNB, where KNB returns the number of neighbours. 
C    IFLAG returns 0 on normal return, 9 if NINDEX is not the index 
C    of an accepted point, and 10 if KNB exceeds KNBMAX. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    NBR(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check NINDEX. 
C 
        IF(NINDEX.LT.1.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN) 1,701,2 
  701   STOP 746 
    1   IFLAG = 9 
        RETURN 
C 
C    Copy the contiguity list after checking its size. 
C 
    2   KNB = L(LOCN+1) 
        IF(KNB.LE.KNBMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
    3   K = 0 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        DO 4  LL = LLO,LHI 
        K = K+1 
    4   NBR(K) = L(LL) 
C 
C    Calculate and save the offsets. 
C 
        XORIG = PT(1,NINDEX) 
        YORIG = PT(2,NINDEX) 
        DO 5  K = 1,KNB 
        N = NBR(K) 
        CALL OFFSET(CN,JCNS,PT,NPTS,XORIG,YORIG,N,U,V) 
        PTOFF(1,K) = U 
    5   PTOFF(2,K) = V 
C 
C    Set IFLAG to 0 and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE AREAS 4.1 DATED 27 OCTOBER 1980 
C 
        SUBROUTINE AREAS(VTXOFF,KVTX,AIN,AREA,AOUT) 
C 
C    Calculates tile area together with the area of the largest disc 
C    centred at the generating point and contained within the tile and 
C    the area of the smallest disc centred at the generating point 
C    and containing the tile.  These values are returned as AREA, AIN, 
C    and AOUT respectively.  (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,..., 
C    KVTX hold the vertex offsets from the generating point for the 
C    KVTX vertices of the tile in anticlockwise cyclic order.  These 
C    values will usually have been loaded by a call to TILVTX or to 
C    TILVTE. 
C 
        DIMENSION VTXOFF(2,KVTX) 
C 
C    The calculations are accomplished trivially in one scan. 
C 
        UCURR = VTXOFF(1,KVTX) 
        VCURR = VTXOFF(2,KVTX) 
        AREA = 0.0 
        AOUT = 0.0 
        AIN = UCURR**2+VCURR**2 
        DO 1  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = VTXOFF(1,K) 
        VCURR = VTXOFF(2,K) 
        A2 = UOLD*VCURR-UCURR*VOLD 
        AREA = AREA+A2 
        AOUT = AMAX1(AOUT,UCURR**2+VCURR**2) 
    1   AIN = AMIN1(AIN,(A2**2)/((UCURR-UOLD)**2+(VCURR-VOLD)**2)) 
        AREA = 0.5*AREA 
        AOUT = 3.141593*AOUT 
        AIN = 3.141593*AIN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRIM 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE TRIM(CN,JCNS,L,LTOP,EPSCN,IFLAG,XPT,YPT,VTXOFF, 
     1    KVTMAX,KVTX,ITRIM,WORK) 
C 
C    The array VTXOFF(2,KVTMAX) holds on entry the vertex offsets from 
C    (XPT,YPT) of a polygon, usually the tile of that point, in 
C    anticlockwise order.  This is modified by TRIM to give the vertex 
C    offsets of the intersection of it with a reduced window moved in 
C    by EPSCN from the original.  The original number of vertices is 
C    KVTX.  This is modified to the new number, or conventionally set 
C    to 1 if the trimmed polygon has empty interior.  IFLAG is zero on 
C    successful return, and returns 10 if at any stage in the 
C    calculation KVTX would exceed KVTMAX, in which case the information 
C    in VTXOFF is lost.  ITRIM returns -1,0,1 as (XPT,YPT) is strictly 
C    within, on, or strictly outside the reduced window.  The array 
C    WORK(2,KVTMAX) must be supplied as workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP),VTXOFF(2,KVTMAX),WORK(2,KVTMAX) 
C 
C    Initialise ITRIM.  Convert offsets to absolute coordinates. 
C 
        ITRIM = -1 
        DO 1  K = 1,KVTX 
        VTXOFF(1,K) = VTXOFF(1,K)+XPT 
    1   VTXOFF(2,K) = VTXOFF(2,K)+YPT 
C 
C    Chop the polygon by each effective constraint in turn, picking them 
C    up from the boundary list in the heap. 
C 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        DO 2  LL = LLO,LHI 
        J = -L(LL) 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        CJ = CN(3,J)+EPSCN*SQRT(AJ**2+BJ**2) 
C 
C    Test (XPT,YPT) and reset ITRIM if necessary. 
C 
        W = AJ*XPT+BJ*YPT+CJ 
        IF(W.EQ.0.0) ITRIM = MAX0(ITRIM,0) 
        IF(W.GT.0.0) ITRIM = 1 
C 
C    Quick scan through the vertices to avoid unnecessary work if the 
C    chop either misses or annihilates the polygon. 
C 
        IF(KVTX.LT.3) GOTO 2 
        KIN = 0 
        KOUT = 0 
        DO 3  K = 1,KVTX 
        IF(AJ*VTXOFF(1,K)+BJ*VTXOFF(2,K)+CJ) 4,3,5 
    4   KIN = KIN+1 
        GOTO 3 
    5   KOUT = KOUT+1 
    3   CONTINUE 
        IF(KOUT.EQ.0) GOTO 2 
        IF(KIN.GT.0) GOTO 6 
        KVTX = 1 
        GOTO 2 
C 
C    Properly chopped case.  Copy to WORK, modifying as we go and 
C    checking space.  Then copy back to VTXOFF. 
C 
    6   KCOPY = 0 
        X = VTXOFF(1,KVTX) 
        Y = VTXOFF(2,KVTX) 
        W = AJ*X+BJ*Y+CJ 
        DO 7  K = 1,KVTX 
        XOLD = X 
        YOLD = Y 
        WOLD = W 
        X = VTXOFF(1,K) 
        Y = VTXOFF(2,K) 
        W = AJ*X+BJ*Y+CJ 
        IF(W.GT.0.0) GOTO 8 
        IF(W.EQ.0.0.OR.WOLD.LE.0.0) GOTO 9 
        IF(KCOPY.LT.KVTMAX) GOTO 10 
        IFLAG = 10 
        RETURN 
   10   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = (WOLD*X-W*XOLD)/(WOLD-W) 
        WORK(2,KCOPY) = (WOLD*Y-W*YOLD)/(WOLD-W) 
    9   IF(KCOPY.LT.KVTMAX) GOTO 11 
        IFLAG = 10 
        RETURN 
   11   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = X 
        WORK(2,KCOPY) = Y 
        GOTO 7 
    8   IF(WOLD.GE.0.0) GOTO 7 
        IF(KCOPY.LT.KVTMAX) GOTO 12 
        IFLAG = 10 
        RETURN 
   12   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = (WOLD*X-W*XOLD)/(WOLD-W) 
        WORK(2,KCOPY) = (WOLD*Y-W*YOLD)/(WOLD-W) 
    7   CONTINUE 
        KVTX = KCOPY 
        DO 13  K = 1,KVTX 
        VTXOFF(1,K) = WORK(1,K) 
   13   VTXOFF(2,K) = WORK(2,K) 
    2   CONTINUE 
C 
C    Return to relative coordinates. 
C 
        DO 14  K = 1,KVTX 
        VTXOFF(1,K) = VTXOFF(1,K)-XPT 
   14   VTXOFF(2,K) = VTXOFF(2,K)-YPT 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDITR 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE INDITR(CN,JCNS,PT,NPTS,NADDR,L,LTOP,EPSCN,IFLAG, 
     1    VTXOFF,KVTMAX,WORK,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    Carries out the same calculation as INDISC, but trims all 
C    tiles by moving the window boundary in by amount EPSCN. 
C    ARETOT returns the area within the trimmed window, and 
C    ANGTOT returns two pi times the number of points within 
C    the trimmed window.  WORK(2,KVTMAX) must be supplied as 
C    workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),WORK(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL), 
     2    AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and checking IFLAG. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Trim the tile. 
C 
        CALL TRIM(CN,JCNS,L,LTOP,EPSCN,IFLAG,PT(1,NINDEX), 
     1    PT(2,NINDEX),VTXOFF,KVTMAX,KVTX,ITRIM,WORK) 
        IF(IFLAG.EQ.10) RETURN 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        IF(ITRIM.EQ.-1) ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDIVT 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE INDIVT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG, 
     1    VTXOFF,KVTMAX,WORK,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    Carries out the same calculation as INDISC, but trims all 
C    tiles by moving the window boundary in by amount RVAL(M) 
C    when making the calculations for that value.  This is done 
C    progressively, so for this routine the values in RVAL must 
C    be in ascending order.  ARETOT(M) returns the area within 
C    the trimmed window, and ANGTOT(M) returns the number of 
C    points within the trimmed window times two pi.  Note that 
C    both these quantities depend on the values in RVAL, and 
C    so in this routine ARETOT and ANGTOT have to be arrays. 
C    WORK(2,KVTMAX) must be supplied as workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),WORK(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL), 
     2    AVAL(MRVAL),ARETOT(MRVAL),ANGTOT(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
        AVAL(M) = 0.0 
        ARETOT(M) = 0.0 
    1   ANGTOT(M) = 0.0 
C 
C    Scan points, finding the tile for each and checking IFLAG. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
C 
C    Trim the tile. 
C 
        CALL TRIM(CN,JCNS,L,LTOP,RVAL(M),IFLAG,PT(1,NINDEX), 
     1    PT(2,NINDEX),VTXOFF,KVTMAX,KVTX,ITRIM,WORK) 
        IF(IFLAG.EQ.10) RETURN 
C 
C    Calculate the quantities. 
C 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
        AVAL(M) = AVAL(M)+ANGLE 
        ARETOT(M) = ARETOT(M)+AREA 
    3   IF(ITRIM.EQ.-1) ANGTOT(M) = ANGTOT(M)+6.2831853 
    2   CONTINUE 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
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
C    DUMMY PLOT INTERFACE ROUTINES
C
      SUBROUTINE PLTON(XMIN,XMAX,YMIN,YMAX)
      RETURN
      END
C
      SUBROUTINE PLTMOV(X,Y)
      RETURN
      END
C
      SUBROUTINE PLTLIN(X,Y)
      RETURN
      END
C
      SUBROUTINE PLTPEN(IPEN)
      RETURN
      END
C
      SUBROUTINE PLTPT(X,Y,LABEL,IMARK)
C
C    PLOTS POINT AT (X,Y), ANNOTATING IT AS FOLLOWS:
C    IMARK = 1  NO CROSS, INTEGER LABEL
C            2  NO CROSS, A4 LABEL
C            3  CROSS
C            4  CROSS, INTEGER LABEL
C            5  CROSS, A4 LABEL
C
      RETURN
      END
C
      SUBROUTINE PLTOFF
      RETURN
      END
C*********************************************************************** 
C
C(2) THE FOLLOWING DUMMY ROUTINES ARE INCLUDED
C    SIMPLY TO SATISFY EXTERNAL REFERENCES IN SOLID
C
      FUNCTION ZED(M,N)
      ZED = 0.0
      RETURN
      END
C
      FUNCTION ZEDX(M,N)
      ZEDX = 0.0
      RETURN
      END
C
      FUNCTION ZEDY(M,N)
      ZEDY = 0.0
      RETURN
      END
C*********************************************************************** 

