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
