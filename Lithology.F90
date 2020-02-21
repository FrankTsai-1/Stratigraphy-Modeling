Program Geologic_Modeling
Implicit Double Precision (A-H,O-Z)

!------------------------------------------------------------------------------------------------------
!GridTopo.csv   : input file, center of cell locations (x,y) and its land surface (ft). Note: at this point, only 500mx500m squared cells are considered. 
!WellLogs.csv   : input file, well log data file.
!------------------------------------------------------------------------------------------------------
!MaxLithoBed    :maximum number of lithology beds in a well log 
!------------------------------------------------------------------------------------------------------
!dipDirection   :dip direction in degrees: 0 deg: East, 90 deg: North, 180 deg: West, 270 deg: South
!dipAngleS      :start dip angle (degrees) !piece-wise line function for dip angle
!dipAngleE      :end dip angle (degrees) 
!ElevP1         :elevation of the first dip angle deflection (feet)
!ElevP2         :elevation of the second dip angle dfelction (feet)
!StartElev      :start elevation (ft)
!EndElev        :bottom elevation (ft)
!DZ             :interval DZ (ft)
!ThickMin       :minimum thickness of a facies (ft)
!cutoff         :for kriging indicator (-) 
!NMinLogs       :minimum number of well logs for indicator kriging (-)
!SearchRange    :Additional range outside the modeing domain (ft)
!------------------------------------------------------------------------------------------------------
Character*80 str1

namelist / geo   / NCELL,NWellLog,MaxLithoBed, &
                   dipDirection,dipAngleS,dipAngleE,ElevP1,ElevP2,StartElev,EndElev, &
                   DZ,ThickMin,cutoff,NMinLogs,SearchRange
Integer :: error

Character*15, Allocatable :: str(:)
Allocatable :: WellLogX(:),WellLogY(:)
Allocatable :: XCELL(:),YCELL(:),TOPO(:)
Real*4, Allocatable :: VCELL(:)
Allocatable :: NLithoBed(:), ElevBedBoundary(:,:), IndexBedLitho(:,:)
Allocatable :: XX(:),YY(:),VV(:),iActiveWell(:)
Allocatable :: indx(:),BB(:)
Allocatable :: Ind(:)

Allocatable :: z2(:,:),I2(:,:),V2(:,:),z3(:,:)

Call cpu_time (program_start)

Open(unit=999,file="geomodel.prn")

!Greetings
write (*,*)
write (*,*) "                            GEOMODEL-Tile                          " 
write (*,*) "                     LOUISIANA STATE UNIVERSITY                    "
write (*,*) "                       Version 1.00 11/26/2017                     "
write (*,*)
write (*,*)

write (999,*) "                            GEOMODEL-Tile                          " 
write (999,*) "                     LOUISIANA STATE UNIVERSITY                    "
write (999,*) "                       Version 1.05 11/26/2017                     "
write (999,*)
!read geo model input parameters
OPEN (UNIT=23, FILE='Geo_model.inp', STATUS='OLD')
READ (23, NML = geo)
CLOSE (23)

write(999,*) "Parameters:"
write(999,*) "NCELL= ",NCELL
write(999,*) "NWellLog= ",NWellLog
write(999,*) "MaxLithoBed= ",MaxLithoBed
write(999,*) "dipDirection= ",dipDirection
write(999,*) "dipAngleS= ",dipAngleS
write(999,*) "dipAngleE= ",dipAngleE
write(999,*) "ElevP1= ",ElevP1
write(999,*) "ElevP2= ",ElevP2
write(999,*) "StartElev= ",StartElev
write(999,*) "EndElev= ",EndElev
write(999,*) "DZ= ",DZ
write(999,*) "ThickMin= ",ThickMin
write(999,*) "cutoff= ",cutoff
write(999,*) "NMinLogs= ",NMinLogs
write(999,*)

if (dipAngleS>dipAngleE) then
    write(*,*) "dipAngleS>dipAngleE"
    write(999,*) "dipAngleS>dipAngleE"
    stop
endif

if (ElevP1<ElevP2) then 
    write(*,*) "ElevP1<ElevP2"
    write(999,*) "ElevP1<ElevP2"
    stop
endif

if (StartElev<EndElev) then
    write(*,*) "StartElev<EndElev"
    write(999,*) "StartElev<EndElev"
    stop
endif

! Read 2D surficial cell location and their elevation
Allocate(XCELL(NCELL),stat=error)
Allocate(YCELL(NCELL),stat=error)
Allocate(TOPO(NCELL),stat=error)
!
Call Grid2D(NCELL,XCELL,YCELL,TOPO)

! Domain of active well logs for interpolation
XSearchMin=Minval(XCELL(1:NCELL))-SearchRange
XSearchMax=Maxval(XCELL(1:NCELL))+SearchRange
YSearchMin=Minval(YCELL(1:NCELL))-SearchRange
YSearchMax=Maxval(YCELL(1:NCELL))+SearchRange

write(999,*) "GridXMin", Minval(XCELL(1:NCELL))
write(999,*) "GridXMax", Maxval(XCELL(1:NCELL))
write(999,*) "GridYMin", Minval(YCELL(1:NCELL))
write(999,*) "GridYMax", Maxval(YCELL(1:NCELL))
write(999,*) "TopoZMin", Minval(TOPO(1:NCELL))
write(999,*) "TopoZMax", Maxval(TOPO(1:NCELL))

! Derive coefficients A,B,C to calcualte distance=A*X+B*Y+C, a distance from a point (X,Y) to the strik (AX+BY+C=0) with a pivot (MinXC,MinYC) 
! Pivot point is Minval(XCELL) and Minval(YCELL)
pi=3.14159265358979323846264338327950288419716939937510d0
A=dcos(dipDirection/180d0*pi)
B=dsin(dipDirection/180d0*pi)
MinXC=Minval(XCELL(1:NCELL)) !Pivot X, at lower left corner
MinYC=Minval(YCELL(1:NCELL)) !Pivot Y, at lower left corner
C=-(A*MinXC+B*MinYC) 

! Calculate the distant to the strike, where two dipped planes intercept.  
Dist1=A*XSearchMin+B*YSearchMin+C !distance from search domain corner to the strike
Dist2=A*XSearchMin+B*YSearchMax+C !distance from search domain corner to the strike
Dist3=A*XSearchMax+B*YSearchMin+C !distance from search domain corner to the strike
Dist4=A*XSearchMax+B*YSearchMax+C !distance from search domain corner to the strike
Distmax=dmin1(Dist1,Dist2,Dist3,Dist4) !max distance from updip region in the search domain
ENew1=ElevP1-Distmax*dipAngleS*3.2808398950131234  !convert meter to foot
ENew2=ElevP2-Distmax*dipAngleE*3.2808398950131234  !convert meter to foot
if (ENew1<ENew2) then ! Program stops if lower dipped plane is higher than the upper dipped plane (two plants intercept within the search domain)
    write(*,*) "The transformed bed elevations inside the search domain may be overturned."    
    write(*,*)    
    dipAngleS1=(ElevP1-ENew2)/(Distmax*3.2808398950131234)
    write(*,'(A31,F10.6,A10)') "Increase dipAngleS greater than", abs(dipAngleS1), "degrees"
    dipAngleE1=(ElevP2-ENew1)/(Distmax*3.2808398950131234)
    write(*,'(A28,F10.6,A10)') "Decrease dipAngleE less than", abs(dipAngleE1), "degrees"       
    ElevP11=ElevP1+(ENew2-ENew1)
    write(*,'(A27,F10.1,A5)') "Increase ElevP1 higher than", ElevP11, "ft" 
    ElevP21=ElevP2-(ENew2-ENew1)
    write(*,*) "OR "
    write(*,'(A26,F10.1,A5)') "Decrease ElevP2 lower than", ElevP21, "ft"       
    write(*,*)
    
    write(999,*) "The transformed bed elevations inside the search domain may be overturned."    
    write(999,*)    
    dipAngleS1=(ElevP1-ENew2)/(Distmax*3.2808398950131234)
    write(999,'(A31,F10.6,A10)') "Increase dipAngleS greater than", abs(dipAngleS1), "degrees"
    dipAngleE1=(ElevP2-ENew1)/(Distmax*3.2808398950131234)
    write(999,'(A28,F10.6,A10)') "Decrease dipAngleE less than", abs(dipAngleE1), "degrees"       
    ElevP11=ElevP1+(ENew2-ENew1)
    write(999,'(A27,F10.1,A5)') "Increase ElevP1 higher than", ElevP11, "ft" 
    ElevP21=ElevP2-(ENew2-ENew1)
    write(999,*) "OR "
    write(999,'(A26,F10.1,A5)') "Decrease ElevP2 lower than", ElevP21, "ft"       
    write(999,*)
    stop
endif

! Calculate distance and elevation to the center strike
if (dipAngleS.eq.dipAngleE) then 
    DistCenter=1.d+10 !dummy value
    ElevCenter=1.d+10 !dummy value 
    C_center=1.d+10 !dummy value
else
    DistCenter=(ElevP1-ElevP2)/((dipAngleS-dipAngleE)*3.2808398950131234) !distance from Center strike to Pivot strike
    ElevCenter=ElevP1-DistCenter*dipAngleS*3.2808398950131234 !elevation of the Center strike
    C_center=C-DistCenter ! C value for Center strike equation AX+BY+C_center
endif

! Read all well log data and select active well logs in the domain
!
Allocate(WellLogX(NWellLog),stat=error)
Allocate(WellLogY(NWellLog),stat=error)
Allocate(NLithoBed(NWellLog),stat=error)
Allocate(ElevBedBoundary(MaxLithoBed,NWellLog),stat=error)
Allocate(IndexBedLitho(MaxLithoBed,NWellLog),stat=error)

NBed=MaxLithoBed
Call WellLogData(XSearchMin,XSearchMax,YSearchMin,YSearchMax,NWellLog,NBed,NLithoBed,WellLogX,WellLogY,ElevBedBoundary,IndexBedLitho)
write(999,*) "Maximum number of beds of boring wells is ", Maxval(NLithoBed(1:NWellLog))+1

write(*,*) "Active wells: ", NWellLog
write(999,*) "Active wells: ", NWellLog

write(999,*) "WellXMin", Minval(WellLogX(1:NWellLog))
write(999,*) "WellXMax", Maxval(WellLogX(1:NWellLog))
write(999,*) "WellYMin", Minval(WellLogY(1:NWellLog))
write(999,*) "WEllYMax", Maxval(WellLogY(1:NWellLog))

! This is for NN interpolation
XMIN2= MIN(Minval(WellLogX(1:NWellLog)),Minval(XCELL(1:NCELL))) -1.d0
XMAX2= MAX(Maxval(WellLogX(1:NWellLog)),Maxval(XCELL(1:NCELL))) +1.d0
YMIN2= MIN(Minval(WellLogY(1:NWellLog)),Minval(YCELL(1:NCELL))) -1.d0
YMAX2= MAX(Maxval(WellLogY(1:NWellLog)),Maxval(YCELL(1:NCELL))) +1.d0

write(999,*) "XMIN2= ",XMIN2
write(999,*) "XMAX2= ",XMAX2
write(999,*) "YMIN2= ",YMIN2
write(999,*) "YMIN2= ",YMAX2

!Consider Dip and re-calculate elevation of bed boudnaries at horizontal plane 
do 15 jj=1,NWellLog
Dist1=A*WellLogX(jj)+B*WellLogY(jj)+C_center !distance to center strike
do 16 i=1,NLithoBed(jj)
    z=ElevBedBoundary(i,jj)
    dipAngle1=dipAngle(ElevCenter,dipAngleS,dipAngleE,Dist1,z)
    Dist=A*WellLogX(jj)+B*WellLogY(jj)+C !distrance to pivot strike can be positive or negative depending on its location at downdip or updip, respectively.     
    E1=Dist*dipAngle1*3.2808398950131234
    znew=z+E1 ! transform original bed elevations to non-dip domain
    ElevBedBoundary(i,jj)=znew
    if (i.ne.1) then
        if (ElevBedBoundary(i,jj).gt.ElevBedBoundary(i-1,jj)) then 
            write(*,*) "transformed bed elevation is over turned: ",ElevBedBoundary(i-1,jj),ElevBedBoundary(i,jj)
            write(999,*) "transformed bed elevation is over turned: ",ElevBedBoundary(i-1,jj),ElevBedBoundary(i,jj)
            stop
        endif
    endif
    
16 enddo
15 enddo

! find max distance of downdip cells to the pivot strike
Distmax=-1.0D+10 !a dummy value
do i=1,NCELL
    Dist=A*XCELL(i)+B*YCELL(i)+C
    if (Dist>Distmax) then 
        Distmax=Dist
        XMax=XCELL(i)
        YMax=YCELL(i)
    endif
enddo
if (Distmax>0) then 
    Dist1=A*XMax+B*YMax+C_center !distance to center strike
    dipAngle1=dipAngle(ElevCenter,dipAngleS,dipAngleE,Dist1,StartElev)   
    Dist=A*XMax+B*YMax+C
    StartElev1=StartElev+Dist*dipAngle1*3.2808398950131234 !for additional elevation
else
    StartElev1=StartElev
endif

write(999,*)
write(999,*) "Additional depth (ft) from ",StartElev, "to ",StartElev1

! find max distance of updip cells to the pivot strike
Distmax=1.0D+10 !a dummy value
do i=1,NCELL
    Dist=A*XCELL(i)+B*YCELL(i)+C
    if (Dist<Distmax) then 
        Distmax=Dist
        XMax=XCELL(i)
        YMax=YCELL(i)
    endif
enddo
if (Distmax<0) then 
    Dist1=A*XMax+B*YMax+C_center !distance to center strike
    dipAngle1=dipAngle(ElevCenter,dipAngleS,dipAngleE,Dist1,EndElev)   
    Dist=A*XMax+B*YMax+C
    EndElev1=EndElev+Dist*dipAngle1*3.2808398950131234 !for additional elevation
else
    EndElev1=EndElev
endif

write(999,*) "Additional depth (ft) from ",EndElev, "to ",EndElev1
write(999,*)

!Active well logs for interpolation
NLay=NINT((StartElev1-EndElev1)/DZ)+1 !number of beds 
!
write(999,*) "Number of planes= ",NLay
write(999,*)
!
! Active wells for each plane
write (999,*) "Plane No.   Elevation    Number of active wells"
do iLay=1,NLAY
z=StartElev1-(iLay-1)*DZ
nActive=0
do jj=1,NWellLog
do i=1,NLithoBed(jj)-1
    if (z.le.ElevBedBoundary(i,jj).and.z.gt.ElevBedBoundary(i+1,jj)) then
        nActive=nActive+1
    endif    
enddo
enddo
write (999,'(1I10,1F12.2,1I20)') ILay,z,nActive
enddo
write (999,*)

!
GBytes=NLay*NCELL*8/1000000000.0d0
write (999,*) "Layers (NLay)= ", NLay, "; Topo cells (NCELL)= ",NCELL
write (*,*) "Layers (NLay)= ", NLay, "; Topo cells (NCELL)= ",NCELL
write (999,*) "Memory for 4 arrays of entire grid is ",GBytes*3.5d0, " GBytes" 

Allocate(z2(NLay,NCELL),stat=error) !for z2
if (error.ne.0) then 
    write (999,*) "Not enough memory for z2. Reduce NLay or NCELL"
    STOP
endif
Allocate(I2(NLay,NCELL),stat=error) !for I2
if (error.ne.0) then 
    write (999,*) "Not enough memory for I2. Reduce NLay or NCELL"
    STOP
endif

Allocate(z3(NLay,NCELL),stat=error) !for z3
if (error.ne.0) then 
    write (999,*) "Not enough memory for z3. Reduce NLay or NCELL"
    STOP
endif

Allocate(XX(NCELL),stat=error)
Allocate(YY(NCELL),stat=error)
Allocate(VV(NCELL),stat=error)
Allocate(iActiveWell(NCELL),stat=error)
Allocate(VCELL(NCELL),stat=error)
Allocate(indx(NWellLog+1),stat=error)
Allocate(BB(NWellLog+1),stat=error)

write (*,*) "      DEPTH (FT)        NO. SAMPLES         NO. PLANES   TOTAL PLANES"

iz=0
do 20 iLay=1,NLAY
z=StartElev1-(iLay-1)*DZ
iz=iz+1
nActive=0
do 21 jj=1,NWellLog
do 22 i=1,NLithoBed(jj)-1
    if (z.le.ElevBedBoundary(i,jj).and.z.gt.ElevBedBoundary(i+1,jj)) then
        nActive=nActive+1
        iActiveWell(nActive)=jj
        XX(nActive)=WellLogX(jj)
        YY(nActive)=WellLogY(jj)
        VV(nActive)=IndexBedLitho(i,jj)
    endif    
22 enddo
21 enddo

!Check for the same location 
write(999,*)
do ii=1,nActive
do jj=ii+1,nActive
    H=dabs(XX(ii)-XX(jj))+dabs(YY(ii)-YY(jj))
    if (H<0.01d0) then ! if distance is less than 0.01m, move jj point to NW 1m from ii. 
    write(999,*) "Well No. ",iActiveWell(ii)," and ",iActiveWell(jj)," are at the same location."    
    WellLogX(jj)=WellLogX(jj)+10.d0
    WellLogY(jj)=WellLogY(jj)+10.d0
    endif
enddo
enddo
       
if (mod(iz,10).eq.1) write (*,*) z, nActive,iz,"out of ",NLay
if (mod(iz,10).eq.1) write (999,*) z, nActive,iz,"out of ",NLay

    if (nActive.ge.NMinLogs) then  !! The MOST expensive block.
    
       Write (*,*) "Plane No.  :", iz
       Write (999,*) "Plane No.  :", iz
       
       Write (*,*) "Matrix size:", nActive
       Write (999,*) "Matrix size:", nActive
          
       CALL NNINT(nActive,SNGL(XX),SNGL(YY),SNGL(VV),SNGL(XMIN2),SNGL(YMIN2),SNGL(XMAX2),SNGL(YMAX2),NCELL,SNGL(XCELL),SNGL(YCELL), &
                  VCELL) !,NearestWell(i,:))
             
       Call cpu_time (t1)
       do 23 i=1,NCELL  ! This do loop is much more expensive than the above matrix inversion due to large NCELL. "CALL ORDINARY_KRIGING" is fast for a large matrix.                
        dipAngle1=dipAngle(ElevCenter,dipAngleS,dipAngleE,-DistCenter,z) !dip at pivot strike plane
        Dist=A*XCELL(i)+B*YCELL(i)+C
        Displacement=Dist*dipAngle1*3.2808398950131234  !convert meter to foot          
        z1=z-Displacement ! Back tranform to dipped domain
        
        if (z1>TOPO(i)) then
            IK=-1 ! indicator for elevation above land surface
            V_IK=-10.d0
            goto 25  !If elevation is above topo (surface elevation)
        endif 
        
        if (z1<EndElev) then
            IK=-3 ! indicator for elevation below intersted depth
            V_IK=-10.d0
            goto 25  !If elevation is below intersted depth
        endif        
                       
            V_IK= dble(VCELL(i))

            if (V_IK<cutoff) then ! cutoff to determine indicator
                IK=0
            else
                IK=1
            endif                      
25      continue 
        z2(iz,i)=z1 !actual elevation
        I2(iz,i)=IK
        z3(iz,i)=z !non-dip elevation
23     enddo 
       if (Minval(I2(iz,1:NCELL)).eq.-1) then
            write (999,*) "Layer ",iz, ",z= ",nint(z1), " contains air cells."
       endif        
       Call cpu_time (t2)
       Write (*,*) " Interpolation time:", nint(t2-t1), "sec ", nint((t2-t1)/60), "min" 
       Write (999,*) " Interpolation time:", nint(t2-t1), "sec ", nint((t2-t1)/60), "min"  
    else
       write(999,'(F10.1,A18,I4,A40)') sngl(z), "ft(z), less than ", NMinLogs, "well logs for interporation the plane."
       do 26 i=1,NCELL !For entire plane
        dipAngle1=dipAngle(ElevCenter,dipAngleS,dipAngleE,-DistCenter,z) !dip at pivot strike plane
        Dist=A*XCELL(i)+B*YCELL(i)+C
        Displacement=Dist*dipAngle1*3.2808398950131234  !convert meter to foot       
        z1=z-Displacement ! Back tranform to dipped domain
        z2(iz,i)=z1 !actual elevation
        I2(iz,i)=-2 
        z3(iz,i)=z !non-dip elevation
26     enddo 
    endif     
20 enddo

Deallocate(WellLogX,stat=error)
Deallocate(WellLogY,stat=error)
Deallocate(NLithoBed,stat=error)
Deallocate(ElevBedBoundary,stat=error)
Deallocate(IndexBedLitho,stat=error)

Deallocate(XX,stat=error)
Deallocate(YY,stat=error)
Deallocate(VV,stat=error)
Deallocate(iActiveWell,stat=error)
Deallocate(VCELL,stat=error)
Deallocate(indx,stat=error)
Deallocate(BB,stat=error)

write(*,*) "Creating Geo-model.dat and Geo-model_no_thin.dat files."

open (unit=41,file='Geo-model.dat')
write(41,*) "Column ID, X(m), Y(m), Z(ft), Indicator"

Allocate(Ind(NLay),stat=error)

do 60 i=1,NCELL
j_new=1 ! index of new vertical column
ibed=I2(1,i) !the first bed
ind(1)=1 !bed index for the first bed of each column

do 61 j=1,NLay-1 !NLay: the number of bed boundaries, not layers.

if (I2(j,i).ne.ibed) then ! bed change
    j_new=j_new+1
    ibed=I2(j,i) !change lithologic index using current j_th bed
    ind(j_new)=j !index of new bed boundaries for cell i after upscaling
endif
61 enddo

n_new=j_new ! the number of upscaled bed boundaries for cell i, not including the last bed boundary

do 62 j=1,n_new ! document upscaled bed boundaries    
    ibed=I2(ind(j),i)
    write(41,fmt='(1i10,3F15.1,1i10)')  i,XCELL(i),YCELL(i),z2(ind(j),i),ibed
62 enddo
    
write(41,fmt='(1i10,3F15.1,1i10)')  i,XCELL(i),YCELL(i),z2(NLay,i),-999 !end of a vertical column

60 enddo

close (41)

Deallocate(Ind,stat=error)
Deallocate(z2,stat=error) 
Deallocate(I2,stat=error)
Deallocate(z3,stat=error)

!
!
Allocate(z2(NLay,1),stat=error) !for z2
IF (error.ne.0) STOP "*** Not enough memory. Reduce NLay or NCELL ***"
Allocate(I2(NLay,1),stat=error) !for I2

open (unit=41,file='Geo-model.dat')
read (41,*)

open (unit=42,file='Geo-model_no_thin.dat')
write (42,*) "Column ID, X(m), Y(m), Z(ft), Indicator"

max_bed=0
j=0
index1 =-10 !a dummy integer
do i=1,20000000
  read (41,*) ID, X, Y, Z, index
    j=j+1
    z2(j,1)= Z
    I2(j,1)= index
    if (index.eq.-999) then    
        do jj=1,j-1
          if ((z2(jj,1)-z2(jj+1,1)).ge.ThickMin) then 
            if(I2(jj,1).ne.index1) then                 
                write(42,fmt='(1i10,3F15.1,1i10)')  ID,X,Y,z2(jj,1),I2(jj,1)
                index1=I2(jj,1)
            endif
          endif
        enddo
      write(42,fmt='(1i10,3F15.1,1i10)')  ID,X,Y,z2(j,1),I2(j,1)   
      if (j.gt.max_bed) max_bed=j  
      j=0
      index1=-10
    endif          
    if (index.eq.-999.and.ID.eq.NCELL) goto 201

if (i.eq.20000000) then 
    write (*,*) "Geo-model.dat has a problem: The last line has neither index as -999 nor ID as NCELL"
    write (999,*) "Geo-model.dat has a problem: The last line has neither index as -999 nor ID as NCELL"
    stop
endif

enddo
201 continue 
close (41)
close (42)
write (999,*) "Maximum beds for Geo-model.dat: ", max_bed
Deallocate(z2,stat=error) 
Deallocate(I2,stat=error)

Deallocate(XCELL,stat=error) 
Deallocate(YCELL,stat=error)
Deallocate(TOPO,stat=error) 

Call cpu_time (program_end)
write (*,*)
Write (999,*) "Runtime: ", nint((program_end-program_start)/60), "min ", nint((program_end-program_start)/3600), "hour"

write (999,*)
write (999,*) "*****************"
write (999,*) "* GEOMODEL-Tile *" 
write (999,*) "* Program End   *" 
write (999,*) "*****************"
write (999,*)

write(*,*)
write(*,*) "Program completed"

close(999)
stop
end 
!
!*******************************************************************************
! SUBROUTINE: 
!*******************************************************************************
! Read surficial cell location and their elevation
! Input:
! Output: NCELL,XCELL,YCELL,TOPO
Subroutine Grid2D(NCELL, XCELL, YCELL, TOPO)
Implicit Double Precision (A-H,O-Z)
DIMENSION XCELL(NCELL),YCELL(NCELL),TOPO(NCELL)

!MaxC=NCELL

open (unit=16,file='GridTopo.csv')
read(16,*)

i = 0
DO   
   READ(16,*,IOSTAT=io) xx,yy,top
   IF (io > 0) THEN
      WRITE (999,*) 'ERROR: Something is wrong at line ',i
      Stop
   ELSE IF (io < 0) THEN !complete reading the file
      WRITE(999,*)  'The total number of cells is ', i
      if (i.lt.NCELL) then
        write (999,*) "ERROR: Reduce NCELL to ", i
        stop
      endif
      EXIT
   ELSE      
      i = i + 1
      if (i.gt.NCELL) then
        write (999,*) "ERROR: NCELL shoud be more than ", NCELL
        stop
      endif
      XCELL(I) = xx
      YCELL(I) = yy
      TOPO(I)  = top
   END IF
END DO
close(16)

return 
end
!
!
!*******************************************************************************
! SUBROUTINE: 
!*******************************************************************************
! Read well log data
Subroutine WellLogData(XSearchMin,XSearchMax,YSearchMin,YSearchMax,NWellLog,NBed,NLithoBed,WellLogX,WellLogY,ElevBedBoundary,IndexBedLitho) 
Implicit Double Precision (A-H,O-Z)

Character*1000000 line
Character*15 str(NWellLog),str1
Dimension WellLogX(NWellLog),WellLogY(NWellLog),NLithoBed(NWellLog)
Dimension Depth(NBed), ElevBed(NBed),indexBed(NBed)
Dimension ElevBedBoundary(NBed,NWellLog), IndexBedLitho(NBed,NWellLog)

Depth=0.d0
ElevBed=0.d0
indexBed=0

!MaxW=NWellLog

Open (unit=16,file='WellLogs.csv')
read(16,*)

ii=0
Do 100 jj=1,NWellLog
    Depth=0.d0 !initialize depth array
    read (16, '(A)', end=998) line
    read (line,*,end=300) str1,X,Y,DT,(Depth(j),j=1,NBed) !It reads blank values in Depth
    300 continue
    i=0
    do j=1,NBed
        if (j.eq.1) then !alway accept the first depth datum
            i=i+1    
        else
            if (Depth(j).gt.0.d0) i=i+1 !accept dapth datum that is >0
        endif
    enddo
    j=i ! number of beds

! check depth 
    do 10 i=1,j-1
        if (Depth(i+1)-Depth(i).lt.0.0d0) then 
            write (*,*) str1,"Depth error", Depth(i+1)
            write (999,*) str1,"Depth error", Depth(i+1)
            stop
        endif        
10  enddo
    do 11 i=2,j-2
        if (Depth(i+1).eq.Depth(i)) then 
            write (*,*) str1,"Zero thickness", Depth(i+1)
            write (999,*) str1,"Zero thickness", Depth(i+1)
            stop
        endif        
11  enddo

! check litho segment
    if (j.eq.NBed) then 
        write (*,*) str1,"Increase MaxLithoBed over", NBed
        write (999,*) str1,"Increase MaxLithoBed over", NBed
        stop
    endif       
    
!if inside the domain
    if (X.ge.XSearchMin.and.X.le.XSearchMax.and.Y.ge.YSearchMin.and.Y.le.YSearchMax) then
        ii=ii+1
        str(ii)=str1
        WellLogX(ii)=X
        WellLogY(ii)=Y
        write(999,fmt='(A20,1I5,A15,2F25.8)') "Selected well log:", ii,trim(str1),X,Y
        
        Do i=1,j
            ElevBed(i)=DT-Depth(i) 
        enddo        
        Do i=1,j,2
            indexBed(i)=0
            indexBed(i+1)=1         
        enddo
!case 1
        If (Depth(1).ne.Depth(2).and.Depth(j-1).ne.Depth(j)) then 
            ElevBedBoundary(:,ii)=ElevBed(:)
            IndexBedLitho(:,ii)=indexBed(:)
            NLithoBed(ii)=j
            goto 100
        endif
!case 2
        If (Depth(1).eq.Depth(2).and.Depth(j-1).ne.Depth(j)) then 
            Do i=1,j-1
                ElevBedBoundary(i,ii)=ElevBed(i+1)
                IndexBedLitho(i,ii)=indexBed(i+1)
            enddo
            NLithoBed(ii)=j-1
            goto 100
        endif
!case 3
        If (Depth(1).ne.Depth(2).and.Depth(j-1).eq.Depth(j)) then 
            ElevBedBoundary(:,ii)=ElevBed(:)
            IndexBedLitho(:,ii)=indexBed(:)
            NLithoBed(ii)=j-1
            goto 100
        endif
!case 4
        If (Depth(1).eq.Depth(2).and.Depth(j-1).eq.Depth(j)) then 
            Do i=1,j-2
                ElevBedBoundary(i,ii)=ElevBed(i+1)
                IndexBedLitho(i,ii)=indexBed(i+1)
            enddo
            NLithoBed(ii)=j-2
            goto 100
        endif
    endif        
100 enddo
998 Continue
if (ii.ne.NWellLog) then 
    write (999,*) "WARNING: The number of borning wells inside the search domain is ", ii
endif
NWellLog=ii
Close(16)

return 
end
!
!
!
!*******************************************************************************
! FUNCTION: 
!*******************************************************************************
FUNCTION dipAngle(ElevCenter,dipAngleS,dipAngleE,Dist,z) !dip angle at pivot strike
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
real*8 dipAngle,E1,E2,dipStart,dipEnd,z

if (dipAngleS.eq.dipAngleE) then 
    dipAngle=dipAngleS
else
    E1=ElevCenter-Dist*dipAngleS*3.2808398950131234
    E2=ElevCenter-Dist*dipAngleE*3.2808398950131234

    if (z.ge.E1) dipAngle=dipAngleS
    if (z.le.E2) dipAngle=dipAngleE ! if E1 and E2 at the same depth (step function), assign dipEnd to E2 and below 
    if (z.lt.E1.and.z.gt.E2) then
        dipAngle=(ElevCenter-z)/(Dist*3.2808398950131234)
    endif
endif  

RETURN
END
!