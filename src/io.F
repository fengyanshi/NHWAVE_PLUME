! ----------------------------------------------------
!    This is subroutine to read input.txt
!  called by 
!        Master program
!   this input file is for both AMR and single grid version 
!   last update 05/24/2013, fyshi
! --------------------------------------------------

SUBROUTINE READ_INPUT
     USE GLOBAL
# if defined (AMR)
     USE NESTING
# endif
     USE Input_Util
# if defined (INTEL)
     use IFPORT, only: rand
# endif
! fyshi add foam
# if defined (FOAM)
     USE MOD_FOAM
# endif
     IMPLICIT NONE
    CHARACTER(LEN=2) Ngrid
    CHARACTER(LEN=80) SubGrid_FILE
    INTEGER :: itmp1,itmp2,itmp3,itmp4,itmp5
     character(len=80) :: FILE_NAME,TMP_NAME
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,  &
                 Fk,Fkdif,Theta_Calc,Wnumy,tmp,tmp1,DFreq, &
                 Freq_Peak,gam,sa,sb,SumInt,A_Jon
     integer :: line,ierr,Iter,i,j,k
 
     ! log and error file
     open(3,file='log.txt')

     ! read from input.txt
     FILE_NAME='input.txt'

     ! title
     CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# else
       write(3,*) 'No TITLE in ', FILE_NAME, 'use default'
# endif
       TITLE='---TEST RUN---'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'---- LOG FILE ---'
     if(myid.eq.0) WRITE(3,*)TITLE
     if(myid.eq.0) WRITE(3,*)'--------------input start --------------'
# else
     WRITE(3,*)'---- LOG FILE ---'
     WRITE(3,*)TITLE
     WRITE(3,*)'--------------input start --------------'
# endif

     ! dimension                                             
     CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
     CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
     CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# else
     WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     WRITE(3,'(A7,I5)')'Kglob= ',Kglob
# endif

# if defined (PARALLEL)
     ! processor number
     CALL GET_INTEGER_VAL(PX,FILE_NAME,'PX',line)
     CALL GET_INTEGER_VAL(PY,FILE_NAME,'PY',line)
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PX= ',PX
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PY= ',PY
     if(PX*PY.ne.NumP) then
       if(myid.eq.0) WRITE(3,'(A6,I5)') 'NumP= ',NumP
       stop
     endif
# endif

# if defined (AMR)
! ----------------NESTING GRID INFO

! number of grid
      CALL GET_INTEGER_VAL(NumGrid,FILE_NAME,'SubGrid',line)
      NumGrid = 1+NumGrid
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A10,I12)')'NumGrid = ', NumGrid
# else
      WRITE(3,'(A10,I12)')'NumGrid = ', NumGrid
# endif
      ALLOCATE(GridDimX(NumGrid),GridDimY(NumGrid),RATIO_SPACING(NumGrid),&
           MboxRef(NumGrid),NboxRef(NumGrid),IsPrint(NumGrid), &
           MboxRef0(NumGrid),NboxRef0(NumGrid), &
           MboxRefFix(NumGrid),NboxRefFix(NumGrid))
      DO K=1,NumGrid
        GridDimX(K)=Mglob
        GridDimY(K)=Nglob
        RATIO_SPACING(K)=1
        MboxRef(K)=Mglob+3*Nghost
        NboxRef(K)=Nglob+3*Nghost
        MboxRef0(K)=Mglob+3*Nghost ! make sure out of domain first step
        NboxRef0(K)=Nglob+3*Nghost
        IsPrint(K)=.FALSE.
      ENDDO

      IF(NumGrid.GT.1)THEN

        CALL GET_STRING_VAL(SubGrid_FILE,FILE_NAME,'SubGrid_FILE',line,ierr)
        OPEN(2,FILE=TRIM(SubGrid_FILE))
          READ(2,*) ! title
          DO K=1,NumGrid-1
           READ(2,*)itmp1,itmp2,itmp3,itmp4,itmp5
           IF(itmp3.LT.1)THEN
# if defined (PARALLEL)
            IF(myid.EQ.0)THEN
             WRITE(*,*)'STOP! Spacing ratio should larger than 1.'
            ENDIF
             call MPI_FINALIZE ( ier )          
# else
            WRITE(*,*)'STOP! Spacing ratio should larger than 1.'
            STOP
# endif
           ENDIF
           GridDimX(K+1)=itmp1
           GridDimY(K+1)=itmp2
           RATIO_SPACING(K+1)=itmp3
           MboxRefFix(K+1)=itmp4
           NboxRefFix(K+1)=itmp5
          ENDDO
        CLOSE(2)

        CALL GET_STRING_VAL(DEPTH_FULL_FILE,FILE_NAME,'DEPTH_FULL_FILE',line,ierr)    
        CALL GET_INTEGER_VAL(Mdep,FILE_NAME,'Mdep',line) 
        CALL GET_INTEGER_VAL(Ndep,FILE_NAME,'Ndep',line)
        CALL GET_INTEGER_VAL(MainGridRatio,FILE_NAME,'MainGridRatio',line)
        ALLOCATE(Depth_Full(Mdep,Ndep)) 

      ENDIF ! end numgrid>1

      ALLOCATE(IS_GRID_MOVE(NumGrid))
      ALLOCATE(  &
        ng_Kstart_EAST(NumGrid),ng_Kend_EAST(NumGrid),  &
                  ng_Kshift_EAST(NumGrid),&
                  ng_Kstart_WEST(NumGrid),ng_Kend_WEST(NumGrid),  &
                  ng_Kshift_WEST(NumGrid),&
                  ng_Kstart_SOUTH(NumGrid),ng_Kend_SOUTH(NumGrid),  &
                  ng_Kshift_SOUTH(NumGrid),&
                  ng_Kstart_NORTH(NumGrid),ng_Kend_NORTH(NumGrid),  &
                  ng_Kshift_NORTH(NumGrid), &
                  ng_IN_DOMAIN(NumGrid), &
                  START_GRID(NumGrid), &
                  Istart_overlap(NumGrid), &
                  Iend_overlap(NumGrid), & 
                  Jstart_overlap(NumGrid), &
                  Jend_overlap(NumGrid), &
                  Ishift_overlap(NumGrid),&
                  Jshift_overlap(NumGrid), &
                  I_IN_DOMAIN_overlap(NumGrid), &
                  J_IN_DOMAIN_overlap(NumGrid) &
                   )

      DO K=1,NumGrid
        IS_GRID_MOVE(K)=.TRUE.
        START_GRID(K)=.FALSE.
      ENDDO
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'NESTING_X1D= ',NESTING_X1D
# else
     WRITE(3,'(A14,L4)')'NESTING_X1D= ',NESTING_X1D
# endif


     CALL GET_LOGICAL_VAL(NESTING_X1D,FILE_NAME,'NESTING_X1D',line,ierr)

! ------------  END NESTING GRID INFO
# endif
  ! end AMR

     ! grid sizes
     CALL GET_Float_VAL(dx,FILE_NAME,'DX',line)
     CALL GET_Float_VAL(dy,FILE_NAME,'DY',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DX= ',dx
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DY= ',dy
# else
     WRITE(3,'(A4,F8.4)')'DX= ',dx
     WRITE(3,'(A4,F8.4)')'DY= ',dy
# endif

     ! vertical grid option
     call GET_INTEGER_VAL(Ivgrd,FILE_NAME,'IVGRD',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
# else
     WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
# endif
     IF(Ivgrd==2)THEN
       CALL GET_Float_VAL(Grd_R,FILE_NAME,'GRD_R',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R
# else
     WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R
# endif
     ENDIF

     ! time step
     CALL GET_Float_VAL(dt_ini,FILE_NAME,'DT_INI',line)
     CALL GET_Float_VAL(dt_min,FILE_NAME,'DT_MIN',line)
     CALL GET_Float_VAL(dt_max,FILE_NAME,'DT_MAX',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max
# else
     WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max
# endif

     ! result folder                                     
     CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# else
     WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER
# endif

     ! simulation steps and time
     call GET_INTEGER_VAL(SIM_STEPS,FILE_NAME,'SIM_STEPS',line)
     CALL GET_Float_VAL(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',line)
     CALL GET_Float_VAL(Plot_Start,FILE_NAME,'PLOT_START',line)
     CALL GET_Float_VAL(Plot_Intv,FILE_NAME,'PLOT_INTV',line)
     CALL GET_Float_VAL(Screen_Intv,FILE_NAME,'SCREEN_INTV',line)
     CALL GET_LOGICAL_VAL(HOTSTART,FILE_NAME,'HOTSTART',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     if(myid.eq.0) WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
     if(myid.eq.0) WRITE(3,'(A11,L4)')'HOTSTART= ',HOTSTART
# else
     WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
     WRITE(3,'(A11,L4)')'HOTSTART= ',HOTSTART
# endif

     ! laminar viscosity
     CALL GET_LOGICAL_VAL(VISCOUS_FLOW,FILE_NAME,'VISCOUS_FLOW',line,ierr)
     CALL GET_INTEGER_VAL(IVturb,FILE_NAME,'IVTURB',line)
     CALL GET_INTEGER_VAL(IHturb,FILE_NAME,'IHTURB',line)
     CALL GET_Float_VAL(Visc,FILE_NAME,'VISCOSITY',line)
     CALL GET_Float_VAL(Schmidt,FILE_NAME,'Schmidt',line)
     CALL GET_Float_VAL(Cvs,FILE_NAME,'Cvs',line)
     CALL GET_Float_VAL(Chs,FILE_NAME,'Chs',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IVTURB= ',IVturb
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IHTURB= ',IHturb
     if(myid.eq.0) WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     if(myid.eq.0) WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Chs= ',Chs
# else
     WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     WRITE(3,'(A7,I2)')'IVTURB= ',IVturb
     WRITE(3,'(A7,I2)')'IHTURB= ',IHturb
     WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     WRITE(3,'(A6,F8.3)')'Chs= ',Chs
# endif


     ! viscous number
     CALL GET_Float_VAL(VISCOUS_NUMBER,FILE_NAME,'VISCOUS_NUMBER',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# else
     WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER
# endif

     IF(HOTSTART)THEN
       CALL GET_STRING_VAL(Eta_HotStart_File,FILE_NAME,'Eta_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(U_HotStart_File,FILE_NAME,'U_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(V_HotStart_File,FILE_NAME,'V_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(W_HotStart_File,FILE_NAME,'W_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(P_HotStart_File,FILE_NAME,'P_HotStart_File',line,ierr)
# if defined (SALINITY)
       CALL GET_STRING_VAL(Sali_HotStart_File,FILE_NAME,'Sali_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(Temp_HotStart_File,FILE_NAME,'Temp_HotStart_File',line,ierr)
# endif
      IF(VISCOUS_FLOW)THEN
       CALL GET_STRING_VAL(Rho_HotStart_File,FILE_NAME,'Rho_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(TKE_HotStart_File,FILE_NAME,'TKE_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(EPS_HotStart_File,FILE_NAME,'EPS_HotStart_File',line,ierr)
      ENDIF
# if defined (AIR_PRESSURE)
       CALL GET_STRING_VAL(Pmb_HotStart_File,FILE_NAME,'Pmb_HotStart_File',line,ierr)
# endif

     ENDIF ! end hotstart

# if defined (FROUDE_CAP)
     CALL GET_Float_VAL(FROUDECAP,FILE_NAME,'FROUDE_CAP',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A13,F8.3)')'FROUDE_CAP= ',FROUDECAP
# else
     WRITE(3,'(A13,F8.3)')'FROUDE_CAP= ',FROUDECAP
# endif

# endif
! end froude cap

! fyshi added bathymetry file name 04/13/2012

     CALL GET_STRING_VAL(Depth_File,FILE_NAME,'DEPTH_FILE',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'Depth File= ', TRIM(Depth_File)
# else
       WRITE(3,'(A14,A20)')'Depth File= ', TRIM(Depth_File)
# endif

     ! courant number
     CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# else
     WRITE(3,'(A5,F8.3)')'CFL= ',CFL
# endif

     ! minimum depth
     CALL GET_Float_VAL(MinDep,FILE_NAME,'MinDep',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# else
     WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep
# endif

     ! bathymetry     
     CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)
     CALL GET_LOGICAL_VAL(ANA_BATHY,FILE_NAME,'ANA_BATHY',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     if(myid.eq.0) WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY
# else
     WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY
# endif

     ! bottom roughness
     CALL GET_INTEGER_VAL(Ibot,FILE_NAME,'Ibot',line)
     CALL GET_Float_VAL(Cd0,FILE_NAME,'Cd0',line)
     CALL GET_Float_VAL(Zob,FILE_NAME,'Zob',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# else
     WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     WRITE(3,'(A5,F8.5)')'Zob= ',Zob
# endif

     ! barotropic or baroclinic
     CALL GET_LOGICAL_VAL(BAROTROPIC,FILE_NAME,'BAROTROPIC',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC
# else
     WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC
# endif

     ! numerical scheme
     CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
     CALL GET_STRING_VAL(TIME_ORDER,FILE_NAME,'TIME_ORDER',line,ierr)
     CALL GET_STRING_VAL(CONVECTION,FILE_NAME,'CONVECTION',line,ierr)
     CALL GET_LOGICAL_VAL(ADV_HLLC,FILE_NAME,'HLLC',line,ierr)
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       HIGH_ORDER='SECOND'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# else
     WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
# endif
     IF(ierr==1)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# else
       WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
# endif
       TIME_ORDER='THIRD'
     ENDIF
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# else
     WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER
# endif

     ! ramp up the simulation
     CALL GET_Float_VAL(TRamp,FILE_NAME,'TRAMP',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp
# else
     WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp
# endif

     ! if non-hydrostatic simulation
     CALL GET_LOGICAL_VAL(NON_HYDRO,FILE_NAME,'NON_HYDRO',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# else
     WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO
# endif

     ! poisson solver
     CALL GET_INTEGER_VAL(isolver,FILE_NAME,'ISOLVER',line)
     CALL GET_INTEGER_VAL(itmax,FILE_NAME,'ITMAX',line)
     CALL GET_Float_VAL(tol,FILE_NAME,'TOL',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     if(myid.eq.0) WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     if(myid.eq.0) WRITE(3,'(A5,E12.3)')'TOL= ',tol
# else
     WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     WRITE(3,'(A5,E12.3)')'TOL= ',tol
# endif

     ! periodic bc
     CALL GET_LOGICAL_VAL(PERIODIC_X,FILE_NAME,'PERIODIC_X',line,ierr)
     CALL GET_LOGICAL_VAL(PERIODIC_Y,FILE_NAME,'PERIODIC_Y',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y
# else
     WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y
# endif

     ! boundary type
     CALL GET_INTEGER_VAL(Bc_X0,FILE_NAME,'BC_X0',line)
     CALL GET_INTEGER_VAL(Bc_Xn,FILE_NAME,'BC_Xn',line)
     CALL GET_INTEGER_VAL(Bc_Y0,FILE_NAME,'BC_Y0',line)
     CALL GET_INTEGER_VAL(Bc_Yn,FILE_NAME,'BC_Yn',line)
     CALL GET_INTEGER_VAL(Bc_Z0,FILE_NAME,'BC_Z0',line)
     CALL GET_INTEGER_VAL(Bc_Zn,FILE_NAME,'BC_Zn',line)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# else
     WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn
# endif

     ! wavemaker 
     CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# else
     WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
# endif

! absorbing generating wavemaker
     IF(WaveMaker(1:3)=='ABS') THEN
       CALL GET_STRING_VAL(Wave_Comp_File,FILE_NAME,'LinearWaveSerFile',line,ierr)
       OPEN(1,FILE=TRIM(Wave_Comp_File))
         READ(1,*) ! title line
         READ(1,*) Dep_Ser
         READ(1,*) U_FLOW_LEFT,U_FLOW_RIGHT
         READ(1,*) NumCompSer
         ALLOCATE(Amp_Ser(NumCompSer),Per_Ser(NumCompSer),&
                  Phase_Ser(NumCompSer),Theta_Ser(NumCompSer),&
                  Segma_Ser(NumCompSer),Wave_Number_Ser(NumCompSer),&
                  Stokes_Drift_Ser(NumCompSer))
         DO I=1,NumCompSer
           READ(1,*)Amp_Ser(I),Per_Ser(I),Phase_Ser(I),Theta_Ser(I)
           IF(Per_Ser(I)==ZERO)THEN
             WRITE(*,*)'input wave frequency is zero, stop'
             STOP
           ELSE
             Per_Ser(I)=1.0_SP/Per_Ser(I)
           ENDIF
         ENDDO
       CLOSE(1)
     
       CALL GET_Float_VAL(WaveMaker_West_Width,FILE_NAME,'WaveMaker_West_Width',line)
       CALL GET_Float_VAL(WaveMaker_East_Width,FILE_NAME,'WaveMaker_East_Width',line)
       CALL GET_Float_VAL(R_WaveMaker,FILE_NAME,'R_WaveMaker',line)
       CALL GET_Float_VAL(A_WaveMaker,FILE_NAME,'A_WaveMaker',line)

   ENDIF  ! end abs wavemaker


     IF(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='INT'.or. &
       WaveMaker(1:3)=='FLU')THEN
       CALL GET_Float_VAL(Amp_Wave,FILE_NAME,'AMP',line)
       CALL GET_Float_VAL(Per_Wave,FILE_NAME,'PER',line)
       CALL GET_Float_VAL(Dep_Wave,FILE_NAME,'DEP',line)
       CALL GET_Float_VAL(Theta_Wave,FILE_NAME,'THETA',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       if(myid.eq.0) WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
# else
       WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
# endif

       IF(WaveMaker(1:3)=='INT') then
         CALL GET_Float_VAL(Xsource_West,FILE_NAME,'Xsource_West',line)
         CALL GET_Float_VAL(Xsource_East,FILE_NAME,'Xsource_East',line)
         CALL GET_Float_VAL(Ysource_Suth,FILE_NAME,'Ysource_Suth',line)
         CALL GET_Float_VAL(Ysource_Nrth,FILE_NAME,'Ysource_Nrth',line)
# if defined (PARALLEL)
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# else
         WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West 
         WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
# endif
       ENDIF

       ! test periodicity
       IF(PERIODIC_Y.and.Theta_Wave.ne.Zero) then
         ! find wave number
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length

         Iter = 0
 75      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85      continue

         if(Theta_Wave>Zero) then
           ! find right angle for periodic bc
           tmp = Large       
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         elseif(Theta_Wave<Zero) then
           ! find right angle for periodic bc 
           tmp = Large
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = -asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         endif

# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave         
           write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
         endif
# else
         write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave                                                      
         write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
# endif
         Theta_Wave = Theta_Calc
       ENDIF
     ENDIF

     ! random wave, read in 2d spectrum
     IF(WaveMaker(5:7)=='SPC') then
       open(14,file='spc2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
         stop
# endif
       endif
       if(NumDir>MaxNumDir) then
# if defined (PARALLEL)
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
# else
         write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
         stop
# endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
       close(14)

       ! random phase for each component
       do j = 1,NumFreq
       do i = 1,NumDir
# if defined (INTEL)
         Random_Phs(i,j) = rand()*2.0*pi
# else
         Random_Phs(i,j) = rand(0)*2.0*pi
# endif
       enddo
       enddo
     ENDIF

     ! JONSWAP spectrum
     if(WaveMaker(5:7)=='JON') then
       CALL GET_Float_VAL(Hm0,FILE_NAME,'Hm0',line)
       CALL GET_Float_VAL(Tp,FILE_NAME,'Tp',line)
       CALL GET_Float_VAL(Freq_Min,FILE_NAME,'Freq_Min',line)
       CALL GET_Float_VAL(Freq_Max,FILE_NAME,'Freq_Max',line)
       CALL GET_INTEGER_VAL(NumFreq,FILE_NAME,'NumFreq',line) 

       ! frequency
       DFreq = (Freq_Max-Freq_Min)/NumFreq
       do i = 1,NumFreq
         Freq(i) = Freq_Min+0.5*DFreq+(i-1)*DFreq
       enddo

       ! jonswap spectrum
       gam = 3.3; sa = 0.07; sb = 0.09
       Freq_Peak = 1.0/Tp
       do i = 1,NumFreq
         if(Freq(i)<Freq_Peak) then
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sa**2)
         else
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sb**2)
         endif
       enddo
         
       ! make sure m0=Hm0**2/16=int S(f)df
       SumInt = Zero
       do i = 1,NumFreq
         SumInt = SumInt+Jon_Spc(i)*DFreq
       enddo
       A_Jon = Hm0**2/16.0/SumInt

       do i = 1,NumFreq
         Jon_Spc(i) = Jon_Spc(i)*A_Jon
         RanPhs(i) = rand()*2.0*pi
       enddo
     endif

     ! sponge layer
     CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# else
     WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
# endif
     IF(SPONGE_ON)THEN
       CALL GET_Float_VAL(Sponge_West_Width,FILE_NAME,'Sponge_West_Width',line)
       CALL GET_Float_VAL(Sponge_East_Width,FILE_NAME,'Sponge_East_Width',line)
       CALL GET_Float_VAL(Sponge_South_Width,FILE_NAME,'Sponge_South_Width',line)
       CALL GET_Float_VAL(Sponge_North_Width,FILE_NAME,'Sponge_North_Width',line)
       CALL GET_Float_VAL(R_Sponge,FILE_NAME,'R_Sponge',line)
       CALL GET_Float_VAL(A_Sponge,FILE_NAME,'A_Sponge',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
       if(myid.eq.0) WRITE(3,'(A10,F6.3)')'R_Sponge= ', R_Sponge
       if(myid.eq.0) WRITE(3,'(A10,F6.3)')'A_Sponge= ', A_Sponge
# else
       WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
       WRITE(3,'(A10,F6.3)')'R_Sponge= ', R_Sponge
       WRITE(3,'(A10,F6.3)')'A_Sponge= ', A_Sponge
# endif
     ENDIF

     ! wave average control
     CALL GET_LOGICAL_VAL(WAVE_AVERAGE_ON,FILE_NAME,'WAVE_AVERAGE_ON',line,ierr)
     CALL GET_Float_VAL(Wave_Ave_Start,FILE_NAME,'WAVE_AVERAGE_START',line)
     CALL GET_Float_VAL(Wave_Ave_End,FILE_NAME,'WAVE_AVERAGE_END',line)
     CALL GET_INTEGER_VAL(WaveheightID,FILE_NAME,'WaveheightID',line)

# if defined (SEDIMENT)
     ! sediment parameters
     CALL GET_STRING_VAL(Sed_Type,FILE_NAME,'Sed_Type',line,ierr)
     CALL GET_STRING_VAL(Sed_Load,FILE_NAME,'Sed_Load',line,ierr)
     CALL GET_LOGICAL_VAL(COUPLE_FS,FILE_NAME,'COUPLE_FS',line) 
     CALL GET_Float_VAL(SD50,FILE_NAME,'D50',line)
     CALL GET_Float_VAL(Shields_c,FILE_NAME,'Shields_c',line)
     CALL GET_Float_VAL(Af,FILE_NAME,'Af',line)
     CALL GET_Float_VAL(Tau_ce,FILE_NAME,'Tau_ce',line)
     CALL GET_Float_VAL(Tau_cd,FILE_NAME,'Tau_cd',line)
     CALL GET_Float_VAL(Erate,FILE_NAME,'Erate',line)
     CALL GET_Float_VAL(Mud_Visc,FILE_NAME,'Mud_Visc',line)
     CALL GET_Float_VAL(Tim_Sedi,FILE_NAME,'Tim_Sedi',line)
     CALL GET_LOGICAL_VAL(BED_CHANGE,FILE_NAME,'BED_CHANGE',line)
# endif

# if defined (LANDSLIDE)
     ! landslide parameters                           
     CALL GET_STRING_VAL(SlideType,FILE_NAME,'SlideType',line,ierr)
     CALL GET_Float_VAL(SlideT,FILE_NAME,'SlideT',line)
     CALL GET_Float_VAL(SlideL,FILE_NAME,'SlideL',line)
     CALL GET_Float_VAL(SlideW,FILE_NAME,'SlideW',line)
     CALL GET_Float_VAL(SlideD,FILE_NAME,'SlideD',line)
     CALL GET_Float_VAL(SlideAngle,FILE_NAME,'SlideAngle',line)
     CALL GET_Float_VAL(SlideX0,FILE_NAME,'SlideX0',line)
     CALL GET_Float_VAL(SlideY0,FILE_NAME,'SlideY0',line)
     if(trim(SlideType)=='RIGID') then
       CALL GET_Float_VAL(SlideUt,FILE_NAME,'SlideUt',line)
       CALL GET_Float_VAL(SlideA0,FILE_NAME,'SlideA0',line)
     endif
     if(trim(SlideType)=='DEFORMABLE') then
       CALL GET_Float_VAL(SlideDens,FILE_NAME,'SlideDens',line)
     endif
     
     ! whether to consider rheology
     CALL GET_LOGICAL_VAL(RHEOLOGY_ON,FILE_NAME,'RHEOLOGY_ON',line)
     CALL GET_Float_VAL(Yield_Stress,FILE_NAME,'Yield_Stress',line)
     CALL GET_Float_VAL(Plastic_Visc,FILE_NAME,'Plastic_Visc',line)
# endif

# if defined (VEGETATION)
     ! vegetation parameters
     CALL GET_Float_VAL(Veg_X0,FILE_NAME,'Veg_X0',line)
     CALL GET_Float_VAL(Veg_Xn,FILE_NAME,'Veg_Xn',line)
     CALL GET_Float_VAL(Veg_Y0,FILE_NAME,'Veg_Y0',line)
     CALL GET_Float_VAL(Veg_Yn,FILE_NAME,'Veg_Yn',line)
     CALL GET_Float_VAL(VegH,FILE_NAME,'VegH',line)
     CALL GET_Float_VAL(VegDens,FILE_NAME,'VegDens',line)
     CALL GET_Float_VAL(VegDrag,FILE_NAME,'VegDrag',line)
# endif

# if defined (OBSTACLE)
     CALL GET_STRING_VAL(MASK3D_File,FILE_NAME,'MASK3D_FILE',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'MASK3D File= ', TRIM(MASK3D_File)
# else
       WRITE(3,'(A14,A20)')'MASK3D File= ', TRIM(MASK3D_File)
# endif
# endif

     ! if there is external forcing
     CALL GET_LOGICAL_VAL(EXTERNAL_FORCING,FILE_NAME,'EXTERNAL_FORCING',line,ierr)

     ! probe output
     CALL GET_INTEGER_VAL(NSTAT,FILE_NAME,'NSTAT',line)
     CALL GET_Float_VAL(Plot_Intv_Stat,FILE_NAME,'PLOT_INTV_STAT',line)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# else
       WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat
# endif

     if(NSTAT>0) then
       open(15,file='stat.txt',status='old')
       do i = 1,NSTAT
         read(15,*) xstat(i),ystat(i)
       enddo
       close(15)
     endif

     ! output parameters  
     CALL GET_LOGICAL_VAL(OUT_DEP,FILE_NAME,'OUT_DEP',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_ETA,FILE_NAME,'OUT_ETA',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'OUT_U',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'OUT_V',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_W,FILE_NAME,'OUT_W',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_P,FILE_NAME,'OUT_P',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_TKE,FILE_NAME,'OUT_TKE',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_EPS,FILE_NAME,'OUT_EPS',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_S,FILE_NAME,'OUT_S',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_MU,FILE_NAME,'OUT_MU',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_BUB,FILE_NAME,'OUT_BUB',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_A,FILE_NAME,'OUT_A',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_F,FILE_NAME,'OUT_F',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_T,FILE_NAME,'OUT_T',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_G,FILE_NAME,'OUT_G',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_SALI,FILE_NAME,'OUT_SALI',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_TEMP,FILE_NAME,'OUT_TEMP',line,ierr)
     CALL GET_LOGICAL_VAL(OUT_RHO,FILE_NAME,'OUT_RHO',line,ierr)
! fyshi add foam
# if defined (FOAM)
     CALL GET_LOGICAL_VAL(OUT_FOAM,FILE_NAME,'OUT_FOAM',line,ierr)
# endif

! fyshi added nesting option 05/15/2013
# if defined (COUPLING)
      CALL GET_STRING_VAL(COUPLING_FILE,FILE_NAME,'COUPLING_FILE',line,ierr)
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A15,A50)')'COUPLING_FILE:', COUPLING_FILE
# else
      WRITE(3,'(A15,A50)')'COUPLING_FILE:', COUPLING_FILE
# endif
# endif

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_DEP= ',OUT_DEP
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_ETA= ',OUT_ETA
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_U= ',OUT_U
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_V= ',OUT_V
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_W= ',OUT_W
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_P= ',OUT_P
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_TKE= ',OUT_TKE
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_EPS= ',OUT_EPS
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_S= ',OUT_S
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_MU= ',OUT_MU
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_BUB= ',OUT_BUB
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_A= ',OUT_A
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_F= ',OUT_F
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_T= ',OUT_T
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_G= ',OUT_G
     if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_SALI= ',OUT_SALI
    if(myid.eq.0) WRITE(3,'(A10,L4)') 'OUT_TEMP= ',OUT_TEMP
# else
     WRITE(3,'(A10,L4)') 'OUT_DEP= ',OUT_DEP
     WRITE(3,'(A10,L4)') 'OUT_ETA= ',OUT_ETA
     WRITE(3,'(A10,L4)') 'OUT_U= ',OUT_U
     WRITE(3,'(A10,L4)') 'OUT_V= ',OUT_V
     WRITE(3,'(A10,L4)') 'OUT_W= ',OUT_W
     WRITE(3,'(A10,L4)') 'OUT_P= ',OUT_P
     WRITE(3,'(A10,L4)') 'OUT_TKE= ',OUT_TKE
     WRITE(3,'(A10,L4)') 'OUT_EPS= ',OUT_EPS
     WRITE(3,'(A10,L4)') 'OUT_S= ',OUT_S
     WRITE(3,'(A10,L4)') 'OUT_MU= ',OUT_MU
     WRITE(3,'(A10,L4)') 'OUT_BUB= ',OUT_BUB
     WRITE(3,'(A10,L4)') 'OUT_A= ',OUT_A
     WRITE(3,'(A10,L4)') 'OUT_F= ',OUT_F
     WRITE(3,'(A10,L4)') 'OUT_T= ',OUT_T
     WRITE(3,'(A10,L4)') 'OUT_G= ',OUT_G
     WRITE(3,'(A10,L4)') 'OUT_SALI= ',OUT_SALI
     WRITE(3,'(A10,L4)') 'OUT_TEMP= ',OUT_TEMP
# endif

# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,*)'--------------input end --------------'
# else
     WRITE(3,*)'--------------input end --------------'
# endif

! fyshi separate boundary condition from wavemaker in order to do both wavemaker
! and boundary conditions
     CALL GET_STRING_VAL(BOUNDARY,FILE_NAME,'BOUNDARY',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'BOUNDARY= ', BOUNDARY
# else
     WRITE(3,'(A11,A50)')'BOUNDARY= ', BOUNDARY
# endif

! fyshi add some must-parameters 12/16/2011
     IF(BOUNDARY(1:10)=='TID_FLX_LR')THEN
       BC_X0 = 3
       BC_Xn = 3
     ELSEIF(BOUNDARY(1:10)=='TID_ELE_LR')THEN
       BC_X0 = 1
       BC_Xn = 1
     ENDIF
     IF(BOUNDARY(1:10)=='TID_FLX_LR'.OR.  &
        BOUNDARY(1:10)=='TID_ELE_LR') THEN

       CALL GET_STRING_VAL(BoundaryFile,FILE_NAME,'BoundaryFile',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A14,A20)')'BoundaryFile= ', TRIM(BoundaryFile)
# else
       WRITE(3,'(A14,A20)')'BoundaryFile= ', TRIM(BoundaryFile)
# endif

       OPEN(2,FILE=TRIM(BoundaryFile))
         READ(2,*)WHAT
         READ(2,*)NumTimeData
         ALLOCATE (TimeData(NumTimeData),  &
                   DataU_L(NumTimeData), DataU_R(NumTimeData), DataSal_R(NumTimeData),&   
                   DataEta_L(NumTimeData), DataEta_R(NumTimeData),DataSal_L(NumTimeData),&
                   DataTem_L(NumTimeData),DataTem_R(NumTimeData) )   
         ALLOCATE (Z_pct_West(Kglob),Z_pct_East(Kglob))
         READ(2,*)WHAT ! z percentage of flux
         READ(2,*)WHAT ! west side
         READ(2,*)(Z_pct_West(I),I=1,Kglob)
         READ(2,*)WHAT ! east side
         READ(2,*)(Z_pct_East(I),I=1,Kglob)   
         READ(2,*)WHAT ! time series
         DO I=1,NumTimeData
           READ(2,*,END=111) TimeData(I)
           READ(2,*)WHAT ! west side    
           READ(2,*,END=111) DataEta_L(I),DataU_L(I),DataSal_L(I), &
                                         DataTem_L(I)
           READ(2,*)WHAT ! east side 
           READ(2,*,END=111) DataEta_R(I),DataU_R(I),DataSal_R(I), &
                                         DataTem_R(I)
         ENDDO
111      CONTINUE
# if defined (PARALLEL)
        IF(myid==0) WRITE(*,*) 'Number of Data: ', I-1
# else
        WRITE(*,*) 'Number of Data: ', I-1
# endif
       CLOSE(2)
     ENDIF ! end flx boundary conditions

! fyshi add initial sali and temp data 11/03/2012

# if defined(SALINITY)
     CALL GET_STRING_VAL(INI_SALI_INPUT,FILE_NAME,'INI_SALI_INPUT',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A17,A50)')'INI_SALI_INPUT= ', INI_SALI_INPUT
# else
     WRITE(3,'(A17,A50)')'INI_SALI_INPUT= ', INI_SALI_INPUT
# endif

     IF(INI_SALI_INPUT(1:4)=='CONS')THEN
        CALL GET_Float_VAL(INI_SALI,FILE_NAME,'INI_SALI',line)
# if defined (PARALLEL)
        if(myid.eq.0) WRITE(3,'(A11,F12.6)')'INI_SALI= ', INI_SALI
# else
         WRITE(3,'(A11,F12.6)')'INI_SALI= ', INI_SALI
# endif
     ELSE IF(INI_SALI_INPUT(1:4)=='DATA')THEN
       CALL GET_STRING_VAL(INI_SALI_FILE,FILE_NAME,'INI_SALI_FILE',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A16,A20)')'INI_SALI_FILE= ', TRIM(INI_SALI_FILE)
# else
       WRITE(3,'(A16,A20)')'INI_SALI_FILE= ', TRIM(INI_SALI_FILE)
# endif
     ENDIF ! end read sali
# endif
!  end salinity

# if defined(TEMPERATURE)
     CALL GET_STRING_VAL(INI_TEMP_INPUT,FILE_NAME,'INI_TEMP_INPUT',line,ierr)
# if defined (PARALLEL)
     if(myid.eq.0) WRITE(3,'(A17,A50)')'INI_TEMP_INPUT= ', INI_TEMP_INPUT
# else
     WRITE(3,'(A17,A50)')'INI_TEMP_INPUT= ', INI_TEMP_INPUT
# endif

     IF(INI_TEMP_INPUT(1:4)=='CONS')THEN
        CALL GET_Float_VAL(INI_TEMP,FILE_NAME,'INI_TEMP',line)
# if defined (PARALLEL)
        if(myid.eq.0) WRITE(3,'(A11,F12.6)')'INI_TEMP= ', INI_TEMP
# else
         WRITE(3,'(A11,F12.6)')'INI_TEMP= ', INI_TEMP
# endif
     ELSE IF(INI_TEMP_INPUT(1:4)=='DATA')THEN
       CALL GET_STRING_VAL(INI_TEMP_FILE,FILE_NAME,'INI_TEMP_FILE',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A16,A20)')'INI_TEMP_FILE= ', TRIM(INI_TEMP_FILE)
# else
       WRITE(3,'(A16,A20)')'INI_TEMP_FILE= ', TRIM(INI_TEMP_FILE)
# endif
     ENDIF ! end read temperatre
# endif
!  end temperature

! fyshi add tidal low pass
     CALL GET_LOGICAL_VAL(TID_LOW_PASS,FILE_NAME,'TID_LOW_PASS',line,ierr)
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A15,L4)') 'TID_LOW_PASS= ',TID_LOW_PASS
# else
       WRITE(3,'(A15,L4)') 'TID_LOW_PASS= ',TID_LOW_PASS
# endif

# if defined (AMR)
      IF(NumGrid.GT.1)THEN
        DO K=1,NumGrid
        itmp1=mod(k/10,10)
        itmp2=mod(k,10)

        write(Ngrid(1:1),'(I1)')itmp1
        write(Ngrid(2:2),'(I1)')itmp2

        TMP_NAME = TRIM(RESULT_FOLDER)//'Grd'//Ngrid//'_'//'track.txt'
        OPEN(K+20,FILE=TRIM(TMP_NAME))
        ENDDO
      ENDIF
# endif
  ! end AMR

END SUBROUTINE READ_INPUT


# ifndef AMR
     subroutine preview
!--------------------------------------------------- 
!    This subroutine is used to preview
!    Called by                         
!       main 
!    Last update: 23/12/2010, Gangfeng Ma 
! fyshi gave an initial string to avoid unexpected characters
!  12/16/2011
!--------------------------------------------------
     use global
! fyshi add foam
# if defined (FOAM)
     USE MOD_FOAM
# endif
     implicit none
     integer :: i,j,k,I1,I2,I3,I4
     character(len=80) :: FDIR=''
     character(len=80) ::FILE_NAME=''
     character(LEN=80) :: file=''
     real(SP), dimension(Mloc,Nloc) :: Int2Flo
     real(SP), dimension(Mloc,Nloc,Kloc) :: Int2Flo3D

     ! file number
     Icount = Icount+1
   
     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     if(myid.eq.0) write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME     
# else
     write(*,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     write(3,102) 'Printing file No.',Icount,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
# endif
102  FORMAT(A20,I4,A14,F16.3,A2,F16.3)
100  FORMAT(15000E16.6)

     I1 = mod(Icount/1000,10)
     I2 = mod(Icount/100,10)
     I3 = mod(Icount/10,10)
     I4 = mod(Icount,10)

     write(FILE_NAME(1:1),'(I1)') I1
     write(FILE_NAME(2:2),'(I1)') I2
     write(FILE_NAME(3:3),'(I1)') I3
     write(FILE_NAME(4:4),'(I1)') I4

# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     open(5,file=TRIM(FDIR)//'time',position="append")
     write(5,*) TIME
     close(5)
# if defined (PARALLEL)
     endif
# endif

     if(Icount==1) then
       if(OUT_DEP) then
         file=TRIM(FDIR)//'depth'
         call putfile2D(file,Hc(1:Mloc,1:Nloc))
       endif
     endif

     if(OUT_ETA) then
       file=TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
       call putfile2D(file,Eta(1:Mloc,1:Nloc))
     endif

     if(OUT_U) then
       file=TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
       call putfile3D(file,U(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_V) then
       file=TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
       call putfile3D(file,V(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_W) then
       file=TRIM(FDIR)//'w_'//TRIM(FILE_NAME)
       call putfile3D(file,W(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_P) then
       file=TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
       call putfile3D(file,P(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_TKE) then
       file=TRIM(FDIR)//'tke_'//TRIM(FILE_NAME)
       call putfile3D(file,Tke(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_EPS) then
       file=TRIM(FDIR)//'eps_'//TRIM(FILE_NAME)
       call putfile3D(file,Eps(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_S) then
       file=TRIM(FDIR)//'prod_'//TRIM(FILE_NAME)
       call putfile3D(file,Prod_s(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_MU) then
       file=TRIM(FDIR)//'mu_'//TRIM(FILE_NAME)
       call putfile3D(file,CmuVt(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_A) then
       file=TRIM(FDIR)//'upwp_'//TRIM(FILE_NAME)
       call putfile3D(file,UpWp(1:Mloc,1:Nloc,1:Kloc))
     endif
# if defined (OBSTACLE)
      file=TRIM(FDIR)//'mask_'//TRIM(FILE_NAME)
       Int2Flo3D=Mask3D
       call putfile3D(file,Int2Flo3D(1:Mloc,1:Nloc,1:Kloc))
# endif

# if defined (BUBBLE)
     if(OUT_BUB) then
       file=TRIM(FDIR)//'b_'//TRIM(FILE_NAME)
       call putfile3D(file,Vbg(1:Mloc,1:Nloc,1:Kloc))
     endif
# endif

# if defined (SEDIMENT)
     if(OUT_F) then
       file=TRIM(FDIR)//'f_'//TRIM(FILE_NAME)
       call putfile3D(file,Conc(1:Mloc,1:Nloc,1:Kloc))
     endif
     if(OUT_T) then
       file=TRIM(FDIR)//'t_'//TRIM(FILE_NAME)
       call putfile2D(file,Taub(1:Mloc,1:Nloc))
     endif
     if(OUT_G) then
       file=TRIM(FDIR)//'g_'//TRIM(FILE_NAME)
       call putfile2D(file,Bed(1:Mloc,1:Nloc))
     endif
# endif

# if defined (SALINITY)
     if(OUT_SALI) then
       file=TRIM(FDIR)//'sali_'//TRIM(FILE_NAME)
       call putfile3D(file,Sali(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_TEMP) then
       file=TRIM(FDIR)//'temp_'//TRIM(FILE_NAME)
       call putfile3D(file,Temp(1:Mloc,1:Nloc,1:Kloc))
     endif

     if(OUT_RHO) then
       file=TRIM(FDIR)//'rho_'//TRIM(FILE_NAME)
       call putfile3D(file,Rho(1:Mloc,1:Nloc,1:Kloc))
     endif
# endif

! fyshi add foam
# if defined (FOAM)
     if(OUT_FOAM) then
       file=TRIM(FDIR)//'foam_'//TRIM(FILE_NAME)
       call putfile2D(file,EtaFoam(1:Mloc,1:Nloc))
     endif
# endif

     end subroutine preview

# endif

! end not define AMR 
# if defined (AMR)
SUBROUTINE PREVIEW_AMR(ng)
!--------------------------------------------------- 
!    This subroutine is used to preview
!    Called by                         
!       AMR_LOOP
! fyshi gave an initial string to avoid unexpected characters
!  12/16/2011
!--------------------------------------------------
     use global
     USE NESTING
! fyshi add foam
# if defined (FOAM)
     USE MOD_FOAM
# endif
     implicit none
     INTEGER,INTENT(IN)::ng
     integer :: i,j,k,I1,I2,I3,I4,itmp1,itmp2
     character(len=80) :: FDIR=''
     character(len=80) ::FILE_NAME=''
     character(LEN=80) :: file=''
     real(SP), dimension(Mloc,Nloc) :: Int2Flo
     CHARACTER(LEN=2)::Ngrid

     ! file number
    IF(ng.EQ.1)THEN
     Ioutput = Ioutput+1
    ENDIF

     ! results directory
     FDIR = TRIM(RESULT_FOLDER)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,102) 'Printing file No.',Ioutput,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     if(myid.eq.0) write(3,102) 'Printing file No.',Ioutput,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME     
# else
     write(*,102) 'Printing file No.',Ioutput,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
     write(3,102) 'Printing file No.',Ioutput,' TIME/TOTAL: ',TIME,'/',TOTAL_TIME
# endif
102  FORMAT(A20,I4,A14,F16.3,A2,F16.3)
100  FORMAT(15000E16.6)

     I1 = mod(Ioutput/1000,10)
     I2 = mod(Ioutput/100,10)
     I3 = mod(Ioutput/10,10)
     I4 = mod(Ioutput,10)

     write(FILE_NAME(1:1),'(I1)') I1
     write(FILE_NAME(2:2),'(I1)') I2
     write(FILE_NAME(3:3),'(I1)') I3
     write(FILE_NAME(4:4),'(I1)') I4

! grid
        itmp1=mod(ng/10,10)
        itmp2=mod(ng,10)

        write(Ngrid(1:1),'(I1)')itmp1
        write(Ngrid(2:2),'(I1)')itmp2

! track
     IF(NumGrid.GT.1)THEN
       WRITE(20+ng,*)TIME,MboxRef(ng),NboxRef(ng)
     ENDIF


# if defined (PARALLEL)
     if(myid.eq.0) then
# endif
     open(5,file=TRIM(FDIR)//'time',position="append")
     write(5,*) TIME
     close(5)
# if defined (PARALLEL)
     endif
# endif

     if(Ioutput==1) then
       if(OUT_DEP) then
         file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'depth'
         call putfile2D_AMR(file,Hc(1:Mloc,1:Nloc),ng)
       endif
     endif

     if(OUT_ETA) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'eta_'//TRIM(FILE_NAME)
       call putfile2D_AMR(file,Eta(1:Mloc,1:Nloc),ng)
     endif

     if(OUT_U) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'u_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,U(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_V) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'v_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,V(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_W) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'w_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,W(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_P) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'p_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,P(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_TKE) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'tke_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Tke(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_EPS) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'eps_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Eps(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_S) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'prod_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Prod_s(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_MU) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'mu_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,CmuVt(1:Mloc,1:Nloc,1:Kloc),ng)
     endif

     if(OUT_A) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'upwp_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,UpWp(1:Mloc,1:Nloc,1:Kloc),ng)
     endif
# if defined (BUBBLE)
     if(OUT_BUB) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'b_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Vbg(1:Mloc,1:Nloc,1:Kloc),ng)
     endif
# endif

# if defined (SEDIMENT)
     if(OUT_F) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'f_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Conc(1:Mloc,1:Nloc,1:Kloc),ng)
     endif
     if(OUT_T) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'t_'//TRIM(FILE_NAME)
       call putfile2D_AMR(file,Taub(1:Mloc,1:Nloc),ng)
     endif
     if(OUT_G) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'g_'//TRIM(FILE_NAME)
       call putfile2D_AMR(file,Bed(1:Mloc,1:Nloc),ng)
     endif
# endif

# if defined (SALINITY)
     if(OUT_SALI) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'sali_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Sali(1:Mloc,1:Nloc,1:Kloc),ng)

       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'temp_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Temp(1:Mloc,1:Nloc,1:Kloc),ng)

       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'rho_'//TRIM(FILE_NAME)
       call putfile3D_AMR(file,Rho(1:Mloc,1:Nloc,1:Kloc),ng)
     endif
# endif

! fyshi add foam
# if defined (FOAM)
     if(OUT_FOAM) then
       file=TRIM(FDIR)//'Grd'//Ngrid//'_'//'foam_'//TRIM(FILE_NAME)
       call putfile2D_AMR(file,EtaFoam(1:Mloc,1:Nloc),ng)
     endif
# endif

END SUBROUTINE PREVIEW_AMR
# endif
  ! end AMR

# if defined (PARALLEL)

# if defined (AMR)
    subroutine putfile2D_AMR(file,phi,ng)
    use global
    USE NESTING,ONLY : GridDimX,GridDimY
    implicit none
    integer,intent(in) :: ng
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(GridDimX(ng),GridDimY(ng)) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,iglob,jglob,len,n

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    ! put the data in master processor into the global var
    if(myid==0) then
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        phiglob(iglob,jglob) = Phi(i,j)
      enddo
      enddo
    endif

    ! collect data from other processors into the master processor
    len = Mloc*Nloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          phiglob(iglob,jglob) = xx(i,j)
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo       

    if(myid==0) then
      open(5,file=TRIM(file))
      do j = 1,GridDimY(ng)
        write(5,100) (phiglob(i,j),i=1,GridDimX(ng))
      enddo
      close(5)
    endif
100 FORMAT(15000f15.5)

    end subroutine putfile2D_AMR

# else
! not AMR

    subroutine putfile2D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(Mglob,Nglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,iglob,jglob,len,n

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    ! put the data in master processor into the global var
    if(myid==0) then
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        phiglob(iglob,jglob) = Phi(i,j)
      enddo
      enddo
    endif

    ! collect data from other processors into the master processor
    len = Mloc*Nloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          phiglob(iglob,jglob) = xx(i,j)
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo       

    if(myid==0) then
      open(5,file=TRIM(file))
      do j = 1,Nglob
        write(5,100) (phiglob(i,j),i=1,Mglob)
      enddo
      close(5)
    endif
100 FORMAT(15000f15.5)

    end subroutine putfile2D

# endif
  ! end if or not AMR

# if defined (AMR)
    subroutine putfile3D_AMR(file,phi,ng)
    use global
    USE NESTING,ONLY : GridDimX,GridDimY
    implicit none
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,philoc
    integer,intent(in) :: ng
    real(SP),dimension(GridDimX(ng),GridDimY(ng),Kglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,k,jk,iglob,jglob,kk,n,len,nreq,NKloc

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)

    NKloc = Nloc*Kloc

    ! put the data in master processor into the global var
    if(myid==0) then
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        kk = k-Nghost
        phiglob(iglob,jglob,kk) = Phi(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(philoc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      philoc(i,jk) = phi(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = Kbeg,Kend
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          kk = k-Nghost
          jk = (k-1)*Nloc+j
          phiglob(iglob,jglob,kk) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(philoc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do k = 1,Kglob
      do j = 1,GridDimY(ng)
        write(5,100) (phiglob(i,j,k),i=1,GridDimX(ng))
      enddo
      enddo
      close(5)
    endif
100 FORMAT(15000f15.5)

    deallocate(philoc)
    deallocate(xx)

    end subroutine putfile3D_AMR

# else
  ! no AMR

    subroutine putfile3D(file,phi)
    use global
    implicit none
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
    character(len=80) :: file
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,philoc
    real(SP),dimension(Mglob,Nglob,Kglob) :: phiglob
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,k,jk,iglob,jglob,kk,n,len,nreq,NKloc

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
          0,MPI_COMM_WORLD,ier)

    NKloc = Nloc*Kloc

    ! put the data in master processor into the global var
    if(myid==0) then
      do k = Kbeg,Kend
      do j = Jbeg,Jend
      do i = Ibeg,Iend
        iglob = i-Nghost
        jglob = j-Nghost
        kk = k-Nghost
        phiglob(iglob,jglob,kk) = Phi(i,j,k)
      enddo
      enddo
      enddo
    endif

    allocate(philoc(Mloc,NKloc))
    allocate(xx(Mloc,NKloc))

    do k = 1,Kloc
    do j = 1,Nloc
    do i = 1,Mloc
      jk = (k-1)*Nloc+j
      philoc(i,jk) = phi(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = Kbeg,Kend
        do j = Jbeg,Jend
        do i = Ibeg,Iend
          iglob = npxs(n+1)*(Iend-Ibeg+1)+i-Nghost
          jglob = npys(n+1)*(Jend-Jbeg+1)+j-Nghost
          kk = k-Nghost
          jk = (k-1)*Nloc+j
          phiglob(iglob,jglob,kk) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(philoc,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo

    if(myid.eq.0) then
      open(5,file=TRIM(file))
      do k = 1,Kglob
      do j = 1,Nglob
        write(5,100) (phiglob(i,j,k),i=1,Mglob)
      enddo
      enddo
      close(5)
    endif
100 FORMAT(15000E16.8)

    deallocate(philoc)
    deallocate(xx)

    end subroutine putfile3D
# endif 
  ! end if or not AMR


# else
  ! not parallel

# if defined (AMR)
     subroutine putfile2D_AMR(file,phi,ng)
# else
     subroutine putfile2D(file,phi)
# endif
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j
# if defined (AMR)
     integer :: ng
# endif

     open(5,file=trim(file))
# if defined(DEBUG)
     do j = 1,Nloc
       write(5,100) (phi(i,j),i=1,Mloc)
     enddo
# else
     do j = Jbeg,Jend
       write(5,100) (phi(i,j),i=Ibeg,Iend)
     enddo
# endif 
! end debug
     close(5)
 100 FORMAT(15000f15.5)

# if defined (AMR)
     end subroutine putfile2D_AMR
# else
     end subroutine putfile2D
# endif


# if defined (AMR)
     subroutine putfile3D_AMR(file,phi,ng)
# else
     subroutine putfile3D(file,phi)
# endif
     use global
     implicit none
     real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi
     character(len=80) :: file
     integer :: i,j,k
# if defined (AMR)
     integer :: ng
# endif

     open(5,file=trim(file))
# if defined (DEBUG)
     do k = 1,Kloc
     do j = 1,Nloc
       write(5,100) (phi(i,j,k),i=1,Mloc)
     enddo
     enddo
# else
     do k = Kbeg,Kend
     do j = Jbeg,Jend
       write(5,100) (phi(i,j,k),i=Ibeg,Iend)
     enddo
     enddo
# endif
! end debug
     close(5)
 100 FORMAT(15000f15.5)


# if defined (AMR)
     end subroutine putfile3D_AMR
# else
     end subroutine putfile3D
# endif


# endif
 ! end parallel

     subroutine read_bathymetry
!------------------------------------------------------ 
!    This subroutine is used to read bathymetry                                                
!    Called by  
!       main    
!    Last update: 21/12/2010, Gangfeng Ma      
!-----------------------------------------------------
     use global
# if defined (AMR)
     use nesting
# endif
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob,II,JJ
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG

# if defined (LANDSLIDE)
     real(SP) :: Hd,alpha0,L0,T,bl,wl,e,kb,kw,x0,x1,x2,xr,  &
                 Hb(Mloc,Nloc),xt,xt1,zt,zt1
# endif

# if defined (LANDSLIDE)
     Hd = 1.5
     T = 0.082
     bl = 0.395
     wl = 0.680
     e = 0.717
     alpha0 = 15.0*pi/180
     kb = 2.0*acosh(1.0/e)/bl
     kw = 2.0*acosh(1.0/e)/wl
     x0 = 0.551
     L0 = x0/cos(alpha0)-T*tan(alpha0)
     xr = 0.1
     x1 = (L0-bl/2.)*cos(alpha0)+xr
     x2 = (L0+bl/2.)*cos(alpha0)+xr
# endif

     ! read bathymetry at grid points
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 0.5
         enddo
         enddo
       else
! fyshi add standard depth file
         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

       ! find pernament dry points
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

# if defined (PARALLEL)
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# else
       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = i-Nghost
         jglob = j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo
# endif

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center
       if(ANA_BATHY) then

# if defined (LANDSLIDE)
         ! base bathymetry                                                                     
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xc(i)<=xr) then
             Hb(i,j) = -(xr-xc(i))*tan(alpha0)
           elseif(xc(i)<=Hd/tan(alpha0)+xr) then
             Hb(i,j) = (xc(i)-xr)*tan(alpha0)
           else
             Hb(i,j) = Hd
           endif
           ! temporarily no runup                                                              
           Hc(i,j) = max(0.01,Hb(i,j))
         enddo
         enddo

         do j = Jbeg,Jend
         do i = Ibeg,Iend
           if(xc(i)<=x1.or.xc(i)>=x2.or.yc(j)>=wl/2.) then
             Hc(i,j) = Hb(i,j)
           else
             iter = 1
             xt = (xc(i)-xr)/cos(alpha0)-L0
             zt = T/(1-e)*(1.0/cosh(kb*xt)/cosh(kw*yc(j))-e)
 60          xt1 = ((xc(i)-xr)/cos(alpha0)-zt*tan(alpha0))-L0
             zt1 = T/(1-e)*(1.0/cosh(kb*xt1)/cosh(kw*yc(j))-e)
             if(abs(zt1-zt)/abs(zt)>1.e-8) then
               iter = iter+1
               if(iter>20) write(*,*) 'too many iterations!'
               zt = zt1
               goto 60
             endif
             Hc(i,j) = Hb(i,j)-max(0.0,zt1)/cos(alpha0)
           endif
         enddo
         enddo
# endif

       else ! not analytical bathymetry, read from depth file

# if defined (AMR)
      IF(NumGrid>1)THEN
         open(5,file=TRIM(Depth_FULL_File),status='old')
         do j = 1,Ndep
           read(5,*) (Depth_full(i,j),i=1,Mdep)
         enddo
         close(5)
        IF(NESTING_X1D)THEN
         DO J=1,Nglob
         DO I=1,Mglob
           II=(I-1)*MainGridRatio+1
           JJ=J
           HCG(I,J)=Depth_full(II,JJ)
         ENDDO
         ENDDO
        ELSE
         DO J=1,Nglob
         DO I=1,Mglob
           II=(I-1)*MainGridRatio+1
           JJ=(J-1)*MainGridRatio+1
           HCG(I,J)=Depth_full(II,JJ)
         ENDDO
         ENDDO
        ENDIF ! end nesting_x1d
      ELSE
# endif
         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo
         close(5)
# if defined (AMR) 
      ENDIF ! end numgrid>1
# endif

# if defined (PARALLEL)
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# else
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = i-Nghost
           jglob = j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
# endif

       endif

     endif

     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-10.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! collect data into ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hc(1,j)
     Hfx(Mloc1,j) = Hc(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hc(i,1)
     Hfy(i,Nloc1) = Hc(i,Nloc)
     enddo

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     ! save initial bathymetry
     Hc0 = Hc


     end subroutine read_bathymetry

  

# if defined (COUPLING)
SUBROUTINE READ_NESTING_FILE
!---------------------------------------------------
!    This subroutine is used to read nesting file at the first time
!    Called by 
!       main
!    Last update: 05/15/2013, fyshi
!---------------------------------------------------
     use global
     use input_util
     implicit none
     INTEGER :: I,J,K
       OPEN(11,FILE=TRIM(COUPLING_FILE))
         READ(11,*)  ! title
         READ(11,*)  ! boundary info
! boundary basic info including point number of coupling, start point, etc
! east
         READ(11,*)  ! east
         READ(11,*) N_COUPLING_EAST,J_START_EAST
! west 
         READ(11,*)  ! west
         READ(11,*) N_COUPLING_WEST,J_START_WEST
! south 
         READ(11,*)  ! south
         READ(11,*) N_COUPLING_SOUTH,I_START_SOUTH
! north 
         READ(11,*)  ! north
         READ(11,*) N_COUPLING_NORTH,I_START_NORTH

! read time and variable at the first level

         READ(11,*) ! time start title
         READ(11,*) TIME_COUPLING_1 
! initialize time_2
         TIME_COUPLING_2 = TIME_COUPLING_1

! east
         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               V_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               W_COUPLING_EAST(N_COUPLING_EAST,Kglob,2),&
               Z_COUPLING_EAST(N_COUPLING_EAST,2), &
               P_COUPLING_EAST(N_COUPLING_EAST,Kglob,2), &
               S_COUPLING_EAST(N_COUPLING_EAST,Kglob,2), &
               T_COUPLING_EAST(N_COUPLING_EAST,Kglob,2) &
               )
             READ(11,*)   ! east
             READ(11,119)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119)((U_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((V_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((W_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((P_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((S_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119)((T_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
!   initialize first step
             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_east
119      FORMAT(5E16.6)

! west
         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               V_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               W_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               Z_COUPLING_WEST(N_COUPLING_WEST,2),&
               P_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               S_COUPLING_WEST(N_COUPLING_WEST,Kglob,2),&
               T_COUPLING_WEST(N_COUPLING_WEST,Kglob,2) &
               )
             READ(11,*)   ! west
             READ(11,119)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119)((U_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((V_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((W_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((P_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((S_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119)((T_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
!   initialize first step
             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_west
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               W_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,2),&
               P_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               S_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2),&
               T_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,2) &
               )
             READ(11,*)   ! south
             READ(11,119)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119)((U_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((V_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((W_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((P_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((S_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119)((T_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
!   initialize first step
             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_south
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               W_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,2), &
               P_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2), &
               S_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2), &
               T_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,2) &
               )
             READ(11,*)   ! north
             READ(11,119)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119)((U_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((V_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((W_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((P_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((S_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119)((T_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)

!   initialize first step
             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
             S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)
         ELSE
             READ(11,*)

         ENDIF ! n_coupling_north


! specify boundary start points

! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
      Kstart_WEST=J_START_WEST+Nghost -npy*Nglob/py
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1 -npy*Nglob/py
      IF((Kstart_WEST>Nghost.AND.Kstart_WEST<Nloc-Nghost+1).OR.&
     (Kend_WEST<Nloc-Nghost+1.AND.Kend_WEST>Nghost))THEN
       IF(Kstart_WEST<Nghost+1)THEN
         Kshift_WEST= -Kstart_WEST+Nghost+1
         Kstart_WEST=Nghost+1
       ELSE
         Kshift_WEST=-(Kstart_WEST-Nghost)+1
       ENDIF
       IF(Kend_WEST>Nloc-Nghost)THEN
         Kend_WEST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_WEST=.TRUE.
      ELSE
       IF(Kstart_WEST<=Nghost.AND.Kend_WEST>=Nloc-Nghost+1)THEN
         Kshift_WEST=-Kstart_WEST+Nghost+1
         Kstart_WEST = Nghost+1
         Kend_WEST = Nloc-Nghost
         IN_DOMAIN_WEST=.TRUE.
       ELSE
         IN_DOMAIN_WEST=.FALSE.
       ENDIF
      ENDIF
! check print*,myid,Kshift_WEST,Kstart_WEST,Kend_WEST
     endif
# else
      Kstart_WEST=J_START_WEST+Nghost
      Kend_WEST = J_START_WEST+Nghost+N_COUPLING_WEST-1
      Kshift_WEST = -(Kstart_WEST-Nghost)+1
      IN_DOMAIN_WEST = .TRUE.
# endif

   ENDIF

! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
      Kstart_EAST=J_START_EAST+Nghost -npy*Nglob/py
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1 -npy*Nglob/py
      IF((Kstart_EAST>Nghost.AND.Kstart_EAST<Nloc-Nghost+1).OR.&
     (Kend_EAST<Nloc-Nghost+1.AND.Kend_EAST>Nghost))THEN
       IF(Kstart_EAST<Nghost+1)THEN
         Kshift_EAST= -Kstart_EAST+Nghost+1
         Kstart_EAST=Nghost+1
       ELSE
         Kshift_EAST=-(Kstart_EAST-Nghost)+1
       ENDIF
       IF(Kend_EAST>Nloc-Nghost)THEN
         Kend_EAST=Nloc-Nghost
       ENDIF
       IN_DOMAIN_EAST=.TRUE.
      ELSE

       IF(Kstart_EAST<=Nghost.AND.Kend_EAST>=Nloc-Nghost+1)THEN
         Kshift_EAST=-Kstart_EAST+Nghost+1
         Kstart_EAST = Nghost+1
         Kend_EAST = Nloc-Nghost
         IN_DOMAIN_EAST=.TRUE.
       ELSE
         IN_DOMAIN_EAST=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_EAST=J_START_EAST+Nghost
      Kend_EAST = J_START_EAST+Nghost+N_COUPLING_EAST-1
      Kshift_EAST = -(Kstart_EAST-Nghost)+1
      IN_DOMAIN_EAST = .TRUE.
# endif
    ENDIF

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
      Kstart_SOUTH=I_START_SOUTH+Nghost -npx*Mglob/px
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1 -npx*Mglob/px
      IF((Kstart_SOUTH>Nghost.AND.Kstart_SOUTH<Mloc-Nghost+1).OR.&
     (Kend_SOUTH<Mloc-Nghost+1.AND.Kend_SOUTH>Nghost))THEN
       IF(Kstart_SOUTH<Nghost+1)THEN
         Kshift_SOUTH= -Kstart_SOUTH+Nghost+1
         Kstart_SOUTH=Nghost+1
       ELSE
         Kshift_SOUTH=-(Kstart_SOUTH-Nghost)+1
       ENDIF
       IF(Kend_SOUTH>Mloc-Nghost)THEN
         Kend_SOUTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_SOUTH=.TRUE.
      ELSE

       IF(Kstart_SOUTH<=Nghost.AND.Kend_SOUTH>=Mloc-Nghost+1)THEN
         Kshift_SOUTH=-Kstart_SOUTH+Nghost+1
         Kstart_SOUTH = Nghost+1
         Kend_SOUTH = Mloc-Nghost
         IN_DOMAIN_SOUTH=.TRUE.
       ELSE
         IN_DOMAIN_SOUTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_SOUTH=I_START_SOUTH+Nghost
      Kend_SOUTH = I_START_SOUTH+Nghost+N_COUPLING_SOUTH-1
      Kshift_SOUTH = -(Kstart_SOUTH-Nghost)+1
      IN_DOMAIN_SOUTH = .TRUE.
# endif
   ENDIF

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
      Kstart_NORTH=I_START_NORTH+Nghost -npx*Mglob/px
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1 -npx*Mglob/px
      IF((Kstart_NORTH>Nghost.AND.Kstart_NORTH<Mloc-Nghost+1).OR.&
     (Kend_NORTH<Mloc-Nghost+1.AND.Kend_NORTH>Nghost))THEN
       IF(Kstart_NORTH<Nghost+1)THEN
         Kshift_NORTH= -Kstart_NORTH+Nghost+1
         Kstart_NORTH=Nghost+1
       ELSE
         Kshift_NORTH=-(Kstart_NORTH-Nghost)+1
       ENDIF
       IF(Kend_NORTH>Mloc-Nghost)THEN
         Kend_NORTH=Mloc-Nghost
       ENDIF
       IN_DOMAIN_NORTH=.TRUE.
      ELSE

       IF(Kstart_NORTH<=Nghost.AND.Kend_NORTH>=Mloc-Nghost+1)THEN
         Kshift_NORTH=-Kstart_NORTH+Nghost+1
         Kstart_NORTH = Nghost+1
         Kend_NORTH = Mloc-Nghost
         IN_DOMAIN_NORTH=.TRUE.
       ELSE
         IN_DOMAIN_NORTH=.FALSE.
       ENDIF
      ENDIF

    endif
# else
      Kstart_NORTH=I_START_NORTH+Nghost
      Kend_NORTH = I_START_NORTH+Nghost+N_COUPLING_NORTH-1
      Kshift_NORTH = -(Kstart_NORTH-Nghost)+1
      IN_DOMAIN_NORTH = .TRUE.
# endif
   ENDIF

END SUBROUTINE READ_NESTING_FILE
# endif


