! -----------------------------------------------
!    This is the MAIN program
!    It calls MAIN (ARM_LOOP) or 
!    Gengfeng Ma's original main program (SINGLE_GRID_LOOP)
!    last update 05/24/2013, fyshi
! ------------------------------------------------

PROGRAM MASTER

     USE GLOBAL
# if defined (AMR)
     USE NESTING
# endif
     IMPLICIT NONE
     REAL(SP) :: tbegin,tend
     INTEGER :: Kgrid,Isubgrid
     INTEGER :: Ksubgrid_internal,NumInternal


# if defined (PARALLEL)
     call MPI_INIT ( ier )
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ier)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,NumP,ier)
# endif

     CALL READ_INPUT

# if defined (AMR)
     CALL ALLOCATE_VAR_TILE      ! for all grids
# endif

# if defined (PARALLEL)
     CALL PARALLEL_CARTESIAN
# endif

     CALL INDEX_LOCAL(1)   ! grid(1): mloc,nloc,ibeg, iend, etc.

# if defined (DEBUG)
print*,'index 1',mloc,nloc,kloc
# endif

     CALL ALLOCATE_VARIABLES ! using maximum dimension

# if defined (DEBUG)
print*,'allocate variables with max mloc,nloc,kloc',mloc,nloc,kloc
# endif

# if defined (FOAM)
     CALL ALLOCATE_FOAM
# endif


     ! read bathymetry
     call read_bathymetry  ! for grid(1) only

# if defined (DEBUG)
print*,'read bathymetry OK'
# endif

     CALL INITIALIZATION   ! for grid(1) only

# if defined (DEBUG)
print*,'init OK with Mloc,Nloc,Kloc=',Mloc,Nloc,Kloc
# endif

     ! generate grids
     call generate_grid   ! for grid(1) only


# if defined (FOAM)
     CALL INITIALIZATION_FOAM
# endif

# if defined (AMR)
     CALL INITIAL_GRID(1)   ! assign initial to grid(1)
# endif

# if defined (COUPLING)
     call read_nesting_file
# endif

# if defined (DEBUG)
print*,'read nesting file OK'
# endif


# if defined (COUPLING)
     TIME=TIME_COUPLING_1
# endif
! time integration

     ! record wall time
     call WALL_TIME_SECS(tbegin)

# if defined(AMR)
   START_GRID(1)=.TRUE.
   START_COUNT=1
# endif

# if defined(AMR)

       CALL AMR_LOOP

# else
       CALL SINGLE_GRID_LOOP

# endif 
! end AMR


# if defined (COUPLING)
     CLOSE(11)
# endif

     ! record wall time at the end
     call WALL_TIME_SECS(tend)

# if defined (PARALLEL)
     if(myid.eq.0) write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     if(myid.eq.0) write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     if (myid.eq.0) WRITE(*,*)'Normal Termination!'
     if (myid.eq.0) WRITE(3,*)'Normal Termination!'
# else
     write(*,*) 'Simulation takes',tend-tbegin,'seconds'
     write(3,*) 'Simulation takes',tend-tbegin,'seconds'
     WRITE(*,*)'Normal Termination!'
     WRITE(3,*)'Normal Termination!'
# endif

# if defined (PARALLEL)
     call MPI_FINALIZE ( ier )
# endif



!     GRID(1) % DX_TILE(1,1) = 1
!     GRID(2) % DX_TILE(1,1) =2

!     PRINT*,GRID(1)%DX_TILE(1,1)
!     PRINT*,GRID(2)%DX_TILE(1,1)

END PROGRAM MASTER

# ifndef AMR
  ! not AMR
! -----------------------------------------------
!  This is Gengfeng's original main program
!  I moved his main from nhwave.F here
!  last update 08/24/2013, fyshi
! -----------------------------------------------
SUBROUTINE SINGLE_GRID_LOOP
     USE GLOBAL
     IMPLICIT NONE
     integer :: j,Istage
     REAL(SP) :: tbegin,tend
     INTEGER :: Kgrid,Isubgrid
     INTEGER :: Ksubgrid_internal,NumInternal

     ! time integration
     do while (TIME<TOTAL_TIME.and.RUN_STEP<SIM_STEPS)


# if defined (COUPLING)
! fyshi added coupling(with a large domain model) 
! condition option 05/15/2013
         CALL OneWayCoupling
# endif

# if defined (OBSTACLE)
         CALL READ_3D_MASK
# endif

       ! time step     
       call estimate_dt

# if defined (LANDSLIDE)
       ! run landslide-generated tsunami
       if(trim(SlideType)=='RIGID') then
         call update_bathymetry
       endif
# endif

       ! update boundary conditions       
       call update_wave_bc

       ! update mask
       call update_mask


       ! update vars
       call update_vars


       ! SSP Runge-Kutta time stepping
       do Istage = 1,It_Order

         ! well-balanced source terms
         call source_terms

         ! fluxes at cell faces
         call fluxes(1)

         ! update all variables
         call eval_duvw(Istage,1)

         ! sponge layer
         if(SPONGE_ON) then
           call sponge_damping
         endif

         ! turbulence model
         if(VISCOUS_FLOW) call eval_turb(Istage,1)


# if defined (SALINITY)
         ! update salinity
         call eval_sali(Istage,1)

         ! update density
         call eval_dens  
# endif


# if defined (BUBBLE)
         if(TIME>=TIM_B) then
           ! bubble rise velocity
           call bslip_velocity

           ! update bubble concentration
           call eval_bub(Istage,1)
         endif
# endif


# if defined (SEDIMENT)
         if(TIME>=TIM_Sedi) then
           ! settling velocity
           call settling_velocity

           ! update sediment concentration
           call eval_sedi(Istage)

           if(BED_CHANGE.and.TIME>=10.0) then
             ! update bed elevation
             call update_bed(Istage)          
           endif

           if(COUPLE_FS) then
             ! update mixture density
             call eval_dens
           endif
         endif
# endif

# if defined (COUPLING)
! fyshi added nesting option 05/15/2013
         CALL OneWayCoupling
# endif

       enddo

# if defined (BALANCE2D)
       ! evaluate momentum balance in cross-shore
       call eval_balance
# endif

       ! wave average quantities
       if(WAVE_AVERAGE_ON) then
         call wave_average
       endif

       ! screen output
       Screen_Count = Screen_Count+dt
       if(Screen_Count>=Screen_Intv) then
         Screen_Count = Screen_Count-Screen_Intv
         call statistics
       endif

       ! probe output to files
       if(NSTAT>0) then
         Plot_Count_Stat = Plot_Count_Stat+dt
         if(Plot_Count_Stat>=Plot_Intv_Stat) then
           Plot_Count_Stat=Plot_Count_Stat-Plot_Intv_Stat
           call probes
         endif
       endif

       ! field output to files
       if(TIME>=Plot_Start) then
         Plot_Count = Plot_Count+dt
         if(Plot_Count>=Plot_Intv) then
           Plot_Count=Plot_Count-Plot_Intv
           call preview
         endif
       endif
      
     end do


END SUBROUTINE SINGLE_GRID_LOOP

# endif
   ! end not AMR

# if defined(AMR)
! -----------------------------------------------
!    This is the whole structure of AMR
!    last update 08/24/2013 fyshi
! -----------------------------------------------
SUBROUTINE AMR_LOOP
     USE GLOBAL
     USE NESTING
     IMPLICIT NONE
     REAL(SP) :: tbegin,tend
     INTEGER :: Kgrid,Isubgrid
     INTEGER :: Ksubgrid_internal,NumInternal
     INTEGER :: I,J

   DO WHILE (TIME<TOTAL_TIME)

# if defined (DEBUG)
print*,'TIME',TIME,'TOTAL',TOTAL_TIME
# endif

! every grid calculation

     DO Kgrid=1,NumGrid
!     DO Kgrid=1,1

      NumInternal=RATIO_SPACING(Kgrid)
      IF(Kgrid.GT.1) DT=DT_grid1/RATIO_SPACING(Kgrid)

      DO Ksubgrid_internal=1,NumInternal

# if defined (DEBUG)
print*, 'GRID ',Kgrid, '----------'
# endif


! subgrid (initial)
       IF(Kgrid.LT.NumGrid)THEN
         IF(START_GRID(Kgrid+1))THEN

           CALL INDEX_LOCAL(Kgrid+1)

# if defined (DEBUG)
print*,'1'
# endif
           
           IF(IS_GRID_MOVE(Kgrid+1))THEN
             CALL MOVING_GRID(Kgrid+1)
           ENDIF

# if defined (DEBUG)
print*,'2'
# endif


           CALL INITIAL_GRID(Kgrid+1)

# if defined (DEBUG)
print*,'3'
# endif

           CALL GET_NESTING_BC(Kgrid+1,1)


# if defined (DEBUG)
print*,'4'
# endif

           IS_GRID_MOVE(Kgrid+1)=.FALSE.

         ENDIF
       ENDIF


! calc grid
      IF(START_GRID(Kgrid))THEN
       CALL INDEX_LOCAL(Kgrid)

# if defined (DEBUG)
print*,'START CALC GRID',kgrid,'Subgrid step',Ksubgrid_internal
print*,'DIM',Mloc,Nloc,Kloc
# endif

!  call re_initialization

       CALL CALC_GRID(Kgrid,Ksubgrid_internal,NumInternal)

!       IF(Kgrid.EQ.1)DT_grid1=DT
      ENDIF

! subgrid (boundary)
       IF(Kgrid.LT.NumGrid)THEN
         IF(START_GRID(Kgrid+1))THEN
           CALL INDEX_LOCAL(Kgrid+1)
           CALL GET_NESTING_BC(Kgrid+1,2)
         ENDIF
       ENDIF


      ENDDO ! end Ksubgrid_internal

! print out
     IF(Kgrid.EQ.1)THEN
      PLOT_COUNT=PLOT_COUNT+DT_grid1

      IF(PLOT_COUNT>=PLOT_INTV)THEN
         PLOT_COUNT=PLOT_COUNT-PLOT_INTV
         DO Isubgrid=1,NumGrid
           IsPrint(Isubgrid)=.TRUE.
         ENDDO
      ENDIF
     ENDIF

      IF(IsPrint(Kgrid))THEN
        CALL INDEX_LOCAL(Kgrid)  ! index back
        CALL PREVIEW_AMR(Kgrid)
        IsPrint(Kgrid)=.FALSE.
      ENDIF

! screen
     IF(Kgrid.EQ.1)THEN
       Screen_Count = Screen_Count+dt_grid1
       if(Screen_Count>=Screen_Intv) then
         Screen_Count = Screen_Count-Screen_Intv
         call statistics
       endif
     ENDIF

     ENDDO  ! end kgrid

     IF(START_COUNT.LT.NumGrid)THEN
       START_COUNT=START_COUNT+1
       START_GRID(START_COUNT)=.TRUE.
      ENDIF

! --------------------------------------------- 

   END DO  ! end time

END SUBROUTINE AMR_LOOP
# endif
! end AMR

