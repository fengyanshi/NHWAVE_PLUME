! -----------------------------------------------
!    This is the main subroutine to calculate everything
!    for a specified grid. It corresponds to
!    Gengfeng's main program. NOTE: nhwave.F is gangfeng's 
!    original program. There're minor changes for AMR
!    CALLED BY MASTER
!    fyshi 05/24/2013
! ------------------------------------------------

# if defined (AMR)

SUBROUTINE CALC_GRID(ng,NestStep,NestTotal)

     USE GLOBAL
     USE NESTING
# if defined (FOAM)
     USE MOD_FOAM
# endif
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng,NestStep,NestTotal
     real(SP), parameter :: Dmin = 0.04,cmiu = 0.09
     ! this should be consistent with in turb.F

     integer :: i,j,k,Istage
     real(SP) :: tbegin,tend

! initialized necessary variables

    DO J=1,Nloc
    DO I=1,Mloc
     Eta(I,J)=Eta_Tile(I,J,ng)
     MASK(I,J)=MASK_Tile(I,J,ng)
     MASK9(I,J)=MASK9_Tile(I,J,ng)
     Hc(I,J)=Depth_Tile(I,J,ng)
     Hfx(I,J)=DepthX_Tile(I,J,ng)
     Hfy(I,J)=DepthY_Tile(I,J,ng)
    ENDDO
    ENDDO
    DO K=1,Kloc
    DO J=1,Nloc
    DO I=1,Mloc
     U(I,J,K)=U_Tile(I,J,K,ng)
     V(I,J,K)=V_Tile(I,J,K,ng)
     W(I,J,K)=W_Tile(I,J,K,ng)
     P(I,J,K)=P_Tile(I,J,K,ng)
# if defined (SALINITY)
     Sali(I,J,K)=S_Tile(I,J,K,ng)
     Temp(I,J,K)=T_Tile(I,J,K,ng)
# endif
     Rho(I,J,K)=Rho_Tile(I,J,K,ng)
    ENDDO
    ENDDO
    ENDDO


    DO K=1,Kloc1
    DO J=1,Nloc
    DO I=1,Mloc
      Omega(I,J,K)=Omega_Tile(I,J,K,ng)
    ENDDO
    ENDDO
    ENDDO

     DX=DX_Tile(ng)
     DY=DY_Tile(ng)

     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

! other needed variables
     D=Hc+Eta

     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       DU(I,J,K)=D(I,J)*U(I,J,K)*MASK(I,J)
       DV(I,J,K)=D(I,J)*V(I,J,K)*MASK(I,J)
       DW(I,J,K)=D(I,J)*W(I,J,K)*MASK(I,J)
# if defined (SALINITY)
       DSali(I,J,K)=D(I,J)*Sali(I,J,K)
       DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif
     ENDDO
     ENDDO
     ENDDO

  IF(VISCOUS_FLOW)THEN
     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       Tke(I,J,K)=Tke_Tile(I,J,K,ng)
       Eps(I,J,K)=Eps_Tile(I,J,K,ng)
       DTke(I,J,K)=D(I,J)*Tke(I,J,K)
       DEps(I,J,K)=D(I,J)*Eps(I,J,K)
     ENDDO
     ENDDO
     ENDDO

    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
    enddo
    enddo
    enddo

    CmuHt = CmuVt   ! for IHturb>=10

   
   ELSE
    do i = 1,Mloc
    do j = 1,Nloc
    do k = 1,Kloc       
       CmuVt(I,J,K) = Cmut_min
       CmuHt(I,J,K) = Cmut_min
    enddo
    enddo
    enddo

   ENDIF ! end viscous_flow


# if defined (DEBUG)
print*,'m1'
# endif


# if defined (COUPLING)
! fyshi added coupling(with a large domain model) 
! condition option 05/15/2013
       IF(ng==1)THEN
         CALL OneWayCoupling
       ENDIF
# endif

       IF(ng.GT.1)THEN
         CALL USE_NESTING_BC(NestStep,NestTotal,ng)
       ENDIF


# if defined (DEBUG)
print*,'m2'
# endif

       ! time step     

     IF(ng==1)THEN
       call estimate_dt
     ENDIF

# if defined (DEBUG)
print*,'m3'
# endif

# if defined (LANDSLIDE)
       ! run landslide-generated tsunami
       call update_bathymetry
# endif

       ! update boundary conditions   
     IF(ng==1)THEN    
       call update_wave_bc
     ENDIF

       ! update mask
       call update_mask

       ! update vars
       call update_vars

       ! SSP Runge-Kutta time stepping
       do Istage = 1,It_Order

         ! well-balanced source terms
         call source_terms

# if defined (DEBUG)
print*,'m4',ng
# endif

         ! fluxes at cell faces
         call fluxes(ng)

# if defined (DEBUG)
print*,'m5',ng
# endif

         ! update all variables
         call eval_duvw(Istage,ng)

       IF(ng==1)THEN
         ! sponge layer
         if(SPONGE_ON) then
           call sponge_damping
         endif
       ENDIF

# if defined (DEBUG)
print*,'m6',ng
# endif

! turbulence model was here before
         ! turbulence model
         if(VISCOUS_FLOW) call eval_turb(Istage,ng)

# if defined (DEBUG)
print*,'m7',ng
# endif

# if defined (SALINITY)
         ! update salinity
         call eval_sali(Istage,ng)

# if defined (DEBUG)
print*,'m8'
# endif

         ! update density
         call eval_dens  
# endif

# if defined (DEBUG)
print*,'m9'
# endif

# if defined (BUBBLE)
         if(TIME>=TIM_B) then
           ! bubble rise velocity
           call bslip_velocity

           ! update bubble concentration
           call eval_bub(Istage,ng)
         endif
# endif

# if defined (SEDIMENT)
         if(TIME>=TIM_Sedi) then
           ! settling velocity
           call settling_velocity

           ! update sediment concentration
           call eval_sedi(Istage)

           if(BED_CHANGE) then
             ! update bed elevation
             call update_bed(Istage)          
           endif

           if(COUPLE_FS) then
             ! update mixture density
             call eval_dens
           endif
         endif
# endif

# if defined (FOAM)
       CALL FOAM_FLUX
       CALL FOAM_UPDATE
       CALL FOAM_BC
# endif

! -------? can we put this outside r-k
# if defined (COUPLING)
! fyshi added nesting option 05/15/2013
       IF(ng==1)THEN
         CALL OneWayCoupling
       ENDIF
# endif

     IF(ng.GT.1)THEN
       CALL USE_NESTING_BC(NestStep,NestTotal,ng)
     ENDIF
! -------? end

       enddo    ! end R-K

! save for next time step
       DO J=1,Nloc
       DO I=1,Mloc
         Eta_TILE(I,J,ng) = Eta(I,J)
         MASK_TILE(I,J,ng) = MASK(I,J)
         MASK9_TILE(I,J,ng) = MASK9(I,J)
       ENDDO
       ENDDO
       DO K=1,Kloc
       DO J=1,Nloc
       DO I=1,Mloc
         U_TILE(I,J,K,ng) = U(I,J,K)
         V_TILE(I,J,K,ng) = V(I,J,K)
         W_TILE(I,J,K,ng) = W(I,J,K)
         P_TILE(I,J,K,ng) = P(I,J,K)
# if defined (SALINITY)
         S_TILE(I,J,K,ng) = Sali(I,J,K)
         T_TILE(I,J,K,ng) = Temp(I,J,K)
# endif
         Rho_TILE(I,J,K,ng) = Rho(I,J,K)
       ENDDO
       ENDDO
       ENDDO
    
       IF(VISCOUS_FLOW)THEN
       DO K=1,Kloc
       DO J=1,Nloc
       DO I=1,Mloc
         Tke_TILE(I,J,K,ng) = Tke(I,J,K)
         Eps_TILE(I,J,K,ng) = Eps(I,J,K)
       ENDDO
       ENDDO
       ENDDO
       ENDIF

       DO K=1,Kloc1
       DO J=1,Nloc
       DO I=1,Mloc
         Omega_TILE(I,J,K,ng) = Omega(I,J,K)
       ENDDO
       ENDDO
       ENDDO


       IF(ng.eq.1)THEN

       ! probe output to files
       if(NSTAT>0) then
         Plot_Count_Stat = Plot_Count_Stat+dt
         if(Plot_Count_Stat>=Plot_Intv_Stat) then
           Plot_Count_Stat=Plot_Count_Stat-Plot_Intv_Stat
           call probes
         endif
       endif

      ENDIF ! end ng=1      

END SUBROUTINE CALC_GRID

# endif
! end AMR


# if defined (COUPLING)
!-------------------------------------------------------------------
!   This subroutine is used to pass coupling variables into ghost cells                                                         
!   Called by
!      MAIN
!   Update: 05/15/2013 Fengyan Shi, University of Delaware                                       
!-------------------------------------------------------------------
SUBROUTINE OneWayCoupling
    USE GLOBAL
    IMPLICIT NONE
    INTEGER :: I,J,K
    REAL(SP) :: tmp1,tmp2
    LOGICAL, SAVE :: FirstReadCoupling = .TRUE.

119      FORMAT(5E16.6)  ! this is a fixed format for I/O
 
! determine time slot

    IF(TIME>TIME_COUPLING_1.AND.TIME>TIME_COUPLING_2) THEN
         TIME_COUPLING_1=TIME_COUPLING_2
         
         READ(11,*,END=120) ! time coupling
         READ(11,*,END=120) TIME_COUPLING_2 
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             READ(11,*,END=120)   ! east

             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)

             READ(11,119,END=120)(Z_COUPLING_EAST(I,2),I=1,N_COUPLING_EAST)
             READ(11,119,END=120)((U_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((S_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)
             READ(11,119,END=120)((T_COUPLING_EAST(I,J,2),I=1,N_COUPLING_EAST),J=1,Kglob)

            IF (FirstReadCoupling)THEN
             U_COUPLING_EAST(:,:,1)=U_COUPLING_EAST(:,:,2)
             V_COUPLING_EAST(:,:,1)=V_COUPLING_EAST(:,:,2)
             W_COUPLING_EAST(:,:,1)=W_COUPLING_EAST(:,:,2)
             Z_COUPLING_EAST(:,1)=Z_COUPLING_EAST(:,2)
             P_COUPLING_EAST(:,:,1)=P_COUPLING_EAST(:,:,2)
             S_COUPLING_EAST(:,:,1)=S_COUPLING_EAST(:,:,2)
             T_COUPLING_EAST(:,:,1)=T_COUPLING_EAST(:,:,2)
            ENDIF

         ELSE
             READ(11,*,END=120)   ! east            
         ENDIF
! west
         IF(N_COUPLING_WEST.GT.0)THEN
             READ(11,*,END=120)   ! west

             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)

             READ(11,119,END=120)(Z_COUPLING_WEST(I,2),I=1,N_COUPLING_WEST)
             READ(11,119,END=120)((U_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((S_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
             READ(11,119,END=120)((T_COUPLING_WEST(I,J,2),I=1,N_COUPLING_WEST),J=1,Kglob)
            IF (FirstReadCoupling)THEN
             U_COUPLING_WEST(:,:,1)=U_COUPLING_WEST(:,:,2)
             V_COUPLING_WEST(:,:,1)=V_COUPLING_WEST(:,:,2)
             W_COUPLING_WEST(:,:,1)=W_COUPLING_WEST(:,:,2)
             Z_COUPLING_WEST(:,1)=Z_COUPLING_WEST(:,2)
             P_COUPLING_WEST(:,:,1)=P_COUPLING_WEST(:,:,2)
             S_COUPLING_WEST(:,:,1)=S_COUPLING_WEST(:,:,2)
             T_COUPLING_WEST(:,:,1)=T_COUPLING_WEST(:,:,2)
            ENDIF

         ELSE
             READ(11,*,END=120)   ! west            
         ENDIF
! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             READ(11,*,END=120)   ! south

             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)

             READ(11,119,END=120)(Z_COUPLING_SOUTH(I,2),I=1,N_COUPLING_SOUTH)
             READ(11,119,END=120)((U_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((S_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             READ(11,119,END=120)((T_COUPLING_SOUTH(I,J,2),I=1,N_COUPLING_SOUTH),J=1,Kglob)
            IF (FirstReadCoupling)THEN
             U_COUPLING_SOUTH(:,:,1)=U_COUPLING_SOUTH(:,:,2)
             V_COUPLING_SOUTH(:,:,1)=V_COUPLING_SOUTH(:,:,2)
             W_COUPLING_SOUTH(:,:,1)=W_COUPLING_SOUTH(:,:,2)
             Z_COUPLING_SOUTH(:,1)=Z_COUPLING_SOUTH(:,2)
             P_COUPLING_SOUTH(:,:,1)=P_COUPLING_SOUTH(:,:,2)
             S_COUPLING_SOUTH(:,:,1)=S_COUPLING_SOUTH(:,:,2)
             T_COUPLING_SOUTH(:,:,1)=T_COUPLING_SOUTH(:,:,2)
            ENDIF

         ELSE
             READ(11,*,END=120)   ! south            
         ENDIF
! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             READ(11,*,END=120)   ! north

             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
             S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)

             READ(11,119,END=120)(Z_COUPLING_NORTH(I,2),I=1,N_COUPLING_NORTH)
             READ(11,119,END=120)((U_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((V_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((W_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((P_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((S_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
             READ(11,119,END=120)((T_COUPLING_NORTH(I,J,2),I=1,N_COUPLING_NORTH),J=1,Kglob)
           IF (FirstReadCoupling)THEN
             U_COUPLING_NORTH(:,:,1)=U_COUPLING_NORTH(:,:,2)
             V_COUPLING_NORTH(:,:,1)=V_COUPLING_NORTH(:,:,2)
             W_COUPLING_NORTH(:,:,1)=W_COUPLING_NORTH(:,:,2)
             Z_COUPLING_NORTH(:,1)=Z_COUPLING_NORTH(:,2)
             P_COUPLING_NORTH(:,:,1)=P_COUPLING_NORTH(:,:,2)
             S_COUPLING_NORTH(:,:,1)=S_COUPLING_NORTH(:,:,2)
             T_COUPLING_NORTH(:,:,1)=T_COUPLING_NORTH(:,:,2)
           ENDIF

         ELSE
             READ(11,*,END=120)   ! north            
         ENDIF

         FirstReadCoupling = .FALSE.

    ENDIF  !stime>time_2 and time_1

120 CONTINUE

    tmp2=1.0
    tmp1=ZERO

    IF(TIME>TIME_COUPLING_1)THEN
      IF(TIME_COUPLING_1.EQ.TIME_COUPLING_2)THEN
        ! no more data
        tmp2=ZERO
        tmp1=1.0
      ELSE
      tmp2=(TIME_COUPLING_2-TIME) &
            /MAX(SMALL, ABS(TIME_COUPLING_2-TIME_COUPLING_1))
      tmp1=1.0_SP - tmp2;
      ENDIF  ! no more data?
    ENDIF ! time>time_1


! west boundary
   IF(N_COUPLING_WEST>0)THEN
# if defined (PARALLEL)
    if ( n_west .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_WEST)THEN

      DO J=Kstart_WEST,Kend_WEST 
      DO I=1,Nghost
        ETA(I,J)=Z_COUPLING_WEST(J-Nghost+Kshift_WEST,2)*tmp1&
                +Z_COUPLING_WEST(J-Nghost+Kshift_WEST,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost
        U(I,J,K)=U_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +U_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +V_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2  
        W(I,J,K)=W_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +W_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2  
        P(I,J,K)=P_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +P_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2 
        Sali(I,J,K)=S_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +S_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2
        Temp(I,J,K)=T_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,2)*tmp1&
                +T_COUPLING_WEST(J-Nghost+Kshift_WEST,K-Nghost,1)*tmp2

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)     
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K)     

      ENDDO
      ENDDO
      ENDDO

!print*,'intp',Sali(1,Jbeg+1,Kend),Sali(2,Jbeg+1,Kend),Sali(3,Jbeg+1,Kend),Sali(4,Jbeg+1,Kend)
!J=Kstart_WEST
!print*,'swet',S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,2),S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,1)
!print*,'calc',S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,2)*tmp1+S_COUPLING_WEST(J-Nghost+Kshift_WEST,kglob,1)*tmp2

     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_west>0


! east boundary
   IF(N_COUPLING_EAST>0)THEN
# if defined (PARALLEL)
    if ( n_east .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_EAST)THEN
      DO J=Kstart_EAST,Kend_EAST  
      DO I=Iend+1,Iend+Nghost
        ETA(I,J)=Z_COUPLING_EAST(J-Nghost+Kshift_EAST,2)*tmp1&
                +Z_COUPLING_EAST(J-Nghost+Kshift_EAST,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost
        U(I,J,K)=U_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +U_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +V_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +W_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +P_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        Sali(I,J,K)=S_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +S_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2
        Temp(I,J,K)=T_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,2)*tmp1&
                +T_COUPLING_EAST(J-Nghost+Kshift_EAST,K-Nghost,1)*tmp2

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)     
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_east>0

! south boundary
   IF(N_COUPLING_SOUTH>0)THEN
# if defined (PARALLEL)
    if ( n_suth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_SOUTH)THEN
      DO I=Kstart_SOUTH,Kend_SOUTH  
      DO J=1,Nghost
        ETA(I,J)=Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,2)*tmp1&
                +Z_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost    
        U(I,J,K)=U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +U_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +V_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +W_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +P_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        Sali(I,J,K)=S_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +S_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2
        Temp(I,J,K)=T_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,2)*tmp1&
                +T_COUPLING_SOUTH(I-Nghost+Kshift_SOUTH,K-Nghost,1)*tmp2

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)     
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_south>0

! north boundary
   IF(N_COUPLING_NORTH>0)THEN
# if defined (PARALLEL)
    if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif
     IF(IN_DOMAIN_NORTH)THEN
      DO I=Kstart_NORTH,Kend_NORTH  
      DO J=Jend+1,Jend+Nghost
        ETA(I,J)=Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,2)*tmp1&
                +Z_COUPLING_NORTH(I-Nghost+Kshift_NORTH,1)*tmp2
        D(I,J)=ETA(I,J)+Hc(I,J)
      DO K=1+Nghost,Kglob+Nghost       
        U(I,J,K)=U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +U_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        V(I,J,K)=V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +V_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        W(I,J,K)=W_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +W_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        P(I,J,K)=P_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +P_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        Sali(I,J,K)=S_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +S_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2
        Temp(I,J,K)=T_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,2)*tmp1&
                +T_COUPLING_NORTH(I-Nghost+Kshift_NORTH,K-Nghost,1)*tmp2

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)     
        Dsali(I,J,K)=D(I,J)*Sali(I,J,K)  
        Dtemp(I,J,K)=D(I,J)*Temp(I,J,K) 
      ENDDO
      ENDDO
      ENDDO
     ENDIF  ! end in domain
# if defined (PARALLEL)
    endif
# endif
    ENDIF ! end of n_coupling_north>0

END SUBROUTINE OneWayCoupling
# endif 



 
