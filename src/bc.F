# if defined (AMR)
! -----------------------------------------------
!   get bc for grid level ng from ng-1
!   for AMR only
!   called by master.F
!   ng - grid level
!   istep - time level 1 and 2
!   fyshi 08/24/2013 
! -----------------------------------------------
SUBROUTINE GET_NESTING_BC(ng,istep)
    USE GLOBAL
    USE NESTING
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ng,istep
    INTEGER :: mbox1,nbox1,isk,K,I,J

       isk=RATIO_SPACING(ng)
       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

! there's bug here. should switch if istep and if move, remember to change funwave-tvd
    IF (istep.eq.1)THEN
     IF(IS_GRID_MOVE(ng)) THEN
! west
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
         DO J=1,Nloc
          DO I=1,Nghost
            Z_NESTING_WEST(J,I,istep)=Eta_Tile(I,J,ng)
            DO K=1,Kloc
            U_NESTING_WEST(J,I,K,istep)=U_Tile(I,J,K,ng)
            V_NESTING_WEST(J,I,K,istep)=V_Tile(I,J,K,ng)
            W_NESTING_WEST(J,I,K,istep)=W_Tile(I,J,K,ng)
            P_NESTING_WEST(J,I,K,istep)=P_Tile(I,J,K,ng)
# if defined (SALINITY)
            S_NESTING_WEST(J,I,K,istep)=S_Tile(I,J,K,ng)
            T_NESTING_WEST(J,I,K,istep)=T_Tile(I,J,K,ng) 
# endif
           IF(VISCOUS_FLOW)THEN
            TKE_NESTING_WEST(J,I,K,istep)=TKE_Tile(I,J,K,ng)
            Eps_NESTING_WEST(J,I,K,istep)=Eps_Tile(I,J,K,ng)
           ENDIF
             ENDDO
          ENDDO
        ENDDO
# if defined (PARALLEL)
     endif
# endif
! east
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
        DO J=1,Nloc
          DO I=Mloc-Nghost+1,Mloc
            Z_NESTING_EAST(J,I-Mloc+Nghost,istep)=Eta_Tile(I,J,ng)
            DO K=1,Kloc
            U_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=U_Tile(I,J,K,ng)
            V_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=V_Tile(I,J,K,ng)
            W_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=W_Tile(I,J,K,ng)
            P_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=P_Tile(I,J,K,ng)
# if defined (SALINITY)
            S_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=S_Tile(I,J,K,ng)
            T_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=T_Tile(I,J,K,ng)
# endif
           IF(VISCOUS_FLOW)THEN
            TKE_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=TKE_Tile(I,J,K,ng)
            Eps_NESTING_EAST(J,I-Mloc+Nghost,K,istep)=Eps_Tile(I,J,K,ng)
           ENDIF
            ENDDO
          ENDDO
         ENDDO
# if defined (PARALLEL)
     endif
# endif

! south
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
# endif
         DO I=1,Mloc
          DO J=1,Nghost
            Z_NESTING_SOUTH(I,J,istep)=Eta_Tile(I,J,ng)
            DO K=1,Kloc
            U_NESTING_SOUTH(I,J,K,istep)=U_Tile(I,J,K,ng)
            V_NESTING_SOUTH(I,J,K,istep)=V_Tile(I,J,K,ng)
            W_NESTING_SOUTH(I,J,K,istep)=W_Tile(I,J,K,ng)
            P_NESTING_SOUTH(I,J,K,istep)=P_Tile(I,J,K,ng)
# if defined (SALINITY)
            S_NESTING_SOUTH(I,J,K,istep)=S_Tile(I,J,K,ng)
            T_NESTING_SOUTH(I,J,K,istep)=T_Tile(I,J,K,ng)
# endif
           IF(VISCOUS_FLOW)THEN
            TKE_NESTING_SOUTH(I,J,K,istep)=TKE_Tile(I,J,K,ng)
            Eps_NESTING_SOUTH(I,J,K,istep)=Eps_Tile(I,J,K,ng)
           ENDIF
            ENDDO
          ENDDO
         ENDDO
# if defined (PARALLEL)
     endif
# endif
! north
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
         DO I=1,Mloc       
          DO J=Nloc-Nghost+1,Nloc
            Z_NESTING_NORTH(I,J-Nloc+Nghost,istep)=Eta_Tile(I,J,ng)
            DO K=1,Kloc
            U_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=U_Tile(I,J,K,ng)
            V_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=V_Tile(I,J,K,ng)
            W_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=W_Tile(I,J,K,ng)
            P_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=P_Tile(I,J,K,ng)
# if defined (SALINITY)
            S_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=S_Tile(I,J,K,ng)
            T_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=T_Tile(I,J,K,ng)
# endif
           IF(VISCOUS_FLOW)THEN
            TKE_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=TKE_Tile(I,J,K,ng)
            Eps_NESTING_NORTH(I,J-Nloc+Nghost,K,istep)=Eps_Tile(I,J,K,ng)
           ENDIF
            ENDDO
          ENDDO
         ENDDO
# if defined (PARALLEL)
     endif
# endif

     ELSE ! if grid not move, interpolate using the coarse grid

         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,U_Tile(:,:,:,ng-1), &
               U_NESTING_WEST(:,:,:,istep),U_NESTING_EAST(:,:,:,istep),  &
               U_NESTING_SOUTH(:,:,:,istep),U_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,V_Tile(:,:,:,ng-1), &
               V_NESTING_WEST(:,:,:,istep),V_NESTING_EAST(:,:,:,istep),  &
               V_NESTING_SOUTH(:,:,:,istep),V_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,W_Tile(:,:,:,ng-1), &
               W_NESTING_WEST(:,:,:,istep),W_NESTING_EAST(:,:,:,istep),  &
               W_NESTING_SOUTH(:,:,:,istep),W_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,P_Tile(:,:,:,ng-1), &
               P_NESTING_WEST(:,:,:,istep),P_NESTING_EAST(:,:,:,istep),  &
               P_NESTING_SOUTH(:,:,:,istep),P_NESTING_NORTH(:,:,:,istep),ng)
# if defined (SALINITY)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,S_Tile(:,:,:,ng-1), &
               S_NESTING_WEST(:,:,:,istep),S_NESTING_EAST(:,:,:,istep),  &
               S_NESTING_SOUTH(:,:,:,istep),S_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,T_Tile(:,:,:,ng-1), &
               T_NESTING_WEST(:,:,:,istep),T_NESTING_EAST(:,:,:,istep),  &
               T_NESTING_SOUTH(:,:,:,istep),T_NESTING_NORTH(:,:,:,istep),ng)
# endif
       IF(VISCOUS_FLOW)THEN
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,TKE_Tile(:,:,:,ng-1), &
               TKE_NESTING_WEST(:,:,:,istep),TKE_NESTING_EAST(:,:,:,istep),  &
               TKE_NESTING_SOUTH(:,:,:,istep),TKE_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,Eps_Tile(:,:,:,ng-1), &
               Eps_NESTING_WEST(:,:,:,istep),Eps_NESTING_EAST(:,:,:,istep),  &
               Eps_NESTING_SOUTH(:,:,:,istep),Eps_NESTING_NORTH(:,:,:,istep),ng)
       ENDIF ! end viscous

         CALL INTERP_BC_2D(MaxDimX,MaxDimY,Mloc,Nloc,Nghost,mbox1,nbox1,isk,Eta_Tile(:,:,ng-1), &
               Z_NESTING_WEST(:,:,istep),Z_NESTING_EAST(:,:,istep),  &
               Z_NESTING_SOUTH(:,:,istep),Z_NESTING_NORTH(:,:,istep),ng)

     ENDIF ! end grid move

   ELSE ! istep=2

         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,U_Tile(:,:,:,ng-1), &
               U_NESTING_WEST(:,:,:,istep),U_NESTING_EAST(:,:,:,istep),  &
               U_NESTING_SOUTH(:,:,:,istep),U_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,V_Tile(:,:,:,ng-1), &
               V_NESTING_WEST(:,:,:,istep),V_NESTING_EAST(:,:,:,istep),  &
               V_NESTING_SOUTH(:,:,:,istep),V_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,W_Tile(:,:,:,ng-1), &
               W_NESTING_WEST(:,:,:,istep),W_NESTING_EAST(:,:,:,istep),  &
               W_NESTING_SOUTH(:,:,:,istep),W_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,P_Tile(:,:,:,ng-1), &
               P_NESTING_WEST(:,:,:,istep),P_NESTING_EAST(:,:,:,istep),  &
               P_NESTING_SOUTH(:,:,:,istep),P_NESTING_NORTH(:,:,:,istep),ng)
# if defined (SALINITY)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,S_Tile(:,:,:,ng-1), &
               S_NESTING_WEST(:,:,:,istep),S_NESTING_EAST(:,:,:,istep),  &
               S_NESTING_SOUTH(:,:,:,istep),S_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,T_Tile(:,:,:,ng-1), &
               T_NESTING_WEST(:,:,:,istep),T_NESTING_EAST(:,:,:,istep),  &
               T_NESTING_SOUTH(:,:,:,istep),T_NESTING_NORTH(:,:,:,istep),ng)
# endif
       IF(VISCOUS_FLOW)THEN
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,TKE_Tile(:,:,:,ng-1), &
               TKE_NESTING_WEST(:,:,:,istep),TKE_NESTING_EAST(:,:,:,istep),  &
               TKE_NESTING_SOUTH(:,:,:,istep),TKE_NESTING_NORTH(:,:,:,istep),ng)
         CALL INTERP_BC(MaxDimX,MaxDimY,Mloc,Nloc,Kloc,Nghost,mbox1,nbox1,isk,Eps_Tile(:,:,:,ng-1), &
               Eps_NESTING_WEST(:,:,:,istep),Eps_NESTING_EAST(:,:,:,istep),  &
               Eps_NESTING_SOUTH(:,:,:,istep),Eps_NESTING_NORTH(:,:,:,istep),ng)
        ENDIF ! end viscous

         CALL INTERP_BC_2D(MaxDimX,MaxDimY,Mloc,Nloc,Nghost,mbox1,nbox1,isk,Eta_Tile(:,:,ng-1), &
               Z_NESTING_WEST(:,:,istep),Z_NESTING_EAST(:,:,istep),  &
               Z_NESTING_SOUTH(:,:,istep),Z_NESTING_NORTH(:,:,istep),ng)

   ENDIF ! end istep


END SUBROUTINE GET_NESTING_BC
# endif
  ! end AMR

# if defined (AMR)
! -----------------------------------------------
!   3D spatial interpolation at ghost cells of grid level ng
!   for AMR only
!   called by GEG_NESTING BC
!   fyshi 08/24/2013 
! -----------------------------------------------
SUBROUTINE INTERP_BC(MaxM,MaxN,M,N,L,Nghost,mb,nb,isk,Fin,  &
                    Fout_west,Fout_east,Fout_south,Fout_north,ng)
!--------------------------------------------------------------
!  linear interpolation for 3D variables
!  fyshi 05/29/2013
!--------------------------------------------------------------
      USE PARAM
# if defined (PARALLEL)
      USE GLOBAL, ONLY : n_west,n_east,n_suth,n_nrth,MPI_PROC_NULL,&
                         Kloc,px,py,npx,npy
      USE NESTING, ONLY : GridDimX,GridDimY
# endif 
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,L,mb,nb,isk,MaxM,MaxN,Nghost,ng
      INTEGER :: II,JJ,K,I,J
      REAL(SP) :: rII,rJJ      
      REAL(SP),DIMENSION(MaxM,MaxN,L),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(MaxN,Nghost,L),INTENT(OUT) :: Fout_west,Fout_east
      REAL(SP),DIMENSION(MaxM,Nghost,L),INTENT(OUT) :: Fout_south,Fout_north
# if defined (PARALLEL)
      REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: VarGrid1  ! global including ghost
      INTEGER :: mm1,nn1,Mloc_grid1,Nloc_grid1
# endif

# if defined (PARALLEL)
      mm1 = GridDimX(ng-1)+2*Nghost
      nn1 = GridDimY(ng-1)+2*Nghost
      ALLOCATE (VarGrid1(mm1,nn1,Kloc) )

      Mloc_grid1=GridDimX(ng-1)/px+2*Nghost
      Nloc_grid1=GridDimY(ng-1)/py+2*Nghost 
      CALL GATHER_GRID_3D(VarGrid1,Fin(1:Mloc_grid1,1:Nloc_grid1,1:Kloc),&
               Mloc_grid1,Nloc_grid1,Kloc,mm1,nn1,Nghost)
# endif

! west
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
      DO K=1,L
       DO J=1,N
        DO I=1,Nghost
         II=mb+(I-1+(M-2*Nghost)*npx)/isk  ! actually npx=0 for west
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_west(J,I,K)=((1.0_SP-rII)*VarGrid1(II,JJ,K)+rII*VarGrid1(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1,K)+rII*VarGrid1(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO
      ENDDO
     endif

# else
    ! not parallel
      DO K=1,L
       DO J=1,N
        DO I=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_west(J,I,K)=((1.0_SP-rII)*Fin(II,JJ,K)+rII*Fin(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1,K)+rII*Fin(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO
      ENDDO
# endif 
! end parallel

! east
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
      DO K=1,L           
       DO J=1,N  
        DO I=M-Nghost+1,M
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_east(J,I-M+Nghost,K)=((1.0_SP-rII)*VarGrid1(II,JJ,K)+rII*VarGrid1(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1,K)+rII*VarGrid1(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO 
      ENDDO
     endif
# else
    ! not parallel
      DO K=1,L           
       DO J=1,N  
        DO I=M-Nghost+1,M
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_east(J,I-M+Nghost,K)=((1.0_SP-rII)*Fin(II,JJ,K)+rII*Fin(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1,K)+rII*Fin(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO 
      ENDDO
# endif 
! end parallel

! south
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
      DO K=1,L
       DO I=1,M
        DO J=1,Nghost
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk  ! actually npy=0 at south
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_south(I,J,K)=((1.0_SP-rII)*VarGrid1(II,JJ,K)+rII*VarGrid1(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1,K)+rII*VarGrid1(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO
      ENDDO
     endif
# else
    ! not parallel
      DO K=1,L
       DO I=1,M
        DO J=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_south(I,J,K)=((1.0_SP-rII)*Fin(II,JJ,K)+rII*Fin(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1,K)+rII*Fin(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO
      ENDDO
# endif 
! end parallel

! north
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
      DO K=1,L
       DO I=1,M
        DO J=N-Nghost+1,N
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_north(I,J-N+Nghost,K)=((1.0_SP-rII)*VarGrid1(II,JJ,K)+rII*VarGrid1(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1,K)+rII*VarGrid1(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO 
      ENDDO
     endif
# else
    ! not parallel
      DO K=1,L
       DO I=1,M
        DO J=N-Nghost+1,N
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_north(I,J-N+Nghost,K)=((1.0_SP-rII)*Fin(II,JJ,K)+rII*Fin(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1,K)+rII*Fin(II+1,JJ+1,K))*rJJ
        ENDDO
       ENDDO 
      ENDDO
# endif 
! end parallel

# if defined (PARALLEL)
      DEALLOCATE (VarGrid1)
# endif

END SUBROUTINE INTERP_BC
# endif
  ! end AMR

# if defined (AMR)
! -----------------------------------------------
!   2D spatial interpolation at ghost cells of grid level ng
!   for AMR only
!   called by GEG_NESTING BC
!   fyshi 08/24/2013 
! -----------------------------------------------
SUBROUTINE INTERP_BC_2D(MaxM,MaxN,M,N,Nghost,mb,nb,isk,Fin,  &
                    Fout_west,Fout_east,Fout_south,Fout_north,ng)
!--------------------------------------------------------------
!  linear interpolation for 2D variables
!  fyshi 05/29/2013
!--------------------------------------------------------------
      USE PARAM
# if defined (PARALLEL)
      USE GLOBAL, ONLY : n_west,n_east,n_suth,n_nrth,MPI_PROC_NULL,&
                         px,py,npx,npy
      USE NESTING, ONLY : GridDimX,GridDimY
# endif 
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: M,N,mb,nb,isk,MaxM,MaxN,Nghost,ng
      INTEGER :: II,JJ,I,J
      REAL(SP) :: rII,rJJ      
      REAL(SP),DIMENSION(MaxM,MaxN),INTENT(IN) :: Fin
      REAL(SP),DIMENSION(MaxN,Nghost),INTENT(OUT) :: Fout_west,Fout_east
      REAL(SP),DIMENSION(MaxM,Nghost),INTENT(OUT) :: Fout_south,Fout_north
# if defined (PARALLEL)
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGrid1  ! global including ghost
      INTEGER :: mm1,nn1,Mloc_grid1,Nloc_grid1
# endif

# if defined (PARALLEL)
      mm1 = GridDimX(ng-1)+2*Nghost
      nn1 = GridDimY(ng-1)+2*Nghost
      ALLOCATE (VarGrid1(mm1,nn1) )

      Mloc_grid1=GridDimX(ng-1)/px+2*Nghost
      Nloc_grid1=GridDimY(ng-1)/py+2*Nghost 
      CALL GATHER_GRID_2D(VarGrid1,Fin(1:Mloc_grid1,1:Nloc_grid1),&
               Mloc_grid1,Nloc_grid1,mm1,nn1,Nghost)
# endif

! west 
# if defined (PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
       DO J=1,N
        DO I=1,Nghost
         II=mb+(I-1+(M-2*Nghost)*npx)/isk  ! actually npx=0 at west
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_west(J,I)=((1.0_SP-rII)*VarGrid1(II,JJ)+rII*VarGrid1(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1)+rII*VarGrid1(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
     endif
# else
! not parallel
       DO J=1,N
        DO I=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_west(J,I)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO
# endif

! east 
# if defined (PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
      DO J=1,N
        DO I=M-Nghost+1,M
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_east(J,I-M+Nghost)=((1.0_SP-rII)*VarGrid1(II,JJ)+rII*VarGrid1(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1)+rII*VarGrid1(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
     endif
# else
! not parallel
       DO J=1,N
        DO I=M-Nghost+1,M
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_east(J,I-M+Nghost)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
# endif

! south
# if defined (PARALLEL)
     if(n_suth.eq.MPI_PROC_NULL) then
       DO I=1,M
        DO J=1,Nghost
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk  ! actually npy=0 at south
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_south(I,J)=((1.0_SP-rII)*VarGrid1(II,JJ)+rII*VarGrid1(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1)+rII*VarGrid1(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO
     endif
# else
! not parallel
       DO I=1,M
        DO J=1,Nghost
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_south(I,J)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO
# endif

! north
# if defined (PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
       DO I=1,M
        DO J=N-Nghost+1,N
         II=mb+(I-1+(M-2*Nghost)*npx)/isk
         JJ=nb+(J-1+(N-2*Nghost)*npy)/isk
         rII=REAL(mb)+REAL(I-1+(M-2*Nghost)*npx)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1+(N-2*Nghost)*npy)/REAL(isk)-REAL(JJ)
         Fout_north(I,J-N+Nghost)=((1.0_SP-rII)*VarGrid1(II,JJ)+rII*VarGrid1(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*VarGrid1(II,JJ+1)+rII*VarGrid1(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
     endif
# else
! not parallel
       DO I=1,M
        DO J=N-Nghost+1,N
         II=mb+(I-1)/isk
         JJ=nb+(J-1)/isk
         rII=REAL(mb)+REAL(I-1)/REAL(isk)-REAL(II)
         rJJ=REAL(nb)+REAL(J-1)/REAL(isk)-REAL(JJ)
         Fout_north(I,J-N+Nghost)=((1.0_SP-rII)*Fin(II,JJ)+rII*Fin(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Fin(II,JJ+1)+rII*Fin(II+1,JJ+1))*rJJ
        ENDDO
       ENDDO 
# endif

END SUBROUTINE INTERP_BC_2D
# endif
  ! end AMR

# if defined (AMR)
! -----------------------------------------------
!   apply nesting bc
!   for AMR only
!   called by main.F
!   istep - time level
!   ng - grid level
!   ratio - useless so far
!   fyshi last update 08/24/2013 
! -----------------------------------------------
SUBROUTINE USE_NESTING_BC(istep,ratio,ng)
    USE GLOBAL
    USE NESTING
    IMPLICIT NONE
    REAL(SP) :: tmp1,tmp2
    INTEGER :: istep,ng
    INTEGER :: ratio,I,J,K
     real(SP), parameter :: Dmin = 0.04,cmiu = 0.09
     ! this should be consistent with in turb.F

! time
      tmp1=REAL(istep)/REAL(ratio)    
      tmp2=1.0_SP - tmp1;

! west boundary
# if defined(PARALLEL)
     if(n_west.eq.MPI_PROC_NULL) then
# endif
      DO J=1,Nloc
      DO I=1,Nghost
        ETA(I,J)=Z_NESTING_WEST(J,I,2)*tmp1&
                +Z_NESTING_WEST(J,I,1)*tmp2
        D(I,J)=Hc(I,J)+ETA(I,J)
        DO K=1,Kloc
        U(I,J,K)=U_NESTING_WEST(J,I,K,2)*tmp1&
                +U_NESTING_WEST(J,I,K,1)*tmp2
        V(I,J,K)=V_NESTING_WEST(J,I,K,2)*tmp1&
                +V_NESTING_WEST(J,I,K,1)*tmp2
        W(I,J,K)=W_NESTING_WEST(J,I,K,2)*tmp1&
                +W_NESTING_WEST(J,I,K,1)*tmp2
        P(I,J,K)=P_NESTING_WEST(J,I,K,2)*tmp1&
                +P_NESTING_WEST(J,I,K,1)*tmp2
# if defined (SALINITY)
        Sali(I,J,K)=S_NESTING_WEST(J,I,K,2)*tmp1&
                +S_NESTING_WEST(J,I,K,1)*tmp2
        Temp(I,J,K)=T_NESTING_WEST(J,I,K,2)*tmp1&
                +T_NESTING_WEST(J,I,K,1)*tmp2
# endif

      IF(VISCOUS_FLOW)THEN
        Tke(I,J,K)=TKE_NESTING_WEST(J,I,K,2)*tmp1&
                +TKE_NESTING_WEST(J,I,K,1)*tmp2
        Eps(I,J,K)=Eps_NESTING_WEST(J,I,K,2)*tmp1&
                +Eps_NESTING_WEST(J,I,K,1)*tmp2
        DTKE(I,J,K)=D(I,J)*Tke(I,J,K)
        DEps(I,J,K)=D(I,J)*Eps(I,J,K)
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
        CmuHt(i,j,k) = CmuVt(i,j,k)

      ELSE ! not viscous
        CmuVt(i,j,k) = Cmut_min
        CmuHt(i,j,k) = CmuVt(i,j,k)
      ENDIF ! end viscous

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)
# if defined (SALINITY)
        DSali(I,J,K)=D(I,J)*Sali(I,J,K)
        DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif
        ENDDO
      ENDDO
      ENDDO
# if defined (PARALLEL)
     ENDIF  ! end west
# endif

! east boundary
# if defined(PARALLEL)
     if(n_east.eq.MPI_PROC_NULL) then
# endif
      DO J=1,Nloc 
      DO I=Iend+1,Iend+Nghost
        ETA(I,J)=Z_NESTING_EAST(J,I-Iend,2)*tmp1&
                +Z_NESTING_EAST(J,I-Iend,1)*tmp2
        D(I,J)=Hc(I,J)+ETA(I,J)
        DO K=1,Kloc
        U(I,J,K)=U_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +U_NESTING_EAST(J,I-Iend,K,1)*tmp2
        V(I,J,K)=V_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +V_NESTING_EAST(J,I-Iend,K,1)*tmp2
        W(I,J,K)=W_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +W_NESTING_EAST(J,I-Iend,K,1)*tmp2
        P(I,J,K)=P_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +P_NESTING_EAST(J,I-Iend,K,1)*tmp2
# if defined (SALINITY)
        Sali(I,J,K)=S_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +S_NESTING_EAST(J,I-Iend,K,1)*tmp2
        Temp(I,J,K)=T_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +T_NESTING_EAST(J,I-Iend,K,1)*tmp2
# endif

      IF(VISCOUS_FLOW)THEN
        TKE(I,J,K)=TKE_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +TKE_NESTING_EAST(J,I-Iend,K,1)*tmp2
        Eps(I,J,K)=Eps_NESTING_EAST(J,I-Iend,K,2)*tmp1&
                +Eps_NESTING_EAST(J,I-Iend,K,1)*tmp2
        DTKE(I,J,K)=D(I,J)*Tke(I,J,K)
        DEps(I,J,K)=D(I,J)*Eps(I,J,K)
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
        CmuHt(i,j,k) = CmuVt(i,j,k)

      ELSE ! not viscous
        CmuVt(i,j,k) = Cmut_min
        CmuHt(i,j,k) = CmuVt(i,j,k)
      ENDIF ! end viscous

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)

# if defined (SALINITY)
        DSali(I,J,K)=D(I,J)*Sali(I,J,K)
        DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif

        ENDDO
      ENDDO
      ENDDO
# if defined (PARALLEL)
     ENDIF  ! end east
# endif


     IF(.NOT.NESTING_X1D) THEN
! south boundary
# if defined(PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
      DO I=1,Mloc
      DO J=1,Nghost
        ETA(I,J)=Z_NESTING_SOUTH(I,J,2)*tmp1&
                +Z_NESTING_SOUTH(I,J,1)*tmp2
        D(I,J)=Hc(I,J)+ETA(I,J)
        DO K=1,Kloc
        U(I,J,K)=U_NESTING_SOUTH(I,J,K,2)*tmp1&
                +U_NESTING_SOUTH(I,J,K,1)*tmp2
        V(I,J,K)=V_NESTING_SOUTH(I,J,K,2)*tmp1&
                +V_NESTING_SOUTH(I,J,K,1)*tmp2
        W(I,J,K)=W_NESTING_SOUTH(I,J,K,2)*tmp1&
                +W_NESTING_SOUTH(I,J,K,1)*tmp2
        P(I,J,K)=P_NESTING_SOUTH(I,J,K,2)*tmp1&
                +P_NESTING_SOUTH(I,J,K,1)*tmp2
# if defined (SALINITY)
        Sali(I,J,K)=S_NESTING_SOUTH(I,J,K,2)*tmp1&
                +S_NESTING_SOUTH(I,J,K,1)*tmp2
        Temp(I,J,K)=T_NESTING_SOUTH(I,J,K,2)*tmp1&
                +T_NESTING_SOUTH(I,J,K,1)*tmp2
# endif

      IF(VISCOUS_FLOW)THEN
        TKE(I,J,K)=TKE_NESTING_SOUTH(I,J,K,2)*tmp1&
                +TKE_NESTING_SOUTH(I,J,K,1)*tmp2
        Eps(I,J,K)=Eps_NESTING_SOUTH(I,J,K,2)*tmp1&
                +Eps_NESTING_SOUTH(I,J,K,1)*tmp2
        DTKE(I,J,K)=D(I,J)*Tke(I,J,K)
        DEps(I,J,K)=D(I,J)*Eps(I,J,K)
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
        CmuHt(i,j,k) = CmuVt(i,j,k)

      ELSE ! no viscous
        CmuVt(i,j,k) = Cmut_min
        CmuHt(i,j,k) = CmuVt(i,j,k)
      ENDIF ! end viscous

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)
# if defined (SALINITY)
        DSali(I,J,K)=D(I,J)*Sali(I,J,K)
        DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif

        ENDDO
      ENDDO
      ENDDO 
# if defined (PARALLEL)
     ENDIF  ! end suth
# endif

! north boundary
# if defined(PARALLEL)
     if(n_nrth.eq.MPI_PROC_NULL) then
# endif
      DO I=1,Mloc
      DO J=Jend+1,Jend+Nghost
        ETA(I,J)=Z_NESTING_NORTH(I,J-Jend,2)*tmp1&
                +Z_NESTING_NORTH(I,J-Jend,1)*tmp2
        D(I,J)=Hc(I,J)+ETA(I,J)
        DO K=1,Kloc
        U(I,J,K)=U_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +U_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        V(I,J,K)=V_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +V_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        W(I,J,K)=W_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +W_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        P(I,J,K)=P_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +P_NESTING_NORTH(I,J-Jend,K,1)*tmp2
# if defined (SALINITY)
        Sali(I,J,K)=S_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +S_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        Temp(I,J,K)=T_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +T_NESTING_NORTH(I,J-Jend,K,1)*tmp2
# endif
      IF(VISCOUS_FLOW)THEN
        TKE(I,J,K)=TKE_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +TKE_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        Eps(I,J,K)=Eps_NESTING_NORTH(I,J-Jend,K,2)*tmp1&
                +Eps_NESTING_NORTH(I,J-Jend,K,1)*tmp2
        DTKE(I,J,K)=D(I,J)*Tke(I,J,K)
        DEps(I,J,K)=D(I,J)*Eps(I,J,K)
      if(D(i,j)>=Dmin.and.Mask(i,j)==1) then
        CmuVt(i,j,k) = Cmiu*Tke(i,j,k)**2/MAX(Eps(i,j,k),Eps_min)
      else
        CmuVt(i,j,k) = Cmut_min
      endif
        CmuHt(i,j,k) = CmuVt(i,j,k)

      ELSE ! no viscous
        CmuVt(i,j,k) = Cmut_min
        CmuHt(i,j,k) = CmuVt(i,j,k)
      ENDIF ! end viscous

        DU(I,J,K)=D(I,J)*U(I,J,K)
        DV(I,J,K)=D(I,J)*V(I,J,K)
        DW(I,J,K)=D(I,J)*W(I,J,K)
# if defined (SALINITY)
        DSali(I,J,K)=D(I,J)*Sali(I,J,K)
        DTemp(I,J,K)=D(I,J)*Temp(I,J,K)
# endif


        ENDDO
      ENDDO
      ENDDO 
# if defined (PARALLEL)
     ENDIF  ! end nrth
# endif

   ENDIF ! end nesting_x1d

END SUBROUTINE USE_NESTING_BC 
# endif
  ! end AMR



