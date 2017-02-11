# if defined (AMR)
! ---------------------------------------------------------
!      This is the module handling nesting of AMR
!      called by MASTER and MAIN
!      last update 08/24/2013, fyshi
!
! ---------------------------------------------------------
MODULE NESTING
        USE PARAM
        IMPLICIT NONE
! NESTING grid numbers

         SAVE

! grid
         REAL(SP),DIMENSION(:),ALLOCATABLE :: DX_TILE,DY_TILE

! bathy
         REAL(SP),DIMENSION(:,:,:), ALLOCATABLE :: Depth_Tile, &
                 DepthX_Tile, DepthY_Tile
         INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: MASK_Tile,MASK9_Tile

! variable
         REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: Eta_Tile
         REAL(SP),DIMENSION(:,:,:,:),ALLOCATABLE :: U_Tile,V_Tile,W_Tile,&
                  P_Tile,Omega_Tile
# if defined (SALINITY)
         REAL(SP),DIMENSION(:,:,:,:),ALLOCATABLE ::S_Tile,T_Tile
# endif
         REAL(SP),DIMENSION(:,:,:,:),ALLOCATABLE ::Rho_Tile,Tke_Tile,Eps_Tile 

       INTEGER :: NumGrid =1
       INTEGER :: START_COUNT
       INTEGER :: MaxDimX,MaxDimY,MaxDimX1,MaxDimY1,MaxDimZ,MaxDimZ1
       INTEGER,DIMENSION(:),ALLOCATABLE :: GridDimX,GridDimY,&
             RATIO_SPACING

       LOGICAL,DIMENSION(:),ALLOCATABLE :: IS_GRID_MOVE
       LOGICAL,DIMENSION(:),ALLOCATABLE  :: IsPrint
       INTEGER,DIMENSION(:),ALLOCATABLE :: MboxRef,NboxRef,&
                                           MboxRef0,NboxRef0, &
                                           MboxRefFix,NboxRefFix ! fix only used in fix grid
 
       INTEGER,DIMENSION(:),ALLOCATABLE :: &
                  ng_Kstart_EAST,ng_Kend_EAST,ng_Kshift_EAST,&
                  ng_Kstart_WEST,ng_Kend_WEST,ng_Kshift_WEST,&
                  ng_Kstart_SOUTH,ng_Kend_SOUTH,ng_Kshift_SOUTH,&
                  ng_Kstart_NORTH,ng_Kend_NORTH,ng_Kshift_NORTH,&
                  Istart_overlap,Iend_overlap,&
                  Jstart_overlap,Jend_overlap, &
                  Ishift_overlap,Jshift_overlap
                      
       LOGICAL,DIMENSION(:),ALLOCATABLE :: &  
                  ng_IN_DOMAIN,I_IN_DOMAIN_overlap,J_IN_DOMAIN_overlap
       LOGICAL,DIMENSION(:),ALLOCATABLE :: START_GRID  
       LOGICAL :: NESTING_X1D = .FALSE.
       CHARACTER(LEN=80) :: DEPTH_FULL_FILE
       INTEGER :: Mdep,Ndep,MainGridRatio
       REAL,DIMENSION(:,:),ALLOCATABLE ::Depth_Full


       CONTAINS

! ---------------------------------------------------------
!      control grid move
!      called by MASTER 
!      last update 08/24/2013, fyshi
!
! ---------------------------------------------------------   
SUBROUTINE MOVING_GRID(ng)  
!---------------------------------------------- 
! to get (MboxRef,NboxRef), 
!  Istart_overlap,Iend_overlap,Jstart_overlap,Jend_overlap
!     fyshi 05/29/2013
!-----------------------------------------------
     USE GLOBAL
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng
     INTEGER :: mbox1,nbox1,isk,m_move,n_move

       MboxRef0(ng)=MboxRef(ng)
       NboxRef0(ng)=NboxRef(ng)
       MboxRef(ng)=MboxRefFix(ng)
       NboxRef(ng)=NboxRefFix(ng)

! mbox is box in a coarser grid, include ghostcells, m represents m+nghost
       isk=RATIO_SPACING(ng)

! adjust if the box moves out of domain
       mbox1=MIN(GridDimX(ng-1)-GridDimX(ng)/isk+1,MboxRef(ng))
       nbox1=MIN(GridDimY(ng-1)-GridDimY(ng)/isk+1,NboxRef(ng))
       MboxRef(ng)=mbox1
       NboxRef(ng)=nbox1


! calculate move m and n in finer domain
       m_move=(MboxRef(ng)-MboxRef0(ng))*RATIO_SPACING(ng)  ! move number in ng
       n_move=(NboxRef(ng)-NboxRef0(ng))*RATIO_SPACING(ng)

! to see if box in last step in domain
!       IF(1+n_move.LT.Nloc.AND.1+m_move.LT.Mloc  &
!          .AND.Nloc+n_move.GT.1.AND.Mloc+m_move.GT.1)THEN
!         ng_IN_DOMAIN(ng) = .TRUE.
!       ENDIF

! X direction
      Istart_overlap(ng)=MAX(1,1+m_move)
      Iend_overlap(ng) = MIN(Mloc,Mloc+m_move)
      Ishift_overlap(ng) = 0
      I_IN_DOMAIN_overlap(ng)=.TRUE.      

! Y direction
      Jstart_overlap(ng)=MAX(1,1+n_move)
      Jend_overlap(ng) = MIN(Nloc,Nloc+n_move)
      Jshift_overlap(ng) = 0
      J_IN_DOMAIN_overlap(ng)=.TRUE.

END SUBROUTINE MOVING_GRID

! ---------------------------------------------------------
!      initialize grid after moving
!      called by MASTER 
!      last update 08/24/2013, fyshi
!
! ---------------------------------------------------------
SUBROUTINE INITIAL_GRID(ng)
! --------------------------------------------
!  To initialize var_tile for ng grid
!   fyshi 05/29/2013
!---------------------------------------------
     USE GLOBAL    
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: ng
     INTEGER :: mbox1,nbox1,isk,mbox0,nbox0,m_move,n_move
     INTEGER :: mbox_to_grid1,nbox_to_grid1
     REAL(SP) :: ratio_to_grid1,ratio_to_grid2
     INTEGER :: I,J,K,II,JJ
     REAL(SP) :: rII,rJJ

! for ng=1
     IF(ng==1)THEN
      DX_Tile(ng)=DX
      DY_Tile(ng)=DY
      DO J=1,Nloc
      DO I=1,Mloc
       Depth_Tile(I,J,ng)=Hc(I,J)
       Eta_Tile(I,J,ng)=Eta(I,J)
       MASK_Tile(I,J,ng)=MASK(I,J)
       MASK9_Tile(I,J,ng)=MASK9(I,J)
      ENDDO
      ENDDO
      DO J=1,Nloc
      DO I=1,Mloc1
       DepthX_Tile(I,J,ng)=Hfx(I,J)
      ENDDO
      ENDDO
      DO J=1,Nloc1
      DO I=1,Mloc
       DepthY_Tile(I,J,ng)=Hfy(I,J)
      ENDDO
      ENDDO

      DO K=1,Kloc
      DO J=1,Nloc
      DO I=1,Mloc
       U_Tile(I,J,K,ng)=U(I,J,K)
       V_Tile(I,J,K,ng)=V(I,J,K)
       W_Tile(I,J,K,ng)=W(I,J,K)
       P_Tile(I,J,K,ng)=P(I,J,K)
# if defined (SALINITY)
       S_Tile(I,J,K,ng)=Sali(I,J,K)
       T_Tile(I,J,K,ng)=Temp(I,J,K)
# endif

     IF(VISCOUS_FLOW)THEN
       Rho_Tile(I,J,K,ng)=Rho(I,J,K)
       Tke_Tile(I,J,K,ng)=Tke(I,J,K)
       Eps_Tile(I,J,K,ng)=Eps(I,J,K)
     ENDIF

      ENDDO
      ENDDO
      ENDDO

      DO K=1,Kloc1
      DO J=1,Nloc
      DO I=1,Mloc
        Omega_Tile(I,J,K,ng)=Omega(I,J,K)
      ENDDO
      ENDDO
      ENDDO

     ELSE
! for ng>1
! mbox is box in a coarser grid, include ghostcells, m represents m+nghost

     IF(IS_GRID_MOVE(ng))THEN

       isk=RATIO_SPACING(ng)

! mbox1 and nbox1 respect to ng-1 grid
       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

       mbox_to_grid1=0
       nbox_to_grid1=0
       DO K=2,ng
         ratio_to_grid1=1.0_SP
         ratio_to_grid2=1.0_SP   ! last level grid
         DO J=2,K
           ratio_to_grid1=ratio_to_grid1*RATIO_SPACING(J)
           ratio_to_grid2=ratio_to_grid2*RATIO_SPACING(J-1)           
         ENDDO
         mbox_to_grid1=mbox_to_grid1+INT(MboxRef(K)/ratio_to_grid2)
         nbox_to_grid1=nbox_to_grid1+INT(NboxRef(K)/ratio_to_grid2)
       ENDDO

!print*,mbox_to_grid1,nbox_to_grid1,ratio_to_grid1,ratio_to_grid2
!stop

! static
       DX_Tile(ng)=DX_Tile(ng-1)/isk
       DY_Tile(ng)=DY_Tile(ng-1)/isk

!       CALL UPDATE_FINER_GRID(Depth_Tile)

  
! depth
# if defined (USE_INTERP)

      CALL INITIALIZE_GRID2_FROM_GRID1_2D(Depth_Tile(:,:,ng-1),Depth_Tile(:,:,ng),ng)

# else 
! use finest resolution given by depth_full
     DO J=1,Nloc
     DO I=1,Mloc
# if defined (PARALLEL)
       II=INT(REAL((mbox_to_grid1-Nghost-1)*MainGridRatio)+ &
          REAL(I)*REAL(MainGridRatio)/REAL(ratio_to_grid1)+ &
          GridDimX(ng)/px*npx*REAL(MainGridRatio)/REAL(ratio_to_grid1))

       JJ=INT(REAL((nbox_to_grid1-Nghost-1)*MainGridRatio)+ &
          REAL(J)*REAL(MainGridRatio)/REAL(ratio_to_grid1)+ &
          GridDimY(ng)/py*npy*REAL(MainGridRatio)/REAL(ratio_to_grid1))

# else
       II=(mbox_to_grid1-Nghost-1)*MainGridRatio+INT(REAL(I)*REAL(MainGridRatio)/REAL(ratio_to_grid1))
       JJ=(nbox_to_grid1-Nghost-1)*MainGridRatio+INT(REAL(J)*REAL(MainGridRatio)/REAL(ratio_to_grid1))
# endif
       IF(NESTING_X1D)THEN
         JJ=1
       ENDIF
       Depth_Tile(I,J,ng)=Depth_Full(II,JJ)
     ENDDO
     ENDDO

# endif 
! end use_interp

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       DepthX_Tile(i,j,ng) = 0.5*(Depth_Tile(i-1,j,ng)+Depth_Tile(i,j,ng))
     enddo
     DepthX_Tile(1,j,ng) = Depth_Tile(1,j,ng)
     DepthX_Tile(Mloc1,j,ng) = Depth_Tile(Mloc,j,ng)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       DepthY_Tile(i,j,ng) = 0.5*(Depth_Tile(i,j-1,ng)+Depth_Tile(i,j,ng))
     enddo
     DepthY_Tile(i,1,ng) = Depth_Tile(i,1,ng)
     DepthY_Tile(i,Nloc1,ng) = Depth_Tile(i,Nloc,ng)
     enddo


! Eta


      CALL INITIALIZE_GRID2_FROM_GRID1_2D(Eta_Tile(:,:,ng-1),Eta_Tile(:,:,ng),ng)


! calculate MASK and MASK9 separately 
       DO J=1,Nloc
       DO I=1,Mloc
        IF(Eta_Tile(I,J,ng)<-DEPTH_Tile(I,J,ng))THEN
         MASK_Tile(I,J,ng)=0
         Eta_Tile(I,J,ng)=MinDep-Depth_Tile(I,J,ng)
        ELSE
         MASK_Tile(I,J,ng)=1
        ENDIF
       ENDDO
       ENDDO

     DO J=Jbeg-1,Jend+1
     DO I=Ibeg-1,Iend+1
      MASK9_Tile(I,J,ng)=MASK_Tile(I,J,ng)*MASK_Tile(I-1,J,ng)*MASK_Tile(I+1,J,ng)  &
                *MASK_Tile(I+1,J+1,ng)*MASK_Tile(I,J+1,ng)*MASK_Tile(I-1,J+1,ng) &
                *MASK_Tile(I+1,J-1,ng)*MASK_Tile(I,J-1,ng)*MASK_Tile(I-1,J-1,ng) 
     ENDDO
     ENDDO

! U
   
       CALL INITIALIZE_GRID2_FROM_GRID1_3D(U_Tile(:,:,:,ng-1),U_Tile(:,:,:,ng),ng)


! V

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(V_Tile(:,:,:,ng-1),V_Tile(:,:,:,ng),ng)

! W

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(W_Tile(:,:,:,ng-1),W_Tile(:,:,:,ng),ng)

! Omega

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(Omega_Tile(:,:,:,ng-1),Omega_Tile(:,:,:,ng),ng)


! P

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(P_Tile(:,:,:,ng-1),P_Tile(:,:,:,ng),ng)

# if defined (SALINITY)
! S

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(S_Tile(:,:,:,ng-1),S_Tile(:,:,:,ng),ng)

! T

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(T_Tile(:,:,:,ng-1),T_Tile(:,:,:,ng),ng)

# endif

! Rho

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(Rho_Tile(:,:,:,ng-1),Rho_Tile(:,:,:,ng),ng)


       IF(VISCOUS_FLOW)THEN

! Tke

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(Tke_Tile(:,:,:,ng-1),Tke_Tile(:,:,:,ng),ng)

! Eps

       CALL INITIALIZE_GRID2_FROM_GRID1_3D(Eps_Tile(:,:,:,ng-1),Eps_Tile(:,:,:,ng),ng)

      ENDIF ! end viscous_flow

! exchange for parallel !!! 

     ELSE ! grid not move
        ! only boundary updating, see get_nesting_bc in bc.F
     ENDIF  ! end grid move

     ENDIF  ! end ng 

END SUBROUTINE INITIAL_GRID

     
! ---------------------------------------------------------
!      allocate variables in nesting module
!      called by MASTER 
!      last update 08/24/2013, fyshi
!
! ---------------------------------------------------------
SUBROUTINE ALLOCATE_VAR_TILE
       USE PARAM
# if defined(PARALLEL)
       USE GLOBAL, ONLY: Nghost,ier,px,py,Kglob    
# else
       USE GLOBAL, ONLY: Nghost,Kglob
# endif

       IMPLICIT NONE

# if defined(PARALLEL)
       MaxDimX = MAXVAL(GridDimX)/px+2*Nghost
       MaxDimY = MAXVAL(GridDimY)/py+2*Nghost
# else
       MaxDimX = MAXVAL(GridDimX)+2*Nghost
       MaxDimY = MAXVAL(GridDimY)+2*Nghost
# endif
       MaxDimX1 = MaxDimX+1
       MaxDimY1 = MaxDimY+1

       MaxDimZ = Kglob +2*Nghost
       MaxDimZ1 = MaxDimZ+1
       

! static
       ALLOCATE ( DX_TILE(NumGrid))
       ALLOCATE ( DY_TILE(NumGrid))
       ALLOCATE ( Depth_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( DepthX_Tile(MaxDimX1,MaxDimY,NumGrid))
       ALLOCATE ( DepthY_Tile(MaxDimX,MaxDimY1,NumGrid))
       ALLOCATE ( MASK_Tile(MaxDimX,MaxDimY,NumGrid))
       ALLOCATE ( MASK9_Tile(MaxDimX,MaxDimY,NumGrid))

! dynamic
       ALLOCATE ( U_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( V_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( W_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( Omega_Tile(MaxDimX,MaxDimY,MaxDimZ1,NumGrid))
       ALLOCATE ( P_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
# if defined (SALINITY)
       ALLOCATE ( S_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( T_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
# endif

       ALLOCATE ( Rho_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( Tke_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))
       ALLOCATE ( Eps_Tile(MaxDimX,MaxDimY,MaxDimZ,NumGrid))

       ALLOCATE ( Eta_Tile(MaxDimX,MaxDimY,NumGrid)) 

! initialization
        MASK_TILE = 1
        MASK9_TILE = 1
        U_TILE = ZERO
        V_TILE = ZERO
        W_TILE = ZERO
        P_TILE = ZERO
# if defined(SALINITY)
        S_TILE = ZERO
        T_TILE = ZERO
# endif
        Eta_TILE = ZERO
        Rho_TILE = ZERO
        Tke_TILE = ZERO
        Eps_TILE = ZERO

END SUBROUTINE ALLOCATE_VAR_TILE


END MODULE NESTING
# endif
  ! end AMR