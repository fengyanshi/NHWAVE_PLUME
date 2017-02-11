! ----------------------------------------
!  This is nesting related subroutine but independent
!  from module nesting. 
!   called by INITIAL_GRID
!  update 08/24/2013, fyshi
! ----------------------------------------
# if defined (AMR)
SUBROUTINE INITIALIZE_GRID2_FROM_GRID1_3D(V_grid1,V_Grid2,ng)
      USE PARAM
      USE GLOBAL,ONLY : Nghost,Mloc,Nloc,Kloc
      USE NESTING, ONLY : MaxDimX,MaxDimY,GridDimX,GridDimY,RATIO_SPACING,&
                          MboxRef,NboxRef,Istart_overlap,Iend_overlap,&
                          Jstart_overlap,Jend_overlap
# if defined (PARALLEL)
      USE GLOBAL,ONLY : npx,npy,px,py
# endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ng
      REAL(SP),DIMENSION(MaxDimX,MaxDimY,Kloc),INTENT(IN)::V_grid1
      REAL(SP),DIMENSION(MaxDimX,MaxDimY,Kloc),INTENT(INOUT)::V_grid2
      REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: VarGrid1,VarGrid2  ! global including ghost

      INTEGER :: Mloc_grid1,Nloc_grid1,mm1,nn1,mm2,nn2,I,J,II,JJ,K
      INTEGER :: mbox1,nbox1,isk,mbox0,nbox0,m_move,n_move     

      mm2 = GridDimX(ng)+2*Nghost
      nn2 = GridDimY(ng)+2*Nghost
      ALLOCATE (VarGrid2(mm2,nn2,Kloc) )
      mm1 = GridDimX(ng-1)+2*Nghost
      nn1 = GridDimY(ng-1)+2*Nghost
      ALLOCATE (VarGrid1(mm1,nn1,Kloc) )

      isk=RATIO_SPACING(ng)

       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

# if defined (PARALLEL)
      Mloc_grid1=GridDimX(ng-1)/px+2*Nghost
      Nloc_grid1=GridDimY(ng-1)/py+2*Nghost 
# endif

# if defined (PARALLEL)
      CALL GATHER_GRID_3D(VarGrid1,V_grid1(1:Mloc_grid1,1:Nloc_grid1,1:Kloc),&
               Mloc_grid1,Nloc_grid1,Kloc,mm1,nn1,Nghost)
      CALL GATHER_GRID_3D(VarGrid2,V_grid2(1:Mloc,1:Nloc,1:Kloc), &
               Mloc,Nloc,Kloc,mm2,nn2,Nghost)
# else
      VarGrid1(1:mm1,1:nn1,1:Kloc) = V_grid1(1:mm1,1:nn1,1:Kloc) 
      VarGrid2(1:mm2,1:nn2,1:Kloc) = V_grid2(1:mm2,1:nn2,1:Kloc)
# endif

      CALL  Grid1_Grid2_3D (mm1,nn1,mm2,nn2, Kloc, &
            MboxRef(ng),NboxRef(ng),RATIO_SPACING(ng), &
            Istart_overlap(ng),Iend_overlap(ng),Jstart_overlap(ng),Jend_overlap(ng),&
            VarGrid1,VarGrid2)

# if defined (PARALLEL)
        DO J=1,Nloc
        DO I=1,Mloc
        DO K=1,Kloc
          II=I+npx*(Mloc-2*Nghost)
          JJ=J+npy*(Nloc-2*Nghost)
          V_grid2(I,J,K)=VarGrid2(II,JJ,K)
        ENDDO
        ENDDO
        ENDDO
# else
        V_grid2(1:Mloc,1:Nloc,1:Kloc)=VarGrid2(1:Mloc,1:Nloc,1:Kloc)
# endif


      DEALLOCATE (VarGrid2)
      DEALLOCATE (VarGrid1)

END SUBROUTINE INITIALIZE_GRID2_FROM_GRID1_3D
# endif
  ! end AMR

! --------------------------------------------
!  Gather 3D variables from all processors
! this gathering includes ghost cells
! 08/20/2013, fyshi
! --------------------------------------------
# if defined (PARALLEL)
SUBROUTINE GATHER_GRID_3D(phi_out,phi_in,&
        Mloc,Nloc,Kloc,mm,nn,Nghost)

! mm and nn are global but include ghost cells

      USE PARAM
      USE GLOBAL, ONLY : NumP,npx,npy,myid
      IMPLICIT NONE
    integer,intent(in) :: Mloc,Nloc,mm,nn,Nghost,Kloc
    integer :: ier
    real(SP),dimension(Mloc,Nloc,Kloc),intent(in) :: phi_in
    real(SP),dimension(mm,nn,Kloc),intent(out) :: phi_out
    integer,dimension(NumP) :: npxs,npys
    real(SP),dimension(NumP) :: xxx
    integer,dimension(1) :: req
    real(SP),dimension(:,:),allocatable :: xx,philoc
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,k,jk,iglob,jglob,kk,n,len,nreq,NKloc,l

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    NKloc = Nloc*Kloc

    ! put the data in master processor into the global var
    if(myid==0) then
! doesn't include ghost cells, OK for interpolation
      do k = 1,Kloc
      do j = 1,Nloc
      do i = 1,Mloc
        iglob = i
        jglob = j
        kk=k
        phi_out(iglob,jglob,kk) = Phi_in(i,j,k)
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
      philoc(i,jk) = phi_in(i,j,k)
    enddo
    enddo
    enddo

    ! collect data from other processors into the master processor
    len = Mloc*NKloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do k = 1,Kloc
        do j = 1,Nloc
        do i = 1,Mloc
          iglob = npxs(n+1)*(Mloc-2*Nghost)+i
          jglob = npys(n+1)*(Nloc-2*Nghost)+j
          kk=k
          jk = (k-1)*Nloc+j
          phi_out(iglob,jglob,kk) = xx(i,jk)
        enddo
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi_in,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo     

! scattering to every processors
        do k = 1,Kloc
        do j = 1,nn
        do i = 1,mm
        if (myid.eq.0) then
           do l=1,NumP
              xxx(l) = phi_out(i,j,k)
           enddo
        endif
        call MPI_Scatter(xxx,1,MPI_SP,&
             phi_out(i,j,k),1,MPI_SP,0,MPI_COMM_WORLD,ier)
        enddo
        enddo
        enddo

!    if(myid==0) then
!      open(5,file='tmp.txt')
!      do j = 1,Nglob
!        write(5,100) (phi_out(i,j),i=1,Mglob)
!      enddo
!      close(5)
!    endif
!100 FORMAT(5000f15.5)

    deallocate(philoc)
    deallocate(xx)

END SUBROUTINE GATHER_GRID_3D
# endif

! --------------------------------------------
!  interpolation for 3D variables from grid1 to grid2
! 
! 08/20/2013, fyshi
! --------------------------------------------
SUBROUTINE Grid1_Grid2_3D(m1,n1,m2,n2,Kloc,mbox1,nbox1,isk,&
           Ist,Ien,Jst,Jen,Var1,Var2)

       USE PARAM
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: m2,n2,m1,n1,mbox1,nbox1,isk,&
           Ist,Ien,Jst,Jen,Kloc
       REAL(SP),DIMENSION(m1,n1,Kloc),INTENT(IN) ::  Var1
       REAL(SP),DIMENSION(m2,n2,Kloc),INTENT(INOUT) :: Var2
       INTEGER :: I,J,K,II,JJ,rII,rJJ
       REAL(SP),DIMENSION(m2,n2,Kloc) :: tmp3D

       tmp3D=Var2


       DO J=1,n2
       DO I=1,m2
       DO K=1,Kloc
        II=mbox1+(I-1)/isk
        JJ=nbox1+(J-1)/isk
        rII=REAL(mbox1)+REAL(I-1)/REAL(isk)-REAL(II)
        rJJ=REAL(nbox1)+REAL(J-1)/REAL(isk)-REAL(JJ)
        Var2(I,J,K)=((1.0_SP-rII)*Var1(II,JJ,K)+rII*Var1(II+1,JJ,K))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Var1(II,JJ+1,K)+rII*Var1(II+1,JJ+1,K))*rJJ
       ENDDO
       ENDDO
       ENDDO 


         DO J=Jst,Jen 
           DO I=Ist,Ien
             DO K=1,Kloc
               Var2(I-Ist+1,J-Jst+1,K)=tmp3D(I,J,K)
             ENDDO
           ENDDO
         ENDDO

END SUBROUTINE Grid1_Grid2_3D

# if defined (AMR)
! ----------------------------------------
!  This is nesting related subroutine but independent
!  from module nesting. 
!   called by INITIAL_GRID
!  update 08/24/2013, fyshi
! ----------------------------------------
SUBROUTINE INITIALIZE_GRID2_FROM_GRID1_2D(V_grid1,V_Grid2,ng)
      USE PARAM
      USE GLOBAL,ONLY : Nghost,Mloc,Nloc
      USE NESTING, ONLY : MaxDimX,MaxDimY,GridDimX,GridDimY,RATIO_SPACING,&
                          MboxRef,NboxRef,Istart_overlap,Iend_overlap,&
                          Jstart_overlap,Jend_overlap
# if defined (PARALLEL)
      USE GLOBAL,ONLY : npx,npy,px,py,myid
# endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ng
      REAL(SP),DIMENSION(MaxDimX,MaxDimY),INTENT(IN)::V_grid1
      REAL(SP),DIMENSION(MaxDimX,MaxDimY),INTENT(INOUT)::V_grid2
      REAL(SP),DIMENSION(:,:),ALLOCATABLE :: VarGrid1,VarGrid2  ! global including ghost

      INTEGER :: Mloc_grid1,Nloc_grid1,mm1,nn1,mm2,nn2,I,J,II,JJ,ier
      INTEGER :: mbox1,nbox1,isk,mbox0,nbox0,m_move,n_move     

      mm2 = GridDimX(ng)+2*Nghost
      nn2 = GridDimY(ng)+2*Nghost
      ALLOCATE (VarGrid2(mm2,nn2) )
      mm1 = GridDimX(ng-1)+2*Nghost
      nn1 = GridDimY(ng-1)+2*Nghost
      ALLOCATE (VarGrid1(mm1,nn1) )

      isk=RATIO_SPACING(ng)

       mbox1=MboxRef(ng)
       nbox1=NboxRef(ng)

# if defined (PARALLEL)
      Mloc_grid1=GridDimX(ng-1)/px+2*Nghost
      Nloc_grid1=GridDimY(ng-1)/py+2*Nghost 
# endif


# if defined (PARALLEL)
      CALL GATHER_GRID_2D(VarGrid1,V_grid1(1:Mloc_grid1,1:Nloc_grid1),&
               Mloc_grid1,Nloc_grid1,mm1,nn1,Nghost)
      CALL GATHER_GRID_2D(VarGrid2,V_grid2(1:Mloc,1:Nloc),Mloc,Nloc,mm2,nn2,Nghost)

# else
      VarGrid1(1:mm1,1:nn1) = V_grid1(1:mm1,1:nn1) 
      VarGrid2(1:mm2,1:nn2) = V_grid2(1:mm2,1:nn2)
# endif

      CALL  Grid1_Grid2_2D (mm1,nn1,mm2,nn2, &
            MboxRef(ng),NboxRef(ng),RATIO_SPACING(ng), &
            Istart_overlap(ng),Iend_overlap(ng),Jstart_overlap(ng),Jend_overlap(ng),&
            VarGrid1,VarGrid2)


# if defined (PARALLEL)
        DO J=1,Nloc
        DO I=1,Mloc
          II=I+npx*(Mloc-2*Nghost)
          JJ=J+npy*(Nloc-2*Nghost)
          V_grid2(I,J)=VarGrid2(II,JJ)
        ENDDO
        ENDDO
# else
        V_grid2(1:Mloc,1:Nloc)=VarGrid2(1:Mloc,1:Nloc)
# endif

      DEALLOCATE (VarGrid2)
      DEALLOCATE (VarGrid1)

END SUBROUTINE INITIALIZE_GRID2_FROM_GRID1_2D
# endif
  ! end AMR

# if defined (PARALLEL)
! --------------------------------------------
!  Gather 2D variables from all processors
! this gathering includes ghost cells
! 08/20/2013, fyshi
! --------------------------------------------
SUBROUTINE GATHER_GRID_2D(phi_out,phi_in,&
        Mloc,Nloc,mm,nn,Nghost)

! mm and nn are global but include ghost cells

      USE PARAM
      USE GLOBAL, ONLY : NumP,npx,npy,myid
      IMPLICIT NONE
    integer,intent(in) :: Mloc,Nloc,mm,nn,Nghost
    integer :: ier
    real(SP),dimension(Mloc,Nloc),intent(in) :: phi_in
    real(SP),dimension(mm,nn),intent(out) :: phi_out
    integer,dimension(NumP) :: npxs,npys
    integer,dimension(1) :: req
    real(SP),dimension(Mloc,Nloc) :: xx
    real(SP),dimension(NumP) :: xxx
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    integer :: i,j,iglob,jglob,len,n,l

    call MPI_GATHER(npx,1,MPI_INTEGER,npxs,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)
    call MPI_GATHER(npy,1,MPI_INTEGER,npys,1,MPI_INTEGER,  &
           0,MPI_COMM_WORLD,ier)

    ! put the data in master processor into the global var
    if(myid==0) then
!  include ghost cells
      do j = 1,Nloc
      do i = 1,Mloc
        iglob = i
        jglob = j
        phi_out(iglob,jglob) = Phi_in(i,j)
      enddo
      enddo
    endif

    ! collect data from other processors into the master processor
    len = Mloc*Nloc

    do n = 1,NumP-1
      if(myid==0) then
        call MPI_IRECV(xx,len,MPI_SP,n,0,MPI_COMM_WORLD,req(1),ier)
        call MPI_WAITALL(1,req,status,ier)
        do j = 1,Nloc
        do i = 1,Mloc
          iglob = npxs(n+1)*(Mloc-2*Nghost)+i
          jglob = npys(n+1)*(Nloc-2*Nghost)+j
          phi_out(iglob,jglob) = xx(i,j)
        enddo
        enddo
      endif

      if(myid==n) then
        call MPI_SEND(phi_in,len,MPI_SP,0,0,MPI_COMM_WORLD,ier)
      endif
    enddo   

! scattering to every processors
        do j = 1,nn
        do i = 1,mm
        if (myid.eq.0) then
           do l=1,NumP
              xxx(l) = phi_out(i,j)
           enddo
        endif
        call MPI_Scatter(xxx,1,MPI_SP,&
             phi_out(i,j),1,MPI_SP,0,MPI_COMM_WORLD,ier)
        enddo
        enddo

!    if(myid==0) then
!      open(5,file='tmp.txt')
!      do j = 1,Nglob
!        write(5,100) (phi_out(i,j),i=1,Mglob)
!      enddo
!      close(5)
!    endif
!100 FORMAT(5000f15.5)



END SUBROUTINE GATHER_GRID_2D
# endif

! --------------------------------------------
!  interpolation for 2D variables from grid1 to grid2
! 
! 08/20/2013, fyshi
! --------------------------------------------
SUBROUTINE Grid1_Grid2_2D(m1,n1,m2,n2,mbox1,nbox1,isk,&
           Ist,Ien,Jst,Jen,Var1,Var2)

       USE PARAM
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: m2,n2,m1,n1,mbox1,nbox1,isk,&
           Ist,Ien,Jst,Jen
       REAL(SP),DIMENSION(m1,n1),INTENT(IN) ::  Var1
       REAL(SP),DIMENSION(m2,n2),INTENT(INOUT) :: Var2
       INTEGER :: I,J,K,II,JJ,rII,rJJ
       REAL(SP),DIMENSION(m2,n2) :: tmp2D

       tmp2D=Var2

       DO J=1,n2
       DO I=1,m2
        II=mbox1+(I-1)/isk
        JJ=nbox1+(J-1)/isk
        rII=REAL(mbox1)+REAL(I-1)/REAL(isk)-REAL(II)
        rJJ=REAL(nbox1)+REAL(J-1)/REAL(isk)-REAL(JJ)
        Var2(I,J)=((1.0_SP-rII)*Var1(II,JJ)+rII*Var1(II+1,JJ))*(1.0_SP-rJJ)+ &
                        ((1.0_SP-rII)*Var1(II,JJ+1)+rII*Var1(II+1,JJ+1))*rJJ
       ENDDO
       ENDDO 


         DO J=Jst,Jen 
           DO I=Ist,Ien
             Var2(I-Ist+1,J-Jst+1)=tmp2D(I,J)
           ENDDO
         ENDDO


END SUBROUTINE Grid1_Grid2_2D


