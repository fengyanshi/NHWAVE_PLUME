        PROGRAM EXTRACT_DATA
        USE GLOBAL

        USE INPUT_UTIL
        IMPLICIT NONE
        INTEGER :: I,J,K,line,Nstart,Nend,Nlayer,kk,ierr,itmp1,itmp2,itmp3,itmp4
        INTEGER :: XYsection,XZsection,YZsection
        LOGICAL :: Sali_OUT = .FALSE.
        LOGICAL :: u_OUT = .FALSE.
        LOGICAL :: v_OUT = .FALSE.
        LOGICAL :: w_OUT = .FALSE.
        LOGICAL :: mu_OUT = .FALSE.
        LOGICAL :: rho_OUT = .FALSE.
        LOGICAL :: p_OUT = .FALSE.

        CHARACTER(LEN=80) File_Name,Out_Name
        INTEGER :: NumData,nwest,nsouth,neast,nnorth

        REAL,DIMENSION(:),ALLOCATABLE :: TIME_DATA
        REAL,DIMENSION(:,:),ALLOCATABLE :: philoc
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: mu
        REAL :: TRamp_DATA = 3600.0
        REAL :: Ramp
        REAL :: T_m2,omega_m2,phase_lag,U_inlet,U_sea,A_inlet,A_sea,&
                S_inlet,S_sea,T_inlet,T_sea,V_river,S_river,A_river

! get data from input.txt

     ! read from input.txt
       FILE_NAME='input.txt'

     ! dimension                                             
       CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
       CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
       CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line) 
       CALL GET_INTEGER_VAL(Nstart,FILE_NAME,'Nstart',line) 
       CALL GET_INTEGER_VAL(Nend,FILE_NAME,'Nend',line) 
       CALL GET_INTEGER_VAL(XYsection,FILE_NAME,'XYsection',line)
       CALL GET_INTEGER_VAL(XZsection,FILE_NAME,'XZsection',line)
       CALL GET_INTEGER_VAL(YZsection,FILE_NAME,'YZsection',line)
       CALL GET_LOGICAL_VAL(Sali_OUT,FILE_NAME,'Sali_OUT',line)
       CALL GET_LOGICAL_VAL(u_OUT,FILE_NAME,'u_OUT',line)
       CALL GET_LOGICAL_VAL(v_OUT,FILE_NAME,'v_OUT',line)
       CALL GET_LOGICAL_VAL(w_OUT,FILE_NAME,'w_OUT',line)
       CALL GET_LOGICAL_VAL(mu_OUT,FILE_NAME,'mu_OUT',line)
       CALL GET_LOGICAL_VAL(rho_OUT,FILE_NAME,'rho_OUT',line)
       CALL GET_LOGICAL_VAL(p_OUT,FILE_NAME,'p_OUT',line)
       CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)

print*,'XYsection: ',XYsection
print*,'XZsection: ',XZsection
print*,'YZsection: ',YZsection


       ALLOCATE(Sali(Mglob,Nglob,Kglob),u(Mglob,Nglob,Kglob),&
                v(Mglob,Nglob,Kglob),w(Mglob,Nglob,Kglob),mu(Mglob,Nglob,Kglob), &
                rho(Mglob,Nglob,Kglob),p(Mglob,Nglob,Kglob))
       ALLOCATE(philoc(Mglob,Nglob*Kglob))

        DO kk=Nstart,Nend

        icount=kk
        itmp1=mod(icount/1000,10)
        itmp2=mod(icount/100,10)
        itmp3=mod(icount/10,10)
        itmp4=mod(icount,10)

        write(out_name(1:1),'(I1)')itmp1
        write(out_name(2:2),'(I1)')itmp2
        write(out_name(3:3),'(I1)')itmp3
        write(out_name(4:4),'(I1)')itmp4

print*,'extract file number: ',TRIM(OUT_NAME)

! SALI

     IF (Sali_OUT) THEN
print*,'READ/WRITE Sali ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'sali_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (sali(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
print*,'XY section: ', XYsection
        open(2,file=TRIM(RESULT_FOLDER)//'saliXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (sali(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
print*,'XZ section: ', XZsection
        open(2,file=TRIM(RESULT_FOLDER)//'saliXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (sali(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
print*,'YZ section: ', YZsection
        open(2,file=TRIM(RESULT_FOLDER)//'saliYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (sali(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF

! U
     IF (u_OUT) THEN
print*,'READ/WRITE u ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'u_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (u(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'uXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (u(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'uXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (u(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'uYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (u(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF

! V 

     IF (v_OUT) THEN
print*,'READ/WRITE v ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'v_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (v(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'vXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (v(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'vXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (v(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'vYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (v(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF

! W

     IF (w_OUT) THEN
print*,'READ/WRITE w ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'w_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (w(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'wXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (w(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'wXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (w(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'wYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (w(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF

! mu

     IF (mu_OUT) THEN
print*,'READ/WRITE mu ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'mu_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (mu(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'muXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (mu(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'muXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (mu(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'muYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (mu(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF

! rho

     IF (rho_OUT) THEN
print*,'READ/WRITE rho ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'rho_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (rho(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'rhoXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (rho(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'rhoXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (rho(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'rhoYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (rho(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF
 
! P

     IF (p_OUT) THEN
print*,'READ/WRITE p ...... '
        open(1,file=TRIM(RESULT_FOLDER)//'p_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (p(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

      IF(XYsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'pXY_'//TRIM(OUT_NAME))
          do k = XYsection,XYsection
          do j = 1,Nglob
            WRITE(2,100) (p(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(XZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'pXZ_'//TRIM(OUT_NAME))
          do j = XZsection,XZsection
          do k = 1,Kglob
            WRITE(2,100) (p(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)
      ENDIF 
      IF(YZsection .GT. 0) THEN
        open(2,file=TRIM(RESULT_FOLDER)//'pYZ_'//TRIM(OUT_NAME))
          do i = YZsection,YZsection
          do k = 1,Kglob
            WRITE(2,100) (p(i,j,k),j=1,Nglob)
          enddo
          enddo    
        close(2)
      ENDIF 
     ENDIF


100 FORMAT(15000E16.8)

        ENDDO  ! end file number

        END









