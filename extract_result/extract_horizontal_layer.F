        PROGRAM EXTRACT_DATA
        USE GLOBAL

        USE INPUT_UTIL
        IMPLICIT NONE
        INTEGER :: I,J,K,line,Nstart,Nend,Nlayer,kk,ierr,itmp1,itmp2,itmp3,itmp4
        CHARACTER(LEN=80) File_Name,Out_Name
        INTEGER :: NumData,nwest,nsouth,neast,nnorth

        REAL,DIMENSION(:),ALLOCATABLE :: TIME_DATA
        REAL,DIMENSION(:,:),ALLOCATABLE :: philoc
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
       CALL GET_INTEGER_VAL(Nlayer,FILE_NAME,'Nlayer',line)  
       CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)

       ALLOCATE(Sali(Mglob,Nglob,Kglob),u(Mglob,Nglob,Kglob),&
                v(Mglob,Nglob,Kglob),w(Mglob,Nglob,Kglob))
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

        open(1,file=TRIM(RESULT_FOLDER)//'sali_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (sali(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

        open(2,file=TRIM(RESULT_FOLDER)//'sali2D_'//TRIM(OUT_NAME))
          do k = Nlayer,Nlayer
          do j = 1,Nglob
            WRITE(2,100) (sali(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)

! uuuuu
        open(1,file=TRIM(RESULT_FOLDER)//'u_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (u(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

        open(2,file=TRIM(RESULT_FOLDER)//'u2D_'//TRIM(OUT_NAME))
          do k = Nlayer,Nlayer
          do j = 1,Nglob
            WRITE(2,100) (u(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)

! vvvvv
        open(1,file=TRIM(RESULT_FOLDER)//'v_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (v(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

        open(2,file=TRIM(RESULT_FOLDER)//'v2D_'//TRIM(OUT_NAME))
          do k = Nlayer,Nlayer
          do j = 1,Nglob
            WRITE(2,100) (v(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)

! wwww
        open(1,file=TRIM(RESULT_FOLDER)//'w_'//TRIM(OUT_NAME))
          do k = 1,Kglob
          do j = 1,Nglob
            READ(1,*) (w(i,j,k),i=1,Mglob)
          enddo
          enddo        
        close(1)

        open(2,file=TRIM(RESULT_FOLDER)//'w2D_'//TRIM(OUT_NAME))
          do k = Nlayer,Nlayer
          do j = 1,Nglob
            WRITE(2,100) (w(i,j,k),i=1,Mglob)
          enddo
          enddo    
        close(2)


100 FORMAT(15000E16.8)

        ENDDO  ! end file number

        END









