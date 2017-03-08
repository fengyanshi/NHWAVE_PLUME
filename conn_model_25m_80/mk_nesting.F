        PROGRAM MAKE_NESTING_FILES
        USE GLOBAL

        USE INPUT_UTIL
        IMPLICIT NONE
        INTEGER :: I,J,K,line
        CHARACTER(LEN=80) File_Name
        INTEGER :: NumData,nwest,nsouth,neast,nnorth

        REAL,DIMENSION(:),ALLOCATABLE :: TIME_DATA
        REAL :: TRamp_DATA = 3600.0
        REAL :: Ramp
        REAL :: T_m2,omega_m2,phase_lag,U_inlet,U_sea,A_inlet,A_sea,&
                S_inlet,S_sea,T_inlet,T_sea,V_river,S_river,A_river

! input data
         NumData = 1.0*4.0   ! 1 hours 
! get data from input.txt

     ! read from input.txt
       FILE_NAME='input.txt'

     ! dimension                                             
       CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
       CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
       CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)   

        N_COUPLING_EAST = Nglob
        N_COUPLING_WEST = Nglob
        N_COUPLING_SOUTH =0 
        N_COUPLING_NORTH =200       ! change 25

        J_START_WEST = 1
        J_START_EAST = 1
        I_START_SOUTH = 1        
        I_START_NORTH = 100     ! change 75


        OPEN(11,FILE='coupling.txt')
! Allocate

         ALLOCATE (TIME_DATA(NumData))
         
           DO K=1,NumData
             TIME_DATA(K)=(K-1.)*900.0
           ENDDO


         IF(N_COUPLING_EAST.GT.0)THEN
           ALLOCATE(U_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               V_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               W_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData),&
               Z_COUPLING_EAST(N_COUPLING_EAST,NumData), &
               P_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData), &
               S_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData), &
               T_COUPLING_EAST(N_COUPLING_EAST,Kglob,NumData) &
               )

           T_m2=12.42*3600.0
           omega_m2=2.0*3.1415926/T_m2
           phase_lag=3.1415926/180.0*45.0
           U_inlet = 1.0                     ! change inlet is on right
           U_sea = 1.0                       ! change
           A_inlet=-0.1                      ! change
           A_sea=0.1                         ! change
           S_inlet=25.0                      ! change
           T_inlet=15.0                      ! change
           S_sea=25.0                        ! change
           T_sea=15.0                        ! change
           V_river = -0.4                    ! change
           S_river = 5.0                     ! change
           A_river = 0.075                     ! change


           DO K=1,NumData
           DO I=1,N_COUPLING_EAST
             Z_COUPLING_EAST(I,K)=A_inlet !*SIN(omega_m2*TIME_DATA(K)+phase_lag)
           DO J=1,Kglob
             U_COUPLING_EAST(I,J,K)=U_inlet !*SIN(omega_m2*TIME_DATA(K)+phase_lag)
             V_COUPLING_EAST(I,J,K)=0.0
             W_COUPLING_EAST(I,J,K)=0.0
             P_COUPLING_EAST(I,J,K)=0.0
             S_COUPLING_EAST(I,J,K)=S_inlet
             T_COUPLING_EAST(I,J,K)=T_inlet
           ENDDO
           ENDDO
           ENDDO
         ENDIF

         IF(N_COUPLING_WEST.GT.0)THEN
           ALLOCATE(U_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               V_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               W_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData),&
               Z_COUPLING_WEST(N_COUPLING_WEST,NumData), &
               P_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData), &
               S_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData), &
               T_COUPLING_WEST(N_COUPLING_WEST,Kglob,NumData) &
               )
           DO K=1,NumData
           DO I=1,N_COUPLING_WEST
             Z_COUPLING_WEST(I,K)=A_sea !*SIN(omega_m2*TIME_DATA(K)+phase_lag)
           DO J=1,Kglob
             U_COUPLING_WEST(I,J,K)=U_sea !*SIN(omega_m2*TIME_DATA(K)+phase_lag)
             V_COUPLING_WEST(I,J,K)=0.0
             W_COUPLING_WEST(I,J,K)=0.0
             P_COUPLING_WEST(I,J,K)=0.0
             S_COUPLING_WEST(I,J,K)=S_sea
             T_COUPLING_WEST(I,J,K)=T_sea
           ENDDO
           ENDDO
           ENDDO
         ENDIF

         IF(N_COUPLING_SOUTH.GT.0)THEN
           ALLOCATE(U_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               V_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               W_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData),&
               Z_COUPLING_SOUTH(N_COUPLING_SOUTH,NumData), &
               P_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData), &
               S_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData), &
               T_COUPLING_SOUTH(N_COUPLING_SOUTH,Kglob,NumData) &
               )
           DO K=1,NumData
           DO I=1,N_COUPLING_SOUTH
             Z_COUPLING_SOUTH(I,K)=0.5
           DO J=1,Kglob
             U_COUPLING_SOUTH(I,J,K)=0.0
             V_COUPLING_SOUTH(I,J,K)=1.0
             W_COUPLING_SOUTH(I,J,K)=0.0
             P_COUPLING_SOUTH(I,J,K)=0.0
             S_COUPLING_SOUTH(I,J,K)=0.0
             T_COUPLING_SOUTH(I,J,K)=0.0
           ENDDO
           ENDDO
           ENDDO
         ENDIF

         IF(N_COUPLING_NORTH.GT.0)THEN
           ALLOCATE(U_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               V_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               W_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData),&
               Z_COUPLING_NORTH(N_COUPLING_NORTH,NumData), &
               P_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData), &
               S_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData), &
               T_COUPLING_NORTH(N_COUPLING_NORTH,Kglob,NumData) &
               )
           DO K=1,NumData
           DO I=1,N_COUPLING_NORTH
             Z_COUPLING_NORTH(I,K)=A_river
           DO J=1,Kglob
             U_COUPLING_NORTH(I,J,K)=0.0
             V_COUPLING_NORTH(I,J,K)=V_river
             W_COUPLING_NORTH(I,J,K)=0.0
             P_COUPLING_NORTH(I,J,K)=0.0
             S_COUPLING_NORTH(I,J,K)=S_river
             T_COUPLING_NORTH(I,J,K)=T_sea
           ENDDO
           ENDDO
           ENDDO
         ENDIF


         WRITE(11,*)  'coupling data'
         WRITE(11,*)  'boundary info'
! boundary basic info including point number of coupling, start point, etc
! east
         WRITE(11,*)  'N_COUPLING_EAST,J_START_EAST'
         WRITE(11,'(2I8)') N_COUPLING_EAST,J_START_EAST
! west 
         WRITE(11,*)  'N_COUPLING_WEST,J_START_WEST'
         WRITE(11,'(2I8)') N_COUPLING_WEST,J_START_WEST
! south 
         WRITE(11,*)  'N_COUPLING_SOUTH,I_START_SOUTH'
         WRITE(11,'(2I8)') N_COUPLING_SOUTH,I_START_SOUTH
! north 
         WRITE(11,*)  'N_COUPLING_NORTH,I_START_NORTH'
         WRITE(11,'(2I8)') N_COUPLING_NORTH,I_START_NORTH

! read time and variable 
119      FORMAT(5E16.6)

       DO K=1,NumData

         Ramp=TANH(TIME_DATA(K)/TRAMP_DATA)

         WRITE(11,*) 'TIME_COUPLING'
         WRITE(11,'(F16.8)') TIME_DATA(K)
      
! east
         IF(N_COUPLING_EAST.GT.0)THEN
             WRITE(11,*)   'east'
             WRITE(11,119)(Z_COUPLING_EAST(I,K),I=1,N_COUPLING_EAST)
             WRITE(11,119)((U_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((S_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
             WRITE(11,119)((T_COUPLING_EAST(I,J,K),I=J_START_EAST,N_COUPLING_EAST+J_START_EAST-1),J=1,Kglob)
         ELSE
             WRITE(11,*) 'east'
         ENDIF

! west
         IF(N_COUPLING_WEST.GT.0)THEN
             WRITE(11,*)   'west'
             WRITE(11,119)(Z_COUPLING_WEST(I,K),I=1,N_COUPLING_WEST)
             WRITE(11,119)((U_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((S_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
             WRITE(11,119)((T_COUPLING_WEST(I,J,K),I=1,N_COUPLING_WEST),J=1,Kglob)
         ELSE
             WRITE(11,*) 'west'
         ENDIF          

! south
         IF(N_COUPLING_SOUTH.GT.0)THEN
             WRITE(11,*)   'south'
             WRITE(11,119)(Z_COUPLING_SOUTH(I,K),I=1,N_COUPLING_SOUTH)
             WRITE(11,119)((U_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((S_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
             WRITE(11,119)((T_COUPLING_SOUTH(I,J,K),I=1,N_COUPLING_SOUTH),J=1,Kglob)
         ELSE
             WRITE(11,*) 'south'
         ENDIF

! north
         IF(N_COUPLING_NORTH.GT.0)THEN
             WRITE(11,*)   'north'
             WRITE(11,119)(Z_COUPLING_NORTH(I,K),I=1,N_COUPLING_NORTH)
             WRITE(11,119)((U_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((V_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((W_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((P_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((S_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
             WRITE(11,119)((T_COUPLING_NORTH(I,J,K),I=1,N_COUPLING_NORTH),J=1,Kglob)
         ELSE
             WRITE(11,*) 'north'
         ENDIF

        ENDDO ! end numdata

        END









