     SUBROUTINE INITIALIZATION
!---------------------------------------------------  
!    This subroutine is used to initialize model run 
!    Called by                                       
!       main 
!    update: 21/12/2010, Gangfeng Ma  
!    update 05/27/2013 fyshi for AMR version
!    this subroutine is only for grid(1)
!    other initializations are done in initial_grid
!---------------------------------------------------
     use global
     implicit none
     integer  :: i,j,k,iter,n,m,nmax,Itmp,Jtmp,Ktmp,&
                 IItmp,JJtmp,KKtmp
     real(SP) :: xsol(80),zsol(80),zmax,xmax,xterp,zterp,tmp,zc(Kloc), &
                 utmp1,wtmp1,utmp2,wtmp2,xk(321,16),zk(321,16),  &
                 uk(321,16),wk(321,16)
     real(SP) :: alpha0,L0,T,bl,wl,e,kb,kw,x0,x1,x2,xr,  &
                 Zslide(Mloc,Nloc),xt,xt1,zt,zt1,Zlev,Slope,Xslide
! fyshi add tmp 3d variable for reading 
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE :: TMP_READ
     INTEGER :: Iglob, Jglob

#if defined (OBSTACLE)
     INTEGER ::  Number_mask
     REAL(SP) :: U_read,V_read,W_read,Ax_read,&
                 Ay_read,Az_read,Surf_bottom_read, &
                 Surf_top_read
#endif

     ALLOCATE (TMP_READ(Mglob,Nglob,Kglob))


     ! simulation time
     TIME = Zero
     RUN_STEP = 0
     dt_grid1 = dt_ini   
     dt=dt_ini
     Screen_Count = Zero
     Plot_Count = Zero
     Plot_Count_Stat = Zero
# if defined (SEDIMENT)
     Update_Bed_T = Zero
# endif
     Icount = 0
     
     ! working arrays
     D = Zero
     U = Zero
     V = Zero
     W = Zero
     P = Zero
     Omega = Zero
     DU = Zero
     DV = Zero
     DW = Zero
     D0 = Zero
     Eta0 = Zero
     DU0 = Zero
     DV0 = Zero
     DW0 = Zero
     Uf = Zero
     Vf = Zero
     Wf = Zero
     Rho = Rho0
     
     ! source terms
     SourceC = Zero
     SourceX = Zero
     SourceY = Zero

     ! fluxes
     DxL = Zero
     DxR = Zero
     DyL = Zero
     DyR = Zero
     UxL = Zero
     UxR = Zero
     UyL = Zero
     UyR = Zero
     UzL = Zero
     UzR = Zero
     VxL = Zero
     VxR = Zero
     VyL = Zero
     VyR = Zero
     VzL = Zero
     VzR = Zero
     WxL = Zero
     WxR = Zero
     WyL = Zero
     WyR = Zero
     WzL = Zero
     WzR = Zero
     DUxL = Zero
     DUxR = Zero
     DUyL = Zero
     DUyR = Zero
     DVxL = Zero
     DVxR = Zero
     DVyL = Zero
     DVyR = Zero
     DWxL = Zero
     DWxR = Zero
     DWyL = Zero
     DWyR = Zero
     OzL = Zero
     OzR = Zero
     SxL = Zero
     SxR = Zero
     SxS = Zero
     SyL = Zero
     SyR = Zero
     SyS = Zero
     Ex = Zero
     Ey = Zero
     Ez = Zero
     Fx = Zero
     Fy = Zero
     Fz = Zero
     Gx = Zero
     Gy = Zero
     Gz = Zero
     Hx = Zero
     Hy = Zero
     Hz = Zero
     EtaxL = Zero
     EtaxR = Zero
     EtayL = Zero
     EtayR = Zero
     DelxEta = Zero
     DelyEta = Zero
     DeltH = Zero
     DeltHo = Zero
     Delt2H = Zero
     DelxD = Zero
     DelyD = Zero
     DelxU = Zero
     DelyU = Zero
     DelzU = Zero
     DelxV = Zero
     DelyV = Zero
     DelzV = Zero
     DelxW = Zero
     DelyW = Zero
     DelzW = Zero
     DelxDU = Zero
     DelyDU = Zero
     DelxDV = Zero
     DelyDV = Zero
     DelzO = Zero
     Sponge = One
     Cmu = Visc
     CmuHt = Zero
     CmuVt = Zero
     CmuR = ZERO
     Richf = ZERO

     Diffxx = Zero
     Diffxy = Zero
     Diffxz = Zero
     Diffyx = Zero
     Diffyy = Zero
     Diffyz = Zero
     Diffzx = Zero
     Diffzy = Zero
     Diffzz = Zero

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     Ein_X0 = Zero
     Din_X0 = Zero
     Uin_Xn = Zero
     Vin_Xn = Zero
     Win_Xn = Zero
     Ein_Xn = Zero
     Din_Xn = Zero  

     Uin_Y0 = Zero
     Vin_Y0 = Zero
     Win_Y0 = Zero
     Ein_Y0 = Zero
     Din_Y0 = Zero
     Uin_Yn = Zero
     Vin_Yn = Zero
     Win_Yn = Zero
     Ein_Yn = Zero
     Din_Yn = Zero  

     Setup = Zero
     WaveHeight = Zero
     Umean = Zero
     Vmean = Zero
     Num_Zero_Up = 0
     Emax = -1000.
     Emin = 1000.

     Lag_Umean = Zero
     Lag_Vmean = Zero
     Lag_Wmean = Zero
     Euler_Umean = Zero
     Euler_Vmean = Zero
     Euler_Wmean = Zero

     ExtForceX = Zero
     ExtForceY = Zero

     ! baroclinic terms
     DRhoX = Zero
     DRhoY = Zero


# if defined (BUBBLE)
     Nbg = Zero
     DNbg = Zero
     DNbg0 = Zero
     Vbg = Zero
# endif

# if defined (SEDIMENT)
     Conc = Zero
     DConc = Zero
     DConc0 = Zero
     Bed = Zero
     Sed_X0 = Zero
     Sed_Xn = Zero
# endif

# if defined(AIR_PRESSURE)
     Pmb = ZERO
# endif

# if defined (SALINITY)
    IF(INI_SALI_INPUT(1:4)=='CONS')THEN
     Sali = INI_SALI

     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       DSali(I,J,K)=Hc(I,J)*Sali(I,J,K)
     ENDDO
     ENDDO
     ENDDO

    ELSE IF(INI_SALI_INPUT(1:4)=='DATA')THEN
     call read_3d(Sali,TRIM(INI_SALI_FILE))

     DO K=1,Kloc
       CALL PHI_2D_COLL(SALI(:,:,K))
     ENDDO

     DO K=1,Kloc
     DO J=1,Nloc
     DO I=1,Mloc
       DSali(I,J,K)=Hc(I,J)*Sali(I,J,K)
     ENDDO
     ENDDO
     ENDDO

     call sali_bc(1)

    ENDIF ! end read salinity

     DSali0 = DSali

     IF(BOUNDARY(1:10)=='TID_FLX_LR'.OR.  &
        BOUNDARY(1:10)=='TID_ELE_LR') THEN     
       Sin_X0 = DataSal_L(1)
       Sin_Xn = DataSal_R(1)
     ENDIF

# endif
! end salinity


# if defined (TEMPERATURE)
    IF(INI_TEMP_INPUT(1:4)=='CONS')THEN
     TEMP = INI_TEMP
    ELSE IF(INI_TEMP_INPUT(1:4)=='DATA')THEN

      call read_3d(Temp,TRIM(INI_TEMP_FILE))
!      call sali_bc(1)

     DO K=1,Kloc
       CALL PHI_2D_COLL(TEMP(:,:,K))
     ENDDO

    ENDIF ! end read temp    
# endif
! end temperature

#if defined (SALINITY)
! gangfeng's code has

     call eval_dens
# endif


     Tke = Zero
     Eps = Zero
     DTke = Zero
     DEps = Zero
     DTke0 = Zero
     DEps0 = Zero

     ! wave breaking mask
     Brks = 0

     ! pressure boundary
     Bc_Prs = Zero

     ! initial surface elevation (user-specified)
     Eta = Zero


     ! collect data into ghost cells
     call phi_2D_coll(Eta)
     Eta0 = Eta

     call vel_bc(1)

# if defined (PARALLEL)
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)
# endif

     ! wetting-drying mask
     ! Mask: 1 - wet; 0 - dry
     ! Mask_Struct: 0 - permanent dry point
     ! Mask9: mask for itself and 8 elements around
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo
     Mask = Mask*Mask_Struct

# if defined (OBSTACLE)
     Mask3D=1
     U_mask=ZERO
     V_mask=ZERO
     W_mask=ZERO
     Ax_mask=ZERO
     Ay_mask=ZERO
     Az_mask=ZERO
! test     Mask3D(100:120,14:20,Kend-3:Kend)=0
     OPEN(4,FILE=TRIM(MASK3D_FILE),STATUS='OLD')  
       READ(4,*) ! title 
       READ(4,*) ! time character
       READ(4,*)TIME_READ_MASK
       READ(4,*,END=190)Number_mask
       READ(4,*,END=190) ! character
       DO I=1,Number_mask
# if defined(MASK3D_TOPLAYER)
          READ(4,*,END=190)Itmp,Jtmp, &
                          Surf_bottom_read,Surf_top_read, &
                          U_read,V_read,W_read,Ax_read,&
                           Ay_read,Az_read
# else
          READ(4,*,END=190)Itmp,Jtmp, Ktmp,&
                          U_read,V_read,W_read,Ax_read,&
                           Ay_read,Az_read
# endif
# if defined(MASK3D_TOPLAYER)
            ! remove masks out of domain
            IF(Itmp>0.AND.Itmp<=Mglob.AND. &
               Jtmp>0.AND.Jtmp<=Nglob)THEN
# else
            IF(Itmp>0.AND.Itmp<=Mglob.AND. &
               Jtmp>0.AND.Jtmp<=Nglob.AND. &
               Ktmp>0.AND.Ktmp<=Kglob)THEN
# endif
# if defined (PARALLEL)
               IItmp=-npx*(Mloc-2*Nghost)+Itmp+Nghost
               JJtmp=-npy*(Nloc-2*Nghost)+Jtmp+Nghost
# else
               IItmp=Itmp+Nghost
               JJtmp=Jtmp+Nghost
# endif
# if defined(MASK3D_TOPLAYER)
               IF(Surf_bottom_read==Surf_top_read)THEN
                 KKtmp=Kend
               ELSE
                 KKtmp=Kend-1
               ENDIF
# else
               KKtmp=Ktmp+Nghost
# endif

              IF(IItmp.GE.1.AND.IItmp.LE.Mloc.AND. &
                 JJtmp.GE.1.AND.JJtmp.LE.Nloc)THEN
                 Mask3D(IItmp,JJtmp,KKtmp)=0
# if defined (MASK3D_TOPLAYER)
                 Surf_Mask(IItmp,JJtmp) = .True.
                 Surf_Mask_bottom(IItmp,JJtmp)=Surf_bottom_read
                 Surf_Mask_top(IItmp,JJtmp)=Surf_top_read
# endif
                 U_mask(IItmp,JJtmp,KKtmp)=U_read
                 V_mask(IItmp,JJtmp,KKtmp)=V_read
                 W_mask(IItmp,JJtmp,KKtmp)=W_read
                 Ax_mask(IItmp,JJtmp,KKtmp)=Ax_read
                 Ay_mask(IItmp,JJtmp,KKtmp)=Ay_read
                 Az_mask(IItmp,JJtmp,KKtmp)=Az_read
              ENDIF
            ENDIF
        ENDDO
  
       READ(4,*,END=190) ! time character
       READ(4,*,END=190)TIME_READ_MASK
 
190    CONTINUE

# endif
! end obstacle

# if defined (PARALLEL)
     ! collect mask into ghost cells
     call phi_int_exch(Mask)
# if defined (OBSTACLE)
     CALL PHI_INT_EXCH3D(Mask3D)
# endif

# endif

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     ! total water depth and flux
     D = max(Hc+Eta, MinDep)

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     if(VISCOUS_FLOW) then
       ! initial seeding values for turbulence
!       Tke_min = 0.5*(1.4e-3)**2
!       Eps_min = 0.09*Tke_min**2/(0.1*Visc)
!       Tke_min = 1.e-12
!       Eps_min = 0.09*Tke_min**2/(1.e-4*Visc)
       Tke_min = 1.e-9
       Eps_min = 1.e-9
       Cmut_min = 9.e-2*Tke_min**2/Eps_min     
       do k = 1,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Tke(i,j,k) = Tke_min
         Eps(i,j,k) = Eps_min
         CmuHt(i,j,k) = Cmut_min
         CmuVt(i,j,k) = Cmut_min
         DTke(i,j,k) = D(i,j)*Tke(i,j,k)*Mask(i,j)
         DEps(i,j,k) = D(i,j)*Eps(i,j,k)*Mask(i,j)
       enddo
       enddo
       enddo
     endif



# if defined (SEDIMENT_lock_case)
     T = 500.
     bl = 4000.0

     ! lock-exchange case
     do i = Ibeg,Iend
     do j = Jbeg,Jend
     do k = Kbeg,Kend
!       if(xc(i)>=0.4) then
!         Sali(i,j,k) = 1.3592
!         DSali(i,j,k) = D(i,j)*Sali(i,j,k)
!       else
!         Sali(i,j,k) = 0.0
!         DSali(i,j,k) = D(i,j)*Sali(i,j,k)
!       endif

        Slope = atan(1./15.)
        if(xc(i)>=115000.-4000.*cos(Slope).and.xc(i)<=115000.+4000.*cos(Slope)) then
          Xslide = (xc(i)-115000.)/cos(Slope)
          Zslide(i,j) = T*sqrt(max(0.0,1.0-Xslide**2/bl**2))
        else
          Zslide(i,j) = 0.0
        endif  

        Zlev = sigc(k)*Hc(i,j)
        if(Zlev<=Zslide(i,j)) then
          Conc(i,j,k) = 1.0
        else
          Conc(i,j,k) = 0.0
        endif
        DConc(i,j,k) = D(i,j)*Conc(i,j,k)   
     enddo
     enddo
     enddo

     call sedi_bc(1)
# if defined (PARALLEL)
     call phi_3D_exch(DConc)
# endif

     DConc0 = DConc

     do i = 1,Mloc
     do j = 1,Nloc
     do k = 1,Kloc
       Conc(i,j,k) = DConc(i,j,k)/D(i,j)
     enddo
     enddo
     enddo

     call eval_dens
# endif

# if defined (LANDSLIDE)
     if(trim(SlideType)=='DEFORMABLE') then
       SlideAngle = SlideAngle*pi/180.
       conc_slide = (SlideDens-Rho0)/(SRho-Rho0)
       ls1 = SlideX0-0.5*SlideL*cos(SlideAngle)
       ls2 = SlideX0+0.5*SlideL*cos(SlideAngle)

       do i = Ibeg,Iend
       do j = Jbeg,Jend
         lsx = (xc(i)-SlideX0)/cos(SlideAngle)
         if(xc(i)>=ls1.and.xc(i)<=ls2) then
           hslide = SlideT*(1.0-(2.0*lsx/SlideL)**2)
           ! The following slide height is an approximation.
           Zslide(i,j) = hslide/cos(SlideAngle)
         else
           Zslide(i,j) = 0.0
         endif

         do k = Kbeg,Kend
            Zlev = sigc(k)*Hc(i,j)
            if(Zlev<=Zslide(i,j)) then
              Conc(i,j,k) = conc_slide
            else
              Conc(i,j,k) = 0.0
            endif
            DConc(i,j,k) = D(i,j)*Conc(i,j,k)
         enddo
       enddo
       enddo

!       ! lab case
!       do i = Ibeg,Iend
!       do j = Jbeg,Jend
!         if(xc(i)>=1.0.and.xc(i)<=1.65) then
!           Zslide(i,j) = xc(i)-1.0
!         else
!           Zslide(i,j) = 0.0
!         endif
!
!         do k = Kbeg,Kend
!            Zlev = sigc(k)*Hc(i,j)
!            if(xc(i)>=1.0.and.xc(i)<=1.65) then
!              if(Zlev<=Zslide(i,j)) then
!                Conc(i,j,k) = conc_slide
!              else
!                Conc(i,j,k) = 0.0
!              endif
!            endif
!            DConc(i,j,k) = D(i,j)*Conc(i,j,k)
!         enddo
!       enddo
!     enddo

       call sedi_bc
# if defined (PARALLEL)
       call phi_3D_exch(DConc)
# endif

       DConc0 = DConc

       do i = 1,Mloc
       do j = 1,Nloc
       do k = 1,Kloc
         Conc(i,j,k) = DConc(i,j,k)/D(i,j)
       enddo
       enddo
       enddo

       ! update mixture density
       call eval_dens

       ! update mixture viscosity
       do i = 1,Mloc
       do j = 1,Nloc
       do k = 1,Kloc
         if(Conc(i,j,k)<0.1) then
           Cmu(i,j,k) = Visc
         else
           Cmu(i,j,k) = Mud_Visc
         endif
       enddo
       enddo
       enddo

       ! rheology
       IsMove = 1
       if(RHEOLOGY_ON) then
         do i = 1,Mloc
         do j = 1,Nloc
         do k = 1,Kloc
           if(Conc(i,j,k)>=0.1) then
             CmuR(i,j,k) = 1.e+10
             IsMove(i,j,k) = 0
           else
             CmuR(i,j,k) = 0.0
             IsMove(i,j,k) = 1
           endif
         enddo
         enddo
         enddo
       endif
     endif
# endif     

if(HOTSTART)then
     call hot_start
endif


     ! SSP Runge-Kutta method parameters
     if(TIME_ORDER(1:3)=='THI') then
       It_Order = 3
       ALPHA(1) = 0.0
       ALPHA(2) = 3.0/4.0
       ALPHA(3) = 1.0/3.0
       BETA(1) = 1.0
       BETA(2) = 1.0/4.0
       BETA(3) = 2.0/3.0
     elseif(TIME_ORDER(1:3)=='SEC') then
       It_Order = 2
       ALPHA(1) = 0.0
       ALPHA(2) = 1.0/2.0
       BETA(1) = 1.0
       BETA(2) = 1.0/2.0
     else
       It_Order = 1
       ALPHA(1) = 0.0
       BETA(1) = 1.0
     endif

     ! wavemaker sponge
     IF(WAVEMAKER(1:3)=='ABS')THEN
       ALLOCATE(SPONGE_WAVEMAKER_LEFT(Mloc,Nloc),SPONGE_WAVEMAKER_RIGHT(Mloc,Nloc))
       CALL CALCULATE_WAVEMAKER_SPONGE
     ENDIF

     ! sponge layer
     if(SPONGE_ON) then
! gangfeng changed the one used in funwave, named calculate_sponge
!  fyshi
!       call calculate_sponge
! I'm trying the one in funwave
       call calculate_sponge_fyshi
       
     endif
    
     DEALLOCATE (TMP_READ)

END SUBROUTINE INITIALIZATION

     subroutine allocate_variables 
!--------------------------------------------------- 
!    This subroutine is used to allocate variables
!    Called by                                    
!       main                                                                                  
!    Last update: 23/12/2010, Gangfeng Ma                                                   
!---------------------------------------------------
     use global
# if defined (AMR)
     USE NESTING
# endif
     implicit none

!   Mloc Nloc Kloc should be maximum defined in index_local
# if defined (AMR)
     Mloc = MaxDimX
     Nloc = MaxDimY
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc = Kglob+2*Nghost
     Kloc1 = Kloc+1
# endif


# if defined (AMR)
! allocate nesting
     IF(NumGrid.GT.1)THEN
       ALLOCATE(U_NESTING_EAST(Nloc,Nghost,Kloc,2),&
                  V_NESTING_EAST(Nloc,Nghost,Kloc,2), &
                  Z_NESTING_EAST(Nloc,Nghost,2),&
                  U_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  V_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  Z_NESTING_WEST(Nloc,Nghost,2),&
                  U_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  V_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  Z_NESTING_SOUTH(Mloc,Nghost,2),&
                  U_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  V_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  Z_NESTING_NORTH(Mloc,Nghost,2), &
! add w p s t etc
                  W_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  P_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  S_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  T_NESTING_WEST(Nloc,Nghost,Kloc,2), &
                  W_NESTING_EAST(Nloc,Nghost,Kloc,2), &
                  P_NESTING_EAST(Nloc,Nghost,Kloc,2), &
                  S_NESTING_EAST(Nloc,Nghost,Kloc,2), &
                  T_NESTING_EAST(Nloc,Nghost,Kloc,2), &
                  W_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  P_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  S_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  T_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  W_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  P_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  S_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  T_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  TKE_NESTING_WEST(Mloc,Nghost,Kloc,2), &
                  TKE_NESTING_EAST(Mloc,Nghost,Kloc,2), &
                  TKE_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  TKE_NESTING_NORTH(Mloc,Nghost,Kloc,2), &
                  Eps_NESTING_WEST(Mloc,Nghost,Kloc,2), &
                  Eps_NESTING_EAST(Mloc,Nghost,Kloc,2), &
                  Eps_NESTING_SOUTH(Mloc,Nghost,Kloc,2), &
                  Eps_NESTING_NORTH(Mloc,Nghost,Kloc,2) &
                   )
     ENDIF
# endif
  ! end AMR


     ! one-dimensional vars
     ALLOCATE(x(Mloc1),xc(Mloc),y(Nloc1),yc(Nloc),sig(Kloc1),dsig(Kloc),sigc(Kloc),  &
              Ein_X0(Nloc),Din_X0(Nloc),Ein_Xn(Nloc),Din_Xn(Nloc),&
              Ein_Y0(Nloc),Din_Y0(Nloc),Ein_Yn(Nloc),Din_Yn(Nloc))

     ! two-dimensional vars
     ALLOCATE(HCG(Mglob,Nglob),Ho(Mloc,Nloc),H(Mloc1,Nloc1),Hc(Mloc,Nloc),Hc0(Mloc,Nloc),&
              Hfx(Mloc1,Nloc),Hfy(Mloc,Nloc1),D(Mloc,Nloc),  &
              Eta0(Mloc,Nloc),Eta00(Mloc,Nloc),D0(Mloc,Nloc),DeltH(Mloc,Nloc),&
              DelxH(Mloc,Nloc),DelyH(Mloc,Nloc),Eta(Mloc,Nloc),Mask(Mloc,Nloc),  &
# if defined(OBSTACLE)
              Mask3D(Mloc,Nloc,Kloc),U_mask(Mloc,Nloc,Kloc), &
              V_mask(Mloc,Nloc,Kloc),W_mask(Mloc,Nloc,Kloc),&
              Ax_mask(Mloc,Nloc,Kloc),Ay_mask(Mloc,Nloc,Kloc),&
              Az_mask(Mloc,Nloc,Kloc),&
              Surf_Mask_top(Mloc,Nloc),Surf_Mask_bottom(Mloc,Nloc), &
              Surf_Mask(Mloc,Nloc), &
# endif
              Mask_Struct(Mloc,Nloc),Mask9(Mloc,Nloc),SourceC(Mloc,Nloc),&
              SourceX(Mloc,Nloc),SourceY(Mloc,Nloc),DeltHo(Mloc,Nloc), &
              Uin_X0(Nloc,Kloc),Vin_X0(Nloc,Kloc),Win_X0(Nloc,Kloc),&
              Uin_Xn(Nloc,Kloc),Delt2H(Mloc,Nloc),  &
              Vin_Xn(Nloc,Kloc),Win_Xn(Nloc,Kloc),  &
              Uin_Y0(Nloc,Kloc),Vin_Y0(Nloc,Kloc),Win_Y0(Nloc,Kloc),&
              Uin_Yn(Nloc,Kloc), &
              Vin_Yn(Nloc,Kloc),Win_Yn(Nloc,Kloc), &
              Bc_Prs(Mloc,Nloc),Brks(Mloc,Nloc))
     ALLOCATE(DxL(Mloc1,Nloc),DxR(Mloc1,Nloc),DyL(Mloc,Nloc1),DyR(Mloc,Nloc1), &
              EtaxL(Mloc1,Nloc),EtaxR(Mloc1,Nloc),EtayL(Mloc,Nloc1),EtayR(Mloc,Nloc1), &
              DelxEta(Mloc,Nloc),DelyEta(Mloc,Nloc),DelxD(Mloc,Nloc),&
              DelyD(Mloc,Nloc),Sponge(Mloc,Nloc), &
              Setup(Mloc,Nloc),WaveHeight(Mloc,Nloc),Umean(Mloc,Nloc),&
              Vmean(Mloc,Nloc),Num_Zero_Up(Mloc,Nloc), &
              Emax(Mloc,Nloc),Emin(Mloc,Nloc))


     ! three-dimensional vars
     ALLOCATE(U(Mloc,Nloc,Kloc),V(Mloc,Nloc,Kloc),W(Mloc,Nloc,Kloc),Omega(Mloc,Nloc,Kloc1), &
              P(Mloc,Nloc,Kloc1),DU(Mloc,Nloc,Kloc),DV(Mloc,Nloc,Kloc),DW(Mloc,Nloc,Kloc),  &
              U0(Mloc,Nloc,Kloc),V0(Mloc,Nloc,Kloc),W0(Mloc,Nloc,Kloc),  &
              Un(Mloc,Nloc,Kloc),Vn(Mloc,Nloc,Kloc),&
                   Wn(Mloc,Nloc,Kloc),  &
              U00(Mloc,Nloc,Kloc),V00(Mloc,Nloc,Kloc),W00(Mloc,Nloc,Kloc),  &
              DU0(Mloc,Nloc,Kloc),DV0(Mloc,Nloc,Kloc),DW0(Mloc,Nloc,Kloc),Uf(Mloc,Nloc,Kloc1), &
              Vf(Mloc,Nloc,Kloc1),Wf(Mloc,Nloc,Kloc1),Cmu(Mloc,Nloc,Kloc), &
              CmuR(Mloc,Nloc,Kloc), &
              Diffxx(Mloc,Nloc,Kloc),Diffxy(Mloc,Nloc,Kloc),Diffxz(Mloc,Nloc,Kloc), &
              Diffyx(Mloc,Nloc,Kloc),Diffyy(Mloc,Nloc,Kloc),Diffyz(Mloc,Nloc,Kloc),  &
              Diffzx(Mloc,Nloc,Kloc), &
              Diffzy(Mloc,Nloc,Kloc),Diffzz(Mloc,Nloc,Kloc),DelxSc(Mloc,Nloc,Kloc),  &
              DelySc(Mloc,Nloc,Kloc), &
              CmuHt(Mloc,Nloc,Kloc),CmuVt(Mloc,Nloc,Kloc),Rho(Mloc,Nloc,Kloc),  &
              Rmean(Mloc,Nloc,Kloc),Tke(Mloc,Nloc,Kloc), &
              Eps(Mloc,Nloc,Kloc),Skl(Mloc,Nloc,Kloc),DTke(Mloc,Nloc,Kloc),  &
              DEps(Mloc,Nloc,Kloc),DTke0(Mloc,Nloc,Kloc), &
              DEps0(Mloc,Nloc,Kloc),Prod_s(Mloc,Nloc,Kloc),Prod_b(Mloc,Nloc,Kloc),  &
              Lag_Umean(Mloc,Nloc,Kloc), &
              Lag_Vmean(Mloc,Nloc,Kloc),Lag_Wmean(Mloc,Nloc,Kloc),  &
              Euler_Umean(Mloc,Nloc,Kloc),Euler_Vmean(Mloc,Nloc,Kloc), &
              Euler_Wmean(Mloc,Nloc,Kloc),DRhoX(Mloc,Nloc,Kloc),  &
              DRhoY(Mloc,Nloc,Kloc),ExtForceX(Mloc,Nloc,Kloc), &
              ExtForceY(Mloc,Nloc,Kloc),UpWp(Mloc,Nloc,Kloc),&
              IsMove(Mloc,Nloc,Kloc),Richf(Mloc,Nloc,Kloc))


     ! fluxes for construction at cell faces    
     ALLOCATE(UxL(Mloc1,Nloc,Kloc),UxR(Mloc1,Nloc,Kloc),VxL(Mloc1,Nloc,Kloc),  &
              VxR(Mloc1,Nloc,Kloc), &
              WxL(Mloc1,Nloc,Kloc),WxR(Mloc1,Nloc,Kloc),DUxL(Mloc1,Nloc,Kloc),  &
              DUxR(Mloc1,Nloc,Kloc), &
              DVxL(Mloc1,Nloc,Kloc),DVxR(Mloc1,Nloc,Kloc),DWxL(Mloc1,Nloc,Kloc),  &
              DWxR(Mloc1,Nloc,Kloc), &
              UyL(Mloc,Nloc1,Kloc),UyR(Mloc,Nloc1,Kloc),VyL(Mloc,Nloc1,Kloc),  &
              VyR(Mloc,Nloc1,Kloc), &
              WyL(Mloc,Nloc1,Kloc),WyR(Mloc,Nloc1,Kloc),DUyL(Mloc,Nloc1,Kloc),  &
              DUyR(Mloc,Nloc1,Kloc), &
              DVyL(Mloc,Nloc1,Kloc),DVyR(Mloc,Nloc1,Kloc),DWyL(Mloc,Nloc1,Kloc),  &
              DWyR(Mloc,Nloc1,Kloc), &
              UzL(Mloc,Nloc,Kloc1),UzR(Mloc,Nloc,Kloc1),VzL(Mloc,Nloc,Kloc1),  &
              VzR(Mloc,Nloc,Kloc1), &
              WzL(Mloc,Nloc,Kloc1),WzR(Mloc,Nloc,Kloc1),OzL(Mloc,Nloc,Kloc1),  &
              OzR(Mloc,Nloc,Kloc1), &
              SxL(Mloc1,Nloc,Kloc),SxR(Mloc1,Nloc,Kloc),SxS(Mloc1,Nloc,Kloc), &
              SyL(Mloc,Nloc1,Kloc),SyR(Mloc,Nloc1,Kloc),SyS(Mloc,Nloc1,Kloc), &
              ExL(Mloc1,Nloc,Kloc),ExR(Mloc1,Nloc,Kloc),FxL(Mloc1,Nloc,Kloc),  &
              FxR(Mloc1,Nloc,Kloc), &
              GxL(Mloc1,Nloc,Kloc),GxR(Mloc1,Nloc,Kloc),HxL(Mloc1,Nloc,Kloc),  &
              HxR(Mloc1,Nloc,Kloc), &
              EyL(Mloc,Nloc1,Kloc),EyR(Mloc,Nloc1,Kloc),FyL(Mloc,Nloc1,Kloc),  &
              FyR(Mloc,Nloc1,Kloc), &
              GyL(Mloc,Nloc1,Kloc),GyR(Mloc,Nloc1,Kloc),HyL(Mloc,Nloc1,Kloc),  &
              HyR(Mloc,Nloc1,Kloc), &
              Ex(Mloc1,Nloc,Kloc),Ey(Mloc,Nloc1,Kloc),&
              Ez(Mloc,Nloc,Kloc1),&
              Fx(Mloc1,Nloc,Kloc),  &
              Fy(Mloc,Nloc1,Kloc), &
              Gx(Mloc1,Nloc,Kloc),Gy(Mloc,Nloc1,Kloc),Hx(Mloc1,Nloc,Kloc),  &
              Hy(Mloc,Nloc1,Kloc), &
              Fz(Mloc,Nloc,Kloc1),Gz(Mloc,Nloc,Kloc1),Hz(Mloc,Nloc,Kloc1),&
              DelxU(Mloc,Nloc,Kloc), &
              DelyU(Mloc,Nloc,Kloc),DelzU(Mloc,Nloc,Kloc),DelxV(Mloc,Nloc,Kloc),  &
              DelyV(Mloc,Nloc,Kloc), &
              DelzV(Mloc,Nloc,Kloc),DelxW(Mloc,Nloc,Kloc),DelyW(Mloc,Nloc,Kloc),  &
              DelzW(Mloc,Nloc,Kloc), &
              DelxDU(Mloc,Nloc,Kloc),DelyDU(Mloc,Nloc,Kloc),DelxDV(Mloc,Nloc,Kloc),  &
              DelyDV(Mloc,Nloc,Kloc), &
              DelxDW(Mloc,Nloc,Kloc),DelyDW(Mloc,Nloc,Kloc),DelzO(Mloc,Nloc,Kloc)) 

# if defined (BUBBLE)
     ALLOCATE(Rbg(Mg),DRbg(Mg),Wbg(Mg),Entrain(Mg))
     ALLOCATE(Nbg(Mloc,Nloc,Kloc,Mg),DNbg(Mloc,Nloc,Kloc,Mg),DNbg0(Mloc,Nloc,Kloc,Mg),  &
              Vbg(Mloc,Nloc,Kloc))
# endif

# if defined (SEDIMENT)
     ALLOCATE(Wset(Mloc,Nloc,Kloc1),Conc(Mloc,Nloc,Kloc), &
            DConc(Mloc,Nloc,Kloc),DConc0(Mloc,Nloc,Kloc),  &
              SSour(Mloc,Nloc),Bed(Mloc,Nloc),Bed0(Mloc,Nloc),&
            Taub(Mloc,Nloc),Sed_X0(Nloc,Kloc),Sed_Xn(Nloc,Kloc))
# endif

# if defined (SALINITY)
     ALLOCATE(Sali(Mloc,Nloc,Kloc),DSali(Mloc,Nloc,Kloc),DSali0(Mloc,Nloc,Kloc))
     ALLOCATE(Sin_X0(Nloc,Kloc),Sin_Xn(Nloc,Kloc))
# endif

# if defined (TEMPERATURE)
     ALLOCATE(Temp(Mloc,Nloc,Kloc),DTemp(Mloc,Nloc,Kloc),DTemp0(Mloc,Nloc,Kloc))                                    
     ALLOCATE(Tin_X0(Nloc,Kloc),Tin_Xn(Nloc,Kloc))
# endif

# if defined (AIR_PRESSURE)
     ALLOCATE(Pmb(Mloc,Nloc))
# endif

# if defined (VEGETATION)
     ALLOCATE(Prod_v(Mloc,Nloc,Kloc))
# endif

! move to index_local
!     ! poisson solver (for NSPCG use)
!     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
!     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))


     end subroutine allocate_variables
     

     subroutine index_dont_use_amr
!---------------------------------------------------
!    This subroutine is used to creat work index
!    Called by                   
!       main                                                                
!    Last update: 20/12/2010, Gangfeng Ma                                     
!---------------------------------------------------
     use global
     implicit none

# if defined (PARALLEL)
     dims(1)=PX
     dims(2)=PY
     periods(1)=.false.
     periods(2)=.false.
     if(PERIODIC_X) periods(1)=.true.
     if(PERIODIC_Y) periods(2)=.true.
     coords(1)=0
     coords(2)=0

     call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims, &
         periods,reorder,comm2d,ier)
     call MPI_CART_COORDS(comm2d,myid,2,coords,ier)

     npx=coords(1)
     npy=coords(2)
 
     call MPI_CART_SHIFT(comm2d,0,1,n_west,n_east,ier)
     call MPI_CART_SHIFT(comm2d,1,1,n_suth,n_nrth,ier)
# else
     PX=1
     PY=1
# endif

     ! local index
     Mloc = Mglob/PX+2*Nghost
     Nloc = Nglob/PY+2*Nghost
     Kloc = Kglob+2*Nghost
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc1 = Kloc+1

     Ibeg = Nghost+1
     Iend = Mloc-Nghost
     Iend1 = Mloc1-Nghost
     Jbeg = Nghost+1
     Jend = Nloc-Nghost
     Jend1 = Nloc1-Nghost
     Kbeg = Nghost+1
     Kend = Kloc-Nghost
     Kend1 = Kloc1-Nghost

     end subroutine index_dont_use_amr

SUBROUTINE INDEX_LOCAL(ng)
# if defined (PARALLEL)
    USE GLOBAL,ONLY: PX,PY,Mloc,Nloc,Mglob,Nglob,Nghost,Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,&
                     Mloc1,Nloc1,Kbeg,Kend,Kend1,Kloc,Kloc1,Kglob
# else
    USE GLOBAL,ONLY: Mloc,Nloc,Mglob,Nglob,Nghost,Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,&
                     Mloc1,Nloc1,Kbeg,Kend,Kend1,Kloc,Kloc1,Kglob
# endif
# if defined (AMR)
    USE NESTING, ONLY: GridDimX,GridDimY
# endif

    USE GLOBAL,ONLY : Coef,JCoef,Rhs,neqns

    IMPLICIT NONE
    INTEGER,INTENT(IN)::ng
    
# if defined(PARALLEL)
# if defined (AMR)
    Mloc=GridDimX(ng)/px+2*Nghost
    Nloc=GridDimY(ng)/py+2*Nghost
# else
    Mloc=Mglob/px+2*Nghost
    Nloc=Nglob/py+2*Nghost
# endif
    Mloc1=Mloc+1
    Nloc1=Nloc+1
# else
# if defined (AMR)
    Mloc=GridDimX(ng)+2*Nghost
    Nloc=GridDimY(ng)+2*Nghost
# else
    Mloc=Mglob+2*Nghost
    Nloc=Nglob+2*Nghost
# endif
    Mloc1=Mloc+1
    Nloc1=Nloc+1
# endif

    Kloc = Kglob+2*Nghost
    Kloc1 = Kloc+1

    Ibeg=Nghost+1
    Iend=Mloc-Nghost
    Iend1=Mloc1-Nghost
    Jbeg=Nghost+1
    Jend=Nloc-Nghost
    Jend1=Nloc1-Nghost

    Kbeg = Nghost+1
    Kend = Kloc-Nghost
    Kend1 = Kloc1-Nghost

     ! poisson solver (for NSPCG use)
     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
     IF(ALLOCATED(Coef))THEN
      DEALLOCATE(Coef,JCoef,Rhs)
     ENDIF
     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))

END SUBROUTINE INDEX_LOCAL
  


# if defined (PARALLEL)
! ----------------------------------------------------
!    This is subroutine to index for MPI
!  called by 
!        MAIN
!    Last Update: 05/06/2010 Fengyan Shi, University of Delaware
! --------------------------------------------------

SUBROUTINE PARALLEL_CARTESIAN
    USE GLOBAL,ONLY: PX,PY,dims,periods,coords,ndims,comm2d,ier,myid,&
                     coords,n_west, n_east, n_suth, n_nrth,&
                     NumP,npx,npy,PERIODIC_X,PERIODIC_Y,reorder
    USE MPI
    IMPLICIT NONE

    NumP = px*py
    dims(1) = px
    dims(2) = py
    periods(1) = .false.
    periods(2) = .false.
     if(PERIODIC_X) periods(1)=.true.
     if(PERIODIC_Y) periods(2)=.true.
    coords(1) = 0
    coords(2) = 0

    call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, &
         periods, reorder, comm2d, ier )
    call MPI_CART_COORDS( comm2d, myid, 2, coords, ier)

    npx = coords(1)
    npy = coords(2)

    call MPI_Cart_shift( comm2d, 0, 1, n_west, n_east, ier )
    call MPI_Cart_shift( comm2d, 1, 1, n_suth, n_nrth, ier )


END SUBROUTINE PARALLEL_CARTESIAN
# endif


     subroutine generate_grid
!------------------------------------------------------
!    This subroutine is used to generate grids
!    Called by
!       master
!    update: 20/12/2010, Gangfeng Ma
!    update 05/27/2013, fyshi
!    for AMR version
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,k

     ! horizontal grid
# if defined (PARALLEL)
     x(Ibeg) = npx*(Mloc-2*Nghost)*dx
# else
     x(Ibeg) = Zero
# endif
     do i = Ibeg+1,Mloc1
       x(i) = x(i-1)+dx
       xc(i-1) = x(i-1)+0.5*dx
     enddo
     do i = Ibeg-1,Ibeg-Nghost,-1
       x(i) = x(i+1)-dx
       xc(i) = x(i+1)-0.5*dx
     enddo

# if defined (PARALLEL)
     y(Jbeg) = npy*(Nloc-2*Nghost)*dy
# else
     y(Jbeg) = Zero
# endif
     do j = Jbeg+1,Nloc1
       y(j) = y(j-1)+dy
       yc(j-1) = y(j-1)+0.5*dy
     enddo
     do j = Jbeg-1,Jbeg-Nghost,-1
       y(j) = y(j+1)-dy
       yc(j) = y(j+1)-0.5*dy
     enddo

     ! vertical grid
     if(Ivgrd==1) then
       do k = 1,Kloc
         dsig(k) = 1.0/float(Kglob)
       enddo
     elseif(Ivgrd==2) then  ! fyshi added gangfeng's new one 05/13/2013
       dsig(Kbeg) = (Grd_R-1.0)/(Grd_R**float(Kglob)-1.0)
       do k = Kbeg+1,Kend
         dsig(k) = dsig(k-1)*Grd_R
       enddo
       do k = 1,Nghost
         dsig(Kbeg-k) = dsig(Kbeg+k-1)
       enddo
       do k = 1,Nghost
         dsig(Kend+k) = dsig(Kend-k+1)         
       enddo
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
       sig(k) = sig(k-1)+dsig(k-1)
       sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
       sig(k) = sig(k+1)-dsig(k)
       sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo

     end subroutine generate_grid

     subroutine CALCULATE_WAVEMAKER_SPONGE
!-------------------------------------------------
!    Calculate sponge function for wavemaker
!    Called by
!      initial
!    Last update: 06/10/2014 fyshi
!------------------------------------------------
     use global
     implicit none
     integer :: i,j,Iwidth
     real(SP) :: ri,lim,Lstart,Lend

     lim=1.0_SP

     if(WaveMaker_West_Width>Zero) then
       Iwidth = int(WaveMaker_West_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_WaveMaker**(50*(i+npx*Mglob/px-1)/(Iwidth-1))
# else
         ri = R_WaveMaker**(50*(i-1)/(Iwidth-1))
# endif
         SPONGE_WAVEMAKER_LEFT(i,j) = max(A_WaveMaker**ri,lim)
       enddo
       enddo
     endif


     if(WaveMaker_East_Width>Zero)then
       Iwidth = int(WaveMaker_East_Width/dx)+Nghost
       do j = 1,Nloc
       do i = 1,Mloc
# if defined(PARALLEL)
         ri = R_WaveMaker**(50*(Mloc-i+(px-npx-1)*Mglob/px)/(Iwidth-1))
# else
         ri = R_WaveMaker**(50*(Mloc-i)/(Iwidth-1))
# endif
         SPONGE_WAVEMAKER_RIGHT(i,j) = max(A_WaveMaker**ri,lim)
       enddo
       enddo
     endif

      SPONGE_WAVEMAKER_LEFT=1.0_SP/SPONGE_WAVEMAKER_LEFT
      SPONGE_WAVEMAKER_RIGHT=1.0_SP/SPONGE_WAVEMAKER_RIGHT

END SUBROUTINE CALCULATE_WAVEMAKER_SPONGE

