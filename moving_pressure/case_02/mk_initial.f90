         parameter(m=512,n=1,m1=720,n1=480,l=30,l1=39)
         real,dimension(m,n) :: var2D
         real,dimension(m,n,l) :: var3D,Zsigma
         real,dimension(m1,n1) :: var2D_fine
         real,dimension(m1,n1,l1)  :: var3D_fine
         real :: eta_left,eta_right,Temp,Sali,Usea,r,dx,dy,x_c,y_c,sigma,a
         real :: S_top,S_mid,S_bot,Lev_mid,dep
         real,dimension(m) :: xx
         real,dimension(n) :: yy
         integer :: i,j,nk,k

         eta_left = 0.0
         eta_right = -0.0
         Temp=15.0
         Sali=25.0
         Usea=11.0

         S_top = 23.0
         S_mid = 23.2
         S_bot = 24.0
         Lev_mid = 3.0
         dep=10.0

         do k=1,l
         do j=1,n
         do i=1,m
           Zsigma(i,j,k)=dep*(1.0-(k-1.0)/(l-1.0))
         enddo
         enddo
         enddo

         dx=5.0
         dy=5.0
         x_c=1000.0
! width (sigma /2, in meters)
         sigma = 200.0
! height (m)
         a=1.0
         

         do i=1,m
           xx(i)=(i-1.0)*dx
         enddo
         do j=1,n
           yy(i)=(j-1.0)*dy
         enddo

         do j=1,n
         do i=1,m
           var2D(i,j) = eta_left +(i-1.0)/(m-1.0)*(eta_right-eta_left)
         enddo
         enddo
         
! eta 

100      format(1000f12.3)

         open(2,file='eta_init.ini')

         do j=1,n
          write(2,100)(var2D(i,j),i=1,m)
         enddo
         close(2)

! Pmb
         do j=1,n
         do i=1,m
           r=xx(i)-x_c
           Var2D(i,j)=a*exp(-r**2/(sigma/4.0)**2)          
         enddo
         enddo

         open(2,file='pmb_init.ini')

         do j=1,n
          write(2,100)(var2D(i,j),i=1,m)
         enddo
         close(2)


         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=Usea
         enddo
         enddo
         enddo

!     

         open(2,file='u_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=0.0
         enddo
         enddo
         enddo


         open(2,file='v_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         open(2,file='w_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         open(2,file='p_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         open(2,file='eps_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         open(2,file='tke_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)



         open(2,file='mu_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         open(2,file='prod_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)


         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=1024.0984
         enddo
         enddo
         enddo

         open(2,file='rho_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=Sali
         enddo
         enddo
         enddo

! sali

         do k=1,l
         do j=1,n
         do i=1,m
           if(Zsigma(i,j,k).ge.Lev_mid)then
! linear     Var3D(i,j,k)=S_bot+(S_mid-S_bot)*(dep-Zsigma(i,j,k))/(dep-Lev_mid)
     Var3D(i,j,k)=S_mid+(S_bot-S_mid)*tanh((Zsigma(i,j,k)-Lev_mid)/dep*10.0)
           else
     Var3D(i,j,k)=S_mid+(S_top-S_mid)*(Lev_mid-Zsigma(i,j,k))/(Lev_mid)    
           endif
         enddo
         enddo
         enddo

         open(2,file='sali_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)


         do k=1,l
         do j=1,n
         do i=1,m
           var3D(i,j,k)=Temp
         enddo
         enddo
         enddo

         open(2,file='temp_init.ini')

         do k=1,l
         do j=1,n
          write(2,100)(var3D(i,j,k),i=1,m)
         enddo
         enddo
         close(2)

         end

          subroutine conv2D(m,n,m1,n1,nk,depc,dep4)
 
          integer, intent(in) :: nk,m,n,m1,n1
         real,dimension(m,n),intent(in) :: depc
         real,dimension(m1,n1),intent(out) :: dep4
         integer :: i,j,ii,jj,kx,ky
         real :: d_down,d_up



          do j=1,n-1 
          do i=1,m-1 
          do kx=1,nk
          do ky=1,nk
            jj=(j-1)*nk+ky
            ii=(i-1)*nk+kx

            if(depc(i,j).le.0.0.and.depc(i+1,j).gt.0.0)then
              if(kx.eq.1)then
                d_down=depc(i,j)
              else
                d_down=depc(i+1,j)
              endif        
            elseif(depc(i,j).gt.0.0.and.depc(i+1,j).le.0.0)then
              d_down=depc(i,j)
            else
            d_down=depc(i,j)+(kx-1.)/(nk)*(depc(i+1,j)-depc(i,j))
            endif

!        up
            if(depc(i,j+1).le.0.0.and.depc(i+1,j+1).gt.0.0)then
              if(kx.eq.1)then
                d_up=depc(i,j+1)
              else
                d_up=depc(i+1,j+1)
              endif             
            elseif(depc(i,j+1).gt.0.0.and.depc(i+1,j+1).le.0.0)then
              d_up=depc(i,j+1)
            else
            d_up=depc(i,j+1)+(kx-1.)/(nk)*(depc(i+1,j+1)-depc(i,j+1))
            endif

!         dep4
            if(d_down.le.0.0.and.d_up.gt.0.0)then
              if(ky.eq.1)then
                dep4(ii,jj)=d_down
              else
                dep4(ii,jj)=d_up
              endif             
            elseif(d_down.gt.0.0.and.d_up.le.0.0)then
              dep4(ii,jj)=d_down
            else
            dep4(ii,jj)=d_down+(ky-1.)/(nk)*(d_up-d_down)
            endif
          enddo
          enddo
          enddo
          enddo

         do j=1,(n-1)*nk+1
         dep4(nk*(m-1)+1,j)=dep4(nk*(m-1),j)
         enddo

         do i=1,(m-1)*nk+1
         dep4(i, nk*(n-1)+1)=dep4(i,nk*(n-1))
         enddo

         do j=1,n1
         dep4(m1,j)=dep4(m1-1,j)
         enddo

         do i=1,m1
         dep4(i, n1)=dep4(i,n1-1)
         enddo

         return

         end subroutine conv2D






