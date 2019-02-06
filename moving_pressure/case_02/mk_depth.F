         parameter(m=512,n=1,m1=720,n1=480,l=40,l1=39)
         real,dimension(m,n) :: var2D
         real,dimension(m,n,l) :: var3D
         real,dimension(m1,n1) :: var2D_fine
         real,dimension(m1,n1,l1)  :: var3D_fine
         real :: eta_left,eta_right,Temp,Sali,Usea
         integer :: i,j,nk

         do j=1,n
         do i=1,m
           var2D(i,j) = 10.0
         enddo
         enddo
         
! eta 

100      format(2000f12.3)

         open(2,file='depth.txt')

         do j=1,n
          write(2,100)(var2D(i,j),i=1,m)
         enddo
         close(2)

         end






