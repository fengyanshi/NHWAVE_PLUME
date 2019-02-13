fdir = '/Volumes/2TB_element/MOVING_PRESSURE/h_100m_u1p2m/';

m=2048;
l=100;
dx=5.0;
dz=1.0;
x=[0:m-1]*dx;
z=[0:l-1]*dz;
xz=meshgrid(x,z);

ns=input('ns=');
ne=input('ne=');

%ns=1;
%ne=1;

% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

wid=8;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

icount=0;
for num=ns:1:ne

icount=icount+1;

fnum=sprintf('%.4d',num);
sali=load([fdir 'sali_' fnum]);
eta=load([fdir 'eta_' fnum]);

subplot(4,1,[1])
plot(x,eta)
grid
axis([min(x) max(x) -10 10])
zlabel(' eta (m) ')
time=['TIME = ', num2str((num-1)*100), ' sec'];
title(time)

subplot(4,1,[2 3 4])
contourf(x,z,sali,10)
axis([min(x) max(x) 40 85])
xlabel(' x(m) ')
zlabel(' z(m) ')

    currframe=getframe(gcf);
    writeVideo(vidObj,currframe); 

pause(0.1)

end
close(vidObj)


