fdir = '/Users/fengyanshi15/tmp2/';

m=512;
l=30;

ns=input('ns=');
ne=input('ne=');

%ns=1;
%ne=1;

% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

wid=6;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

icount=0;
for num=ns:1:ne

icount=icount+1;

fnum=sprintf('%.4d',num);
sali=load([fdir 'sali_' fnum]);
eta=load([fdir 'eta_' fnum]);

subplot(211)
plot(eta')
grid
axis([0 512 -2 6])
subplot(212)
contourf(sali,50)
axis([0 512 10 30])

    currframe=getframe(gcf);
    writeVideo(vidObj,currframe); 

pause(0.1)

end
close(vidObj)


