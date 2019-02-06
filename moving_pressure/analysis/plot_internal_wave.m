fdir = '/Users/fengyanshi15/tmp2/';

m=1024;
l=100;

ns=input('ns=');
ne=input('ne=');

%ns=1;
%ne=1;

% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

wid=8;
len=10;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

icount=0;
for num=ns:10:ne

icount=icount+1;

fnum=sprintf('%.4d',num);
sali=load([fdir 'sali_' fnum]);
clf
contourf(sali,100)
axis([0 1024 73 83])

    currframe=getframe(gcf);
    writeVideo(vidObj,currframe); 

pause(0.1)

end
close(vidObj)


