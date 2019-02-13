fdir = '/Volumes/2TB_element/MOVING_PRESSURE/h_100m_u1p2m/';

m=2048;
l=100;
dx=5.0;
dz=1.0;
x=[0:m-1]*dx;
z=[0:l-1]*dz;
xz=meshgrid(x,z);

%ns=input('ns=');
%ne=input('ne=');

files=[21 41 61 81];

wid=8;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

icount=0;
for num=1:length(files)

icount=icount+1;

fnum=sprintf('%.4d',files(num));
sali=load([fdir 'sali_' fnum]);

subplot(length(files),1,num)


contourf(x,z,sali,5)
axis([min(x) max(x) 50 80])
if(num==length(files))
xlabel(' x(m) ')
end

time=['TIME = ', num2str((files(num)-1)*100), ' sec'];
title(time)

ylabel(' z(m) ')

end


