function lsdifffig
% LSDIFFFIG  make least-squares differencing illustration for staggered grid
% (appearance is tuned for Octave "print -dpdf")

hp = 0.1;

s = 0.3;
x = -1-s:hp:2+s;
xzero=zeros(size(x));
xone=ones(size(x));

y = -1-s:hp:1+s;
yzero=zeros(size(y));
yone=ones(size(y));

figure(1), clf
hold on
plot(x,-1*xone,'k',x,xzero,'k',x,1*xone,'k')
plot(-1*yone,y,'k',yzero,y,'k',1*yone,y,'k',2*yone,y,'k')
dot(-1,0)
dot(0,0)
dot(1,0)
dot(2,0)
dot(0,-1)
dot(0,1)
dot(1,-1)
dot(1,1)
dot(0.5,0,'^k',[16])
Q = 0.2;
R = 0.5;
txt(-1-Q,-1-R,'i-1')
txt(0-Q,-1-R,' i ')
txt(1-Q,-1-R,'i+1')
txt(2-Q,-1-R,'i+2')
R = 0.8;
txt(-1-R,-1,'j-1')
txt(-1-R,0,' j')
txt(-1-R,1,'j+1')
num(0,-1,1)
num(1,-1,2)
num(-1,0,3)
num(0,0,4)
num(1,0,5)
num(2,0,6)
num(0,1,7)
num(1,1,8)
hold off
axis([-2 3 -2 2]), axis off
end

  function dot(x,y,shape,sizes)
    if nargin < 3, shape='ok'; end
    if nargin < 4, sizes=1:1:12; end
    for k=1:length(sizes)
      plot(x,y,shape,'markersize',sizes(k))
    end
  end

  function txt(x,y,mystr)
    text(x,y,mystr,'fontsize',18,'fontname','Courier','fontweight','bold')
  end

  function num(i,j,N)
    r = 0.1; u = 0.1;
    text(i+r,j+u,num2str(N),'fontsize',12,'fontname','Helvetica')
  end

