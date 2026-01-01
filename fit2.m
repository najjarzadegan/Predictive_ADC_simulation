
x = 2:4096;
y = Xr7(2:end);
y = dac';
yavg = 2368;

Py = sum((y-yavg).^2);

amplitude=1750;
freq=2*pi/50;
phase=2*pi/2;
offset=2150;

LB = [1300, 2*pi/10, 0, 2250];
UB = [1400, 2*pi/10, 2*pi, 2350];

initialparameter=[amplitude,freq,phase,offset];
mx=@(initialparameter)fita(initialparameter,x,y,Py)
%outputparameters=fminsearchbnd(mx,initialparameter,LB,UB)

outputparameters=fminsearch(mx,initialparameter)

function E=fita(initialparameter,xm,ym,P)
gs=initialparameter(1)*sin(initialparameter(2)*xm+initialparameter(3))+initialparameter(4);
E=sum(abs(gs-ym).^2)
E2=sqrt(mean((gs-ym).^2))
SNR=10*log10(P/E)

ya=initialparameter(1)*sin(initialparameter(2)*xm+initialparameter(3))+initialparameter(4);
clf;
plot(xm,ym,'.','linewidth',200);
hold on;
plot(xm,ya,'r');
axis([min(xm) max(xm) min(ym) max(ym)]);
drawnow;
end