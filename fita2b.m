% Fit a signal to another
close all
Input = readtable('Input (CH1) vs Prediction(CH2)_40ksample WO Noise Sine wave Stanford.csv');
l2 = height(Input)
Input = table2array(Input);
%time = Input(:,1);
Value = Input(:,1);

T0 = 20e-6;
Tin = 2.5e-6;

time = Tin*[1:l2];
% f = 1/(time(2)-time(1))



plot(time,Value)
title('Input signal')

figure


In1 = Xr7;
l1 = length(In1);

plot(T0*[1:l1],In1)
title('sampled signal')

In2 = Value;
l2 = length(In2);
% t1 = 0:0.1:10;
% t2 = 0:0.01:12;
% 
% In1 = 2000*sin(2*pi*1*t1)+2000;
% In2 = 1.5*sin(2*pi*1*t2);

G0 = (max(In1)-min(In1))/(max(In2)-min(In2))


G  = G0
DC = max(In1)-max(G0*In2)

G  = 119.439/1000
DC = 2244.30

shift = 509;
ratio = 8;

param=[DC,G,shift];

arg = shift+ratio*[1:l1];
Value_c=G*Value(arg)+DC;



figure
plot(T0*[1:l1],In1,T0*[1:l1],Value_c)
title('In vs sampled')
legend('Sampled','Input')

E2=sqrt(mean((Value_c-In1).^2))
% 
% a = 1
% b= 0
% 
% mx=@(param)fita(a,b,shift,In1,Value_c,ratio)
% 
% outputparameters=fminsearch(mx,[a,b])
% 
% function E=fita(a,b,shift,xm,ym,R)
% 
% l1 = length(xm);
% l2 = length(ym);
% 
% 
% 
% gs=a*ym+b;
% E=sum(abs(gs-xm).^2)
% E2=sqrt(mean((gs-xm).^2))
% 
% 
% clf;
% plot(gs,'.','linewidth',200);
% hold on;
% plot(xm,'r');
% %axis([min(xm) max(xm) min(ym) max(ym)]);
% drawnow;
% end
