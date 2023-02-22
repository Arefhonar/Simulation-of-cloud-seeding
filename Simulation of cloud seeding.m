clc
clear
%% const
 
nx=10;
nt=10;
tf=3000;
XF=1;
tol=1e-6;
x=linspace(0,XF,nx);
t=linspace(0,tf,nt);
dt=max(diff(t));
dz=max(diff(x));
for i=1:nx
alfa(i)=-0.0155+0.0085*x(i);
end
Tf=zeros(nx,nt);
N=zeros(nx,nt);
Tu=zeros(nx,nt);
 
%%%% I.C for Temp
Tf(1,:)=-14.3238;
Tf(:,1)=-;
 
%%%% I.C for N
N0=zeros(nx)+2.097+0.2433*Tf(1,1)+0.006765*Tf(1,1)^2;
%% regrasion 
x1=[-31.0889 -29.2316 -28.2718 -27.3120 -26.4176 -25.7774 -24.8193 -23.8599 -22.3262 -21.4313 -19.7057 -17.9788 -17.1489 -16.3814 -14.6558 -13.7605 -12.2905];
y=[1.16712 0.85681 0.62901 0.46177 0.26886 0.23036 0.12415 0.08437 0.03606 0.02268 0.00897 0.00448 0.00241 0.00177 0.00070 0.00048 0.00022];
 
real_time=[0 43.4355 45.1113 53.6464 48.2617 50.0045 67.41 62.7625 91.8089 301.479 727.679 1056.1 1496.22 1908.37 2320.45 2536.51 2667.78 2735.83 2803.76 2844.18 2891.19 2931.67]
real_N=[0 13.831 19.9806 28.0303 39.8977 58.4915 123.903 207.346 262.493 255.668 235.538 223.152 202.621 183.904 164.477 140.345 98.783 65.4943 42.1642 30.0772 19.6442 14.2206]
 
Msserr=0.007358;
tcrit=2.145;
n=17;
m=3;
 
z=zeros(n,m);
z(:,1)=1;
z(:,2)=x1;
z(:,3)=x1.*x1;
a=(z'*z)^-1*z'*y'
C=(z'*z)^-1;
 
sigmaa0=(Msserr*C(1,1))^0.5
sigmaa1=(Msserr*C(2,2))^0.5
sigmaa2=(Msserr*C(3,3))^0.5
tcala0=a(1)/sigmaa0
tcala1=a(2)/sigmaa1
tcala2=a(3)/sigmaa2
 
tmina0=a(1)-tcrit*sigmaa0;tmaxa0=a(1)+tcrit*sigmaa0;
tmina1=a(2)-tcrit*sigmaa1;tmaxa1=a(2)+tcrit*sigmaa1;
tmina2=a(3)-tcrit*sigmaa2;tmaxa2=a(3)+tcrit*sigmaa2;
 
fprintf('%f <a0<%f\n',tmina0,tmaxa0)
fprintf('%f <a1<%f\n',tmina1,tmaxa1)
fprintf('%f <a2<%f\n',tmina2,tmaxa2)

%% method
%%% temp Calculation
for n=1:nt-1
    for i=1:nx-2
        if i==1
       A(i,i)=1;
       A(i,i+1)=-alfa(i)*dt/(2*dz*(2*a(3)*Tu(i)+a(2)));
       B(i,1)=Tu(i)-alfa(i)*dt/(2*dz*(2*a(3)*Tu(i)+a(2)));
        else
            A(i,i)=1;
            A(i,i+1)=-alfa(i)*dt/(2*dz*(2*a(3)*Tu(i)+a(2)));
            A(i,i-1)=alfa(i)*dt/(2*dz*(2*a(3)*Tu(i)+a(2)));
            B(i,1)=Tu(i);
        end
    end
   y=A\B;
   Tu=y;
   Tf(2:nx,n+1)=y;
   
end
%%%% N calculation
for i=1:nt
    for j=1:nx-1
        N(j,i)=alfa(i)*(Tf(j+1,i)-Tf(j,i))/dz*t(j)+N0(j);
    end
end
%% Plotting Result
subplot(2,1,1)
[Time,X]=meshgrid(t,x);
surf(Time,X,N')
title('N(Time,Z)when -31.0889<Temp<-12.2905')
xlabel('Time')
ylabel('Z')
zlabel('N')
 
%%%% if (N(k,:) & k=2n)
subplot(2,1,2) 
plot(t,N(7,:))
hold on
plot(real_time,real_N)
 
legend('N','real N')
title('N(Time) when z=0.6667 & -31.0889<Temp<-12.2905')
xlabel('Time')
ylabel('N')