
clear
clc

r2d=180/pi;
N=64;
R=10.;
sx=zeros(N+1,1);
cx=zeros(N,  1);
th = linspace(0,2*pi,N+1)';

dth=2*pi/N;
sx=R*cos(th-0.5*dth);
sy=R*sin(th-0.5*dth);

for i = 1:N
    cx(i)=0.5*(sx(i)+sx(i+1));
    cy(i)=0.5*(sy(i)+sy(i+1));
end

plot(cx,cy,'o','Linewidth',2); hold on; 
plot(sx,sy,'Linewidth',2); axis equal;
axis([-11, 11, -11, 11]);
t=[0:0.001:2*pi]; x=R*cos(t); y=R*sin(t);
plot(x,y,'-');
xlabel('x'); ylabel('y'); axis off
for i=1:N
    text(cx(i)+0.5, cy(i), num2str(i));
end
Xtext=sprintf('Hit a key if # of panel (%d) is good',N);
disp(Xtext);
pause;


windvel_x=1.0;
windvel_y=0.0;
U_inf=sqrt(windvel_x^2+windvel_y^2);
wind_ang=atan2( windvel_y,windvel_x);
tcp_ang=atan2(diff(sy),diff(sx));
ncp_ang=tcp_ang-pi/2;
beta=ncp_ang-wind_ang;


I_ij=zeros(N,N); A_ij=zeros(N,N); S_ij=zeros(N,N);
for i = 1:N
    for j=[1:(i-1) (i+1):N]
        A=-(cx(i)-sx(j))*cos(tcp_ang(j))-(cy(i)-sy(j))*sin(tcp_ang(j));
        B=(cx(i)-sx(j))^2+(cy(i)-sy(j))^2;
        C=sin(tcp_ang(i)-tcp_ang(j));
        D=(cy(i)-sy(j))*cos(tcp_ang(i))-(cx(i)-sx(j))*sin(tcp_ang(i));
        S=sqrt((sx(j+1)-sx(j))^2+(sy(j+1)-sy(j))^2);
        E=(cx(i)-sx(j))*sin(tcp_ang(j))-(cy(i)-sy(j))*cos(tcp_ang(j));
        I_ij(i,j)=0.5*C*log((S^2+2*A*S+B)/B)+((D-A*C)/E)*(atan((S+A)/E)-atan(A/E));
    end
end

for i=1:N
    A_ij(i,i)=pi;
end

A_ij=A_ij-I_ij;
b_vec=-2*pi*U_inf*cos(beta);
Lambda=A_ij\b_vec;


for i = 1:N
    for j=[1:(i-1) (i+1):N]
        A=-(cx(i)-sx(j))*cos(tcp_ang(j))-(cy(i)-sy(j))*sin(tcp_ang(j));
        B=(cx(i)-sx(j))^2+(cy(i)-sy(j))^2;
        C=sin(tcp_ang(i)-tcp_ang(j));
        D=(cy(i)-sy(j))*cos(tcp_ang(i))-(cx(i)-sx(j))*sin(tcp_ang(i));
        S=sqrt((sx(j+1)-sx(j))^2+(sy(j+1)-sy(j))^2);
        E=(cx(i)-sx(j))*sin(tcp_ang(j))-(cy(i)-sy(j))*cos(tcp_ang(j));
        S_ij(i,j)=((D-A*C)/(2*E))*log((S^2+2*A*S+B)/B)-C*(atan((S+A)/E)-atan(A/E));
    end
end

vel_cp=U_inf*sin(beta)-S_ij*Lambda/(2*pi);
Cp=1. - (vel_cp/U_inf).^2;
figure; plot(th(1:N+1)*r2d,[Cp;Cp(1)],'d'); hold on; 
plot(t*r2d,1-4*sin(t).^2)
legend('Panel','Potential');
xlabel('\theta [deg.]'); ylabel('Cp'); grid on;

Ns=100;
Ng_y=-2*R:R/Ns:2*R;
Ng_x=-4*R:R/Ns:4*R;
[Xs, Ys]=meshgrid(Ng_x,Ng_y); [Ng_i,Ng_j]=size(Xs);

Us=zeros(size(Xs)); Vs=zeros(size(Xs)); uL=Us; vL=Us;

j=1:N; Sj=sqrt((sx(j+1)-sx(j)).^2+(sy(j+1)-sy(j)).^2);

for i=1:Ng_i
    for j=1:Ng_j
        for k=1:N
            xL= cos(tcp_ang(k))*(Xs(i,j)-cx(k))+sin(tcp_ang(k))*(Ys(i,j)-cy(k));
            yL=-sin(tcp_ang(k))*(Xs(i,j)-cx(k))+cos(tcp_ang(k))*(Ys(i,j)-cy(k));
            uL=Lambda(k)/(2.*pi)*0.5*log( ((xL+0.5*Sj(k))^2+yL^2) ...
                /((xL-0.5*Sj(k))^2+yL^2));
            vL=Lambda(k)/(2.*pi)*(atan( (xL+0.5*Sj(k))/yL) ...
                - atan( (xL-0.5*Sj(k))/yL));
            Us(i,j)=uL*cos(tcp_ang(k))-vL.*sin(tcp_ang(k)) +Us(i,j);
            Vs(i,j)=uL*sin(tcp_ang(k))+vL.*cos(tcp_ang(k))+Vs(i,j);

        end
    end
end

Us = Us+windvel_x;
Vs = Vs+windvel_y;

figure(3);
Sx = Xs(35:4:Ng_i-34,1); Sy=Ys(35:4:Ng_i-34,1);
h = streamline(Xs,Ys,Us,Vs,Sx,Sy);
set(h,'Color','red'); hold on; plot(sx,sy,'Linewidth',2);
plot(x,y,'-');
view(2);
axis equal;

figure(4);
Cp=1.-(Us.^2+Vs.^2)/(windvel_x^2+windvel_y^2);
contour(Xs,Ys,Cp); colorbar;
title('Cp');
axis equal;