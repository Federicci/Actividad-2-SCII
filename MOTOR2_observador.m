clc, clear all, close all;

%Motor con carga
%Variables de estado: x1=ia, x2=wr, x3=titat

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

Laa=366e-6;
J=5e-9;
Ra=55.6;
Bm=0;
Ki=6.49e-3;
Km=6.53e-3;

Tl=1.15e-3; %Solo para la referencia de pi/2

A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];
D=[0];

%Agrego un integrador para mejorar el error en estado estacionario por la
%perturbacion de torque:
%Amplio el sistema

AA=[A zeros(3,1); -C 0];
BB=[B; 0];
CC=[C 0];

%Diseño con LQR
QQ=diag([1 1/10000 1/40 1000000]);    RR=0.1;
%Hamiltoniano
HH=[AA -BB*inv(RR)*BB'; -QQ -AA'];
[vects1,autovals1]=eig(HH);  %columnas de vects: autovectores
%Debo extraer solo los autovectores cuyos autovalores son negativos:
autovects_neg1=[];
for i=1:1:length(autovals1)
    if (real(autovals1(i,i)))<0
        autovects_neg1=[autovects_neg1 vects1(:,i)];
    end
end    

%divido la matriz de autovectores en 2 matrices:
[filas1,colums1]=size(autovects_neg1);
M1=autovects_neg1(1:(filas1/2),:);
PM1=autovects_neg1((filas1/2+1):filas1,:);
P1=real(PM1*inv(M1));

%Con la matriz P construyo el controlador
KK=inv(RR)*BB'*P1;
%KK=[K -Ki]

%Calculo del observador
A_o=AA';
B_o=CC';
C_o=BB';

M_cont=[B_o A_o*B_o A_o*A_o*B_o A_o*A_o*A_o*B_o];
rank(M_cont)  %No es controlable?

Q_o=1e5*diag([1 1 1 1]);    R_o=1;
Q_o=20*diag([1 1/10000 1/40 1000000]);    R_o=0.1/20;

H_o=[A_o -B_o*inv(R_o)*B_o'; -Q_o -A_o'];

[vects_o,autovals_o]=eig(H_o);  %columnas de vects: autovectores
%Debo extraer solo los autovectores cuyos autovalores son negativos:
autovects_neg_o=[];
for i=1:1:length(autovals_o)
    if (real(autovals_o(i,i)))<0
        autovects_neg_o=[autovects_neg_o vects_o(:,i)];
    end
end    

%divido la matriz de autovectores en 2 matrices:
[filas_o,colums_o]=size(autovects_neg_o);
M_o=autovects_neg_o(1:(filas_o/2),:);
PM_o=autovects_neg_o((filas_o/2+1):filas_o,:);
P_o=real(PM_o*inv(M_o));

%Con la matriz P construyo el controlador del observador
K_o=inv(R_o)*B_o'*P_o;

%Simulación del control:
deltat=1e-5;
ts=1.2;
pasos=round(ts/deltat);
Ci=[0 0 0 0];
t=0:deltat:(ts-deltat);
%Funciones de referencia y torque mientras va variando el tiempo:
ref=(pi/2)*square(2*pi*t/0.6);
fTl=(Tl/2)*square(2*pi*t/0.6)+Tl/2;

x_hat=zeros(4,pasos);
x_hat(1,1)=0;
x_hat(2,1)=0;
x_hat(3,1)=0;
x_hat(4,1)=0;
x=zeros(4,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
u(1)=0;

for i=2:1:pasos
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    integracion=x(4,i-1)+deltat*(ref(1,i-1)-CC*x_actual);
    u_actual=-KK(1:3)*x_hat_actual(1:3)-integracion*KK(4);
    u=[u u_actual];
    
    x1_p=-Ra*x_actual(1)/Laa-Km*x_actual(2)/Laa+u_actual/Laa;
    x2_p=Ki*x_actual(1)/J-Bm*x_actual(2)/J-fTl(i-1)/J;
    x3_p=x_actual(2);
    x_p_actual=[x1_p; x2_p; x3_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(1,i)=x_sig(1);
    x(2,i)=x_sig(2);
    x(3,i)=x_sig(3);
    x(4,i)=integracion;
    
    y_actual=CC*x_actual;
    y_hat_actual=CC*x_hat_actual;
    e=y_actual-y_hat_actual;
    
    x_hat_p=e*K_o'+A*x_hat_actual+x_hat_actual*u_actual;
    
    x_hat_sig=x_hat_actual+deltat*x_hat_p;
    x_hat(1,i)=x_hat_sig(1);
    x_hat(2,i)=x_hat_sig(2);
    x_hat(3,i)=x_hat_sig(3);
end

figure
plot(t,x(3,:));
hold on;
plot(t,ref);
figure
plot(t,u);






