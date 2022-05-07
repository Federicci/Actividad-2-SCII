clear all; close all;

m=.1;
Fricc=0.1; 
l=0.6;
g=9.8;
M=.5;

%Linealizado en el equilibrio estable:
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; 1/(l*M)];
C=[1 0 0 0]; 

%Amplio el sistema
AA=[A zeros(4,1); -C 0];
BB=[B; 0];
CC=[C 0];

%Diseño con LQR
Q=1*diag([0.0001 0.001 10000 10 100]);    R=1000;
%Hamiltoniano
H=[AA -BB*inv(R)*BB'; -Q -AA'];
[vects1,autovals1]=eig(H);  %columnas de vects: autovectores
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
K=inv(R)*BB'*P1;

%Cálculo del observador
%Diseño con LQR para el observador
Q_o=1e4*diag([1 1 10 1]);    R_o=0.1;

A_o=A';
B_o=C';
C_o=B';

%Hamiltoniano
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
deltat=10^-3;
ts=50;
pasos=round(ts/deltat);
Ci=[0 0 pi 0 0];
t=0:deltat:(ts-deltat);

%Funciones de referencia y masa mientras va variando el tiempo:
ref_dist=10*square(2*pi*t/50);
m=ones(1,pasos);
m=m*0.1;
m((pasos/2):end)=m((pasos/2):end)*10;

x_hat=zeros(4,pasos);
x_hat(1,1)=0;
x_hat(2,1)=0;
x_hat(3,1)=0;
x_hat(4,1)=0;
x=zeros(5,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
x(5,1)=Ci(5);
ua(1)=0;

x_OP=[0;0;pi;0;0];

for i=2:1:pasos
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    integracion=x_actual(5)+deltat*(ref_dist(i-1)-CC*x_actual);
    u_actual=-K(1:4)*x_hat_actual(1:4)-integracion*K(5); %El - va por -Ki
    ua=[ua u_actual];
    
    x_1_p=x_actual(2);
    x_2_p=-Fricc*x_actual(2)/M-m(i-1)*g*(x_actual(3)-x_OP(3))/M+u_actual/M;
    x_3_p=x_actual(4);
    x_4_p=-Fricc*x_actual(2)/(l*M)-g*(m(i-1)+M)*(x_actual(3)-x_OP(3))/(l*M)+u_actual/(l*M);
    x_5_p=0;
    x_p_actual=[x_1_p;x_2_p;x_3_p;x_4_p;x_5_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(:,i)=x_sig;
    x(5,i)=integracion;
    
    y_actual=CC*x_actual;
    y_hat_actual=C*x_hat_actual;
    e=y_actual-y_hat_actual;
    
    x_hat_p=e*K_o'+A*x_hat_actual+B*u_actual;
    
    x_hat_sig=x_hat_actual+deltat*x_hat_p;
    x_hat(:,i)=x_hat_sig;
end

figure
plot(t,x(1,:)); 
hold on;
plot(t,ref_dist);
figure
plot(t,x(3,:));
figure
plot(t,ua); %Acc. de control.
figure
plot(x(3,:),x(4,:));
figure
plot(x(1,:),x(2,:));