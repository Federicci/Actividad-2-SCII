clc; clear all; close all;

m=.1;
Fricc=0.1; 
l=0.6;
g=9.8;
M=.5;

%Linealizado en el equilibrio inestable:
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(l*M) g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; -1/(l*M)];
C=[1 0 0 0]; 

%Dise?o con LQR
Q=10*diag([1 1 1 1]);    R=1;
%Hamiltoniano
H=[A -B*inv(R)*B'; -Q -A'];
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
K=inv(R)*B'*P1;

%Ganancia de prealimentacion para ref distinta de 0
G=-inv(C*inv(A-B*K)*B);

%Simulaci?n del control:
deltat=10^-4;
ts=5;
pasos=round(ts/deltat);
Ci=[0 0 0.3 0];
t=0:deltat:(ts-deltat);
%Funciones de referencia y torque mientras va variando el tiempo:
ref_dist=10;
ref_ang=0;

x=zeros(4,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
ua(1)=0;

for i=2:1:pasos
    x_actual=x(:,i-1);
    u_actual=-K*x_actual+ref_dist*G;
    ua=[ua u_actual];
    
    x_p_actual=A*x_actual+B*u_actual;
    
    x_sig=x_actual+deltat*x_p_actual;
    x(:,i)=x_sig;
end

figure
plot(t,x(1,:)); 
hold on;
plot(t,ref_dist);
figure
plot(t,ua); %Acc. de control.
figure
plot(x(3,:),x(1,:));
figure
plot(x(1,:),x(2,:));