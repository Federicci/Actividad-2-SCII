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

%Dise?o con LQR
Q=1*diag([1000 1 1000 10]);    R=.1;
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
deltat=10^-3;
ts=50;
pasos=round(ts/deltat);
Ci=[0 0 pi 0];
t=0:deltat:(ts-deltat);

%Funciones de referencia y masa mientras va variando el tiempo:
ref_dist=5*square(2*pi*t/ts)+5;
m=ones(1,pasos);
m=m*0.1;
m((pasos/2):end)=m((pasos/2):end)*10;

x=zeros(4,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
ua(1)=0;

x_OP=[0;0;pi;0];

for i=2:1:pasos
    x_actual=x(:,i-1);
    u_actual=-K*x_actual+ref_dist(i-1)*G;
    ua=[ua u_actual];
    
    x_1_p=x_actual(2);
    x_2_p=-Fricc*x_actual(2)/M-m(i-1)*g*(x_actual(3)-x_OP(3))/M+u_actual/M;
    x_3_p=x_actual(4);
    x_4_p=-Fricc*x_actual(2)/(l*M)-g*(m(i-1)+M)*(x_actual(3)-x_OP(3))/(l*M)+u_actual/(l*M);
    x_p_actual=[x_1_p;x_2_p;x_3_p;x_4_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(:,i)=x_sig;
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
plot(x(3,:),x(1,:));
figure
plot(x(1,:),x(2,:));