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

%Diseño con LQR
Q=diag([0.01 1 5]);    R=0.15;
%Hamiltoniano
H=[A -B*inv(R)*B'; -Q -A'];
[vects,autovals]=eig(H);  %columnas de vects: autovectores
%Debo extraer solo los autovectores cuyos autovalores son negativos:
autovects_neg=[];
for i=1:1:length(autovals)
    if (real(autovals(i,i)))<0
        autovects_neg=[autovects_neg vects(:,i)];
    end
end    

%divido la matriz de autovectores en 2 matrices:
[filas,colums]=size(autovects_neg);
M=autovects_neg(1:(filas/2),:);
PM=autovects_neg((filas/2+1):filas,:);
P=real(PM*inv(M));

%Con la matriz P construyo el controlador
K=inv(R)*B'*P;

%G para referencia distinta de 0:
G=-inv(C*inv(A-B*K)*B);

%Calculo del observador
A_o=A';
B_o=C';
C_o=B';
Q_o=diag([1000000 1000000 1000000]);    R_o=0.000000001;

H_o=[A_o -B_o*inv(R)*B_o'; -Q_o -A_o'];

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

%Con la matriz P construyo el controlador
K_o=inv(R_o)*B_o'*P_o;

%Simulación del control:
deltat=1e-5;
ts=0.7;
pasos=round(ts/deltat);
Ci=[0 0 0];
t=0:deltat:(ts-deltat);
%Funciones de referencia y torque mientras va variando el tiempo:
ref=(pi/2)*square(2*pi*t/0.2);
figure
plot(t,ref);
fTl=(Tl/2)*square(2*pi*t/0.2)+Tl/2;
figure
plot(t,fTl);

x_hat=zeros(3,pasos);
x_hat(1,1)=0;
x_hat(2,1)=0;
x_hat(3,1)=0;
x=zeros(3,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
u(1)=0;

for i=2:1:pasos
    u_actual=-K*x_hat(:,i-1)+ref(i-1)*G;
    u=[u u_actual];
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    
    x1_p=-Ra*x_actual(1)/Laa-Km*x_actual(2)/Laa+u_actual/Laa;
    x2_p=Ki*x_actual(1)/J-Bm*x_actual(2)/J-fTl(i-1)/J;
    x3_p=x_actual(2);
    x_p_actual=[x1_p; x2_p; x3_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(1,i)=x_sig(1);
    x(2,i)=x_sig(2);
    x(3,i)=x_sig(3);
    
    y_actual=x_actual(3);
    y_hat_actual=x_hat_actual(3);
    e=y_actual-y_hat_actual;
    
    x_hat_p=e*K_o+A*x_hat_actual+x_hat_actual*K*B;
    
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






