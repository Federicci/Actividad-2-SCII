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
%Tl=0;

A=[-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B=[1/Laa; 0; 0];
C=[0 0 1];
D=[0];

M_cont=[B A*B A*A*B];
rank(M_cont) %->3, el sistema es controlable

%Diseño con LQR
Q=diag([1 1/10000 1/40]);    R=0.1;
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

%G para referencia distinta de 0, corrimiento de origen
G=-inv(C*inv(A-B*K)*B);

%Simulación del control:
deltat=1e-5;
ts=1.2;
pasos=round(ts/deltat);
Ci=[0 0 0];
t=0:deltat:(ts-deltat);
%Funciones de referencia y torque mientras va variando el tiempo:
ref=(pi/2)*square(2*pi*t/0.6);
fTl=(Tl/2)*square(2*pi*t/0.6)+Tl/2;

x=zeros(3,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
u(1)=0;

for i=2:1:pasos
    x_actual=[x(1,i-1); x(2,i-1); x(3,i-1)];
    u_actual=-K*x_actual+ref(i-1)*G;
    u_actual=-K*x_actual+10*(ref(i-1)-C*x_actual);
    u=[u u_actual];
    
    x1_p=-Ra*x_actual(1)/Laa-Km*x_actual(2)/Laa+u_actual/Laa;
    x2_p=Ki*x_actual(1)/J-Bm*x_actual(2)/J-fTl(i-1)/J;
    x3_p=x_actual(2);
    x_p_actual=[x1_p; x2_p; x3_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(1,i)=x_sig(1);
    x(2,i)=x_sig(2);
    x(3,i)=x_sig(3);
end

figure
subplot(1,2,1)
hold on;
grid on;
plot(t,ref,'k');
title('Referencia');
xlabel('Tiempo');
ylabel('Ángulo');
subplot(1,2,2)
hold on;
grid on;
plot(t,fTl,'k');
title('Torque de perturbación');
xlabel('Tiempo');
ylabel('Torque');

figure
grid on;
hold on;
plot(t,x(3,:),'r');
plot(t,ref,'k');
legend({'Salida','Referencia'},'Location','southeast');
title('Salida del sistema');
xlabel('Tiempo');
ylabel('Ángulo');
