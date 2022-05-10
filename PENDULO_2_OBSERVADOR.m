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

%Diseño con LQR
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

%Diseño con LQR para el observador
Q_o=diag([1 1 1 1]);    R_o=1;
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

%Ganancia de prealimentacion para ref distinta de 0
G=-inv(C*inv(A-B*K)*B);

%Simulación del control:
deltat=10^-4;
ts=10;
pasos=round(ts/deltat);
Ci=[0 0 0.7 0];
t=0:deltat:(ts-deltat);
%Funciones de referencia y torque mientras va variando el tiempo:
ref_dist=10;
ref_ang=0;

x=zeros(4,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
x_compar=zeros(4,pasos);
x_compar(1,1)=Ci(1);
x_compar(2,1)=Ci(2);
x_compar(3,1)=Ci(3);
x_compar(4,1)=Ci(4);
x_hat=zeros(4,pasos);
x_hat(1,1)=0;
x_hat(2,1)=0;
x_hat(3,1)=0;
x_hat(4,1)=0;
ua(1)=0;
ua_compar(1)=0;

for i=2:1:pasos
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    u_actual=-K*x_hat_actual+ref_dist*G;
    ua=[ua u_actual];
    
    x_p_actual=A*x_actual+B*u_actual;
    x_sig=x_actual+deltat*x_p_actual;
    x(:,i)=x_sig;
    
    y_actual=C*x_actual;
    y_hat_actual=C*x_hat_actual;
    e=y_actual-y_hat_actual;
    
    x_hat_p=e*K_o'+A*x_hat_actual+B*u_actual;
    
    x_hat_sig=x_hat_actual+deltat*x_hat_p;
    x_hat(:,i)=x_hat_sig;
    
    
    %Sist. sin observador para comparar
    x_actual_compar=x_compar(:,i-1);
    u_actual_compar=-K*x_actual_compar+ref_dist*G;
    ua_compar=[ua_compar u_actual_compar];
    x_p_actual_compar=A*x_actual_compar+B*u_actual_compar;
    x_sig_compar=x_actual_compar+deltat*x_p_actual_compar;
    x_compar(:,i)=x_sig_compar;
end

ref_dist=ref_dist*ones(1,pasos);
figure(1)
subplot(2,2,1);
plot(t,x(1,:),'color','r');
hold on;
plot(t,x_compar(1,:),'color',[0.4660 0.6740 0.1880]);
plot(t,ref_dist,'k');
grid on;
title('Distancia');
xlabel('Tiempo');
legend({'Con observador','Sin observador','Referencia'},'Location','southeast')
subplot(2,2,2);
plot(t,x(3,:),'color','r');
hold on;
plot(t,x_compar(3,:),'color',[0.4660 0.6740 0.1880]);
grid on;
title('Ángulo');
xlabel('Tiempo');
legend({'Con observador','Sin observador'},'Location','southeast')
subplot(2,2,[3,4]);
plot(t,ua,'color','r');
hold on;
plot(t,ua_compar,'color',[0.4660 0.6740 0.1880]);
grid on;
title('Acción de control');
xlabel('Tiempo');
legend({'Con observador','Sin observador'},'Location','southeast')

%Planos de fase con scatter para variar colores
figure(2)
subplot(2,1,1);
grid on;
hold on;
plot(x(1,1),x(2,1),'b');
plot(x_compar(1,1),x_compar(2,1),'r');
%f=-t/ts+1;
f=(t-ts).^6/(ts^6);
f(end)=1;
C=zeros(pasos,3);
C(:,3)=f;
scatter(x(1,:),x(2,:),3,C,'filled')
C1=zeros(pasos,3);
C1(:,1)=f;
scatter(x_compar(1,:),x_compar(2,:),3,C1,'filled')
title('Distancia vs velocidad');
xlabel('Distancia');
ylabel('Velocidad');
legend({'Con observador','Sin observador'},'Location','southeast')

subplot(2,1,2);
grid on;
hold on;
plot(x(3,1),x(4,1),'b');
plot(x_compar(3,1),x_compar(4,1),'r');
%f=-t/ts+1;
f=(t-ts).^6/(ts^6);
f(end)=1;
C=zeros(pasos,3);
C(:,3)=f;
scatter(x(3,:),x(4,:),3,C,'filled')
C1=zeros(pasos,3);
C1(:,1)=f;
scatter(x_compar(3,:),x_compar(4,:),3,C1,'filled')
title('Ángulo vs velocidad angular');
xlabel('Ángulo');
ylabel('Velocidad angular');
legend({'Con observador','Sin observador'},'Location','southeast')