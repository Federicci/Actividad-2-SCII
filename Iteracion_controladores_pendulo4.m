min_ang=10;
min_dist=20;
for p=1:1:200
    
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

%Dise?o con LQR
%Q=1*diag([0.0001 0.001 10000 10 100]);    R=1000;
%Q=1*diag([0.1 1 10 10 1]);    R=100;
Q=rand(1)*0.01*diag([0.1*rand(1) 1*rand(1) 10*rand(1) 10*rand(1) 1*rand(1)]);    R=100; Qa=Q;
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



%Controlador para cambio de masa:
m=m*10;
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
AA=[A zeros(4,1); -C 0];

%Dise?o con LQR
%Q=1*diag([0.0001 0.001 10000 10 100]);    R=1000;
Q=1*diag([0.1 1 10 10 1]);    R=100;
Q=rand(1)*1*diag([0.1*rand(1) 1*rand(1) 10*rand(1) 10*rand(1) 1*rand(1)]);    R=100; Qb=Q;
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
K1=inv(R)*BB'*P1;




%C?lculo del observador
m=.1*10;
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];

%Dise?o con LQR para el observador
Q_o=1e4*diag([1 1 10 1]);    R_o=0.1;
Q_o=rand(1)*1e4*diag([1*rand(1) 1*rand(1) 10*rand(1) 1*rand(1)]);    R_o=0.1;  Qc=Q_o;
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

%Simulaci?n del control:
deltat=10^-3;
ts=100;
pasos=round(ts/deltat);
Ci=[0 0 pi 0 0];
t=0:deltat:(ts-deltat);

%Funciones de referencia y masa mientras va variando el tiempo:
ref_dist=10*square(2*pi*t/ts);
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
x_compar=zeros(5,pasos);
x_compar(1,1)=Ci(1);
x_compar(2,1)=Ci(2);
x_compar(3,1)=Ci(3);
x_compar(4,1)=Ci(4);
x_compar(5,1)=Ci(5);

x_OP=[0;0;pi;0;0];

for i=2:1:pasos
    if m(i-1)>.5
        K=K1;
    end
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    integracion=x_actual(5)+deltat*(ref_dist(i-1)-CC*x_actual);
    u_actual=-K(1)*x_actual(1)-K(2:4)*x_hat_actual(2:4)-integracion*K(5); %El - va por -Ki
        
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
    
    x_hat_p=e*K_o'+A*(x_hat_actual-x_OP(1:4))+B*u_actual;
    
    x_hat_sig=x_hat_actual+deltat*x_hat_p;
    x_hat(:,i)=x_hat_sig;
    
    %Sistema sin observador para comparar
    x_actual_compar=x_compar(:,i-1);
    integracion_compar=x_actual_compar(5)+deltat*(ref_dist(i-1)-CC*x_actual_compar);
    u_actual_compar=-K(1:4)*x_actual_compar(1:4)-integracion_compar*K(5);
    
    x_1_p_com=x_actual_compar(2);
    x_2_p_com=-Fricc*x_actual_compar(2)/M-m(i-1)*g*(x_actual_compar(3)-x_OP(3))/M+u_actual_compar/M;
    x_3_p_com=x_actual_compar(4);
    x_4_p_com=-Fricc*x_actual_compar(2)/(l*M)-g*(m(i-1)+M)*(x_actual_compar(3)-x_OP(3))/(l*M)+u_actual_compar/(l*M);
    x_5_p_com=0;
    x_p_actual_compar=[x_1_p_com;x_2_p_com;x_3_p_com;x_4_p_com;x_5_p_com];
    
    x_sig_compar=x_actual_compar+deltat*x_p_actual_compar;
    x_compar(:,i)=x_sig_compar;
    x_compar(5,i)=integracion_compar;
end

ang_max=max(x(3,:));
dist_max=max(x(1,:));
if ang_max<min_ang && dist_max<min_dist
    Qopt1=[Qa(1,1) Qa(2,2), Qa(3,3), Qa(4,4), Qa(5,5)];
    Qopt2=[Qb(1,1) Qb(2,2), Qb(3,3), Qb(4,4), Qb(5,5)];
    Qopto=[Qc(1,1) Qc(2,2), Qc(3,3), Qc(4,4)];
    min_ang=ang_max
    min_dist=dist_max
end
end

Qopt1
Qopt2
Qopto