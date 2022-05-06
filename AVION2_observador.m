clc, clear all, close all;

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

w=3;
a=0.05;
b=5;
c=100;

A=[-a a 0 0; 0 0 1 0; w^2 -w^2 0 0; c 0 0 0];
B=[0; 0; b*w^2; 0];
C=[0 0 0 1];
%x1=alpha x2=phi x3=phi_p x4=h

%Calculo del controlador. Diseño por Ackerman, posicionamiento de polos
p1=-15+15i;
p2=-15-15i;
p3=-0.5+0.5i;
p4=-0.5-0.5i;

%syms p1 p2 p3 p4 s
%expand((s-p1)*(s-p2)*(s-p3)*(s-p4))
a0=1;
a1=-p1-p2-p3-p4;
a2=p1*p2+p1*p3+p1*p4+p2*p3+p2*p4+p3*p4;
a3=-p1*p2*p3-p1*p2*p4-p1*p3*p4-p2*p3*p4;
a4=p1*p2*p3*p4;
phi_A=a0*A^4+a1*A^3+a2*A^2+a1*A^1+a0*A^0;

AUX=[B A*B A*A*B A*A*A*B];
K=[0 0 0 1]*inv(AUX)*phi_A;
K=place(A,B,[p1;p2;p3;p4]);

%G para referencia distinta de 0:
G=-inv(C*inv(A-B*K)*B);

%Cálculo del observador
A_o=A';
B_o=C';
C_o=B';

K_o=place(A_o,B_o,[-5 -3 -1-i -1+i]);

%Simulacion
deltat=10^-3;
ts=70;
pasos=ts/deltat;
Ci=[0 0 0 500];
t=0:deltat:(ts-deltat);
ref=-100;
x=zeros(4,pasos);
x(1,1)=Ci(1);
x(2,1)=Ci(2);
x(3,1)=Ci(3);
x(4,1)=Ci(4);
x_hat=zeros(4,pasos);
x_hat(1,1)=0;
x_hat(2,1)=0;
x_hat(3,1)=0;
x_hat(4,1)=0;
u(1)=0;

for i=2:1:pasos
    x_actual=x(:,i-1);
    x_hat_actual=x_hat(:,i-1);
    u_actual=-K*x_hat_actual+0*ref*G;
    u=[u u_actual];
    
    x1_p=-a*x_actual(1)+a*x_actual(2);
    x2_p=x_actual(3);
    x3_p=w^2*x_actual(1)-w^2*x_actual(2)+b*w^2*u_actual;
    x4_p=c*x_actual(1);
    x_p_actual=[x1_p; x2_p; x3_p; x4_p];
    
    x_sig=x_actual+deltat*x_p_actual;
    x(1,i)=x_sig(1);
    x(2,i)=x_sig(2);
    x(3,i)=x_sig(3);
    x(4,i)=x_sig(4);
    
    y_actual=C*x_actual;
    y_hat_actual=C*x_hat_actual;
    e=y_actual-y_hat_actual;
    
    x_hat_p=K_o*e+A*x_hat_actual+B*-K*x_hat_actual;
    
    x_hat_sig=x_hat_actual+deltat*x_hat_p;
    x_hat(1,i)=x_hat_sig(1);
    x_hat(2,i)=x_hat_sig(2);
    x_hat(3,i)=x_hat_sig(3);
end

figure
plot(t,x(4,:));
figure
plot(t,u);
