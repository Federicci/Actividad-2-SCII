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
B=[1/Laa 0; 0 -1/J; 0 0];
C=[0 0 1];
D=[0 0];

delta=0.000001;
ts=0.7;
pasos=ts/delta;

%Agrego un integrador para mejorar el error en estado estacionario por la
%perturbacion de torque:
%Amplio el sistema

AA=[A zeros(3,1); -C 0];
BB=[B;0 0];
%Se obtendrá un sistema ampliado con controlador ampliado
MM=[BB AA*BB AA*AA*BB AA*AA*AA*BB];

%Lugar de polos deseados:
p(1)=-3; p(2)=-3; p(3)=-5; p(4)=-2;





