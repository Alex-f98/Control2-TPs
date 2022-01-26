clc; clear all; close all;

A= [1 0; 2 0]; B= [1; 0];D= 0;
C1= [0 1]; C2= [1 0];
S1= [0 0; 0 0];
S2= [2 0; 0 2];

P_ctrb= ctrb(A,B); rank(P_ctrb) %ES CONTROLABLE->Estabilizable.

Q1= obsv(A, C1); rank(Q1) % ES OBSERVABLE
Q2= obsv(A, C2); rank(Q2) % NO ES OBSERVABLE
%Matriz de ponderacion.
Q= C1'*C1
%Q= C2'*C2
R= 1
%Resolvemos la ecuacion dif. de ricatti.
P= care(A,B,Q,R)


%u(t)= {-inv(R)*B'*P} X(t)= K.X(t)
K= -inv(R)*B'*P

%Los autovalores de la realimentacion a lazo cerrado seran:
eig(A + B*K)  %sqrt(5)/2 +- j*sqrt(3)/2
%Lo cual indican que son asintoticamente estables.
%%  A'*P + P*A + C'*C - P*B*inv(R)*B'*P = 0
syms P11 P12 P12 P22 
A= [0 1; 0 0]; B= [0; 1];   
C1= [1 0]; C2= [0 1];
C=C2
R= 1

P= [P11 P12; P12 P22]
A'*P + P*A + C'*C - P*B*inv(R)*B'*P
