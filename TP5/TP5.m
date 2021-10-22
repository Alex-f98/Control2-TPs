clear all;
close all;
clc
s=tf('s')



R = 20  % Ohm
C = 5e-6  % F
L = 10e-6  % H
% Ts = ... Completá con el valor que creas adecuado
A_ = [[0  1]; [-1/(L*C)  -2/(R*C)]]
B_ = [[0]; [1]]
C_ = [[-1/(L*C) -1/(R*C)]]
D_ = [1]

numerador = [1 (1/(R*C))];
denominador = [1 (2/(R*C)) (1/(L*C))];
SIST = tf(numerador,denominador);

SS_C= ss(A_,B_,C_,D_)
figure()
step(SS_C)

figure()
step(SIST)

wd= 141067.324351708  %frecuencia natural amortiguada.(la parte imaginaria de los polos del sistema)
T= 2*pi/(wd)
Ts=T/4

[V,D]= eig(A_)

Ad_ZOH= expm(A_*Ts)
Bd_ZOH= inv(A_)*( Ad_ZOH - eye(2) )*B_