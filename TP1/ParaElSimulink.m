%TP1
clc; clear all; close all;

w= 90 %rad/s

M = 0.3 % kg
K = 0.03 % N/m
B = 0.3 % Ns/m
Ma = 0.03  % kg

Ka= w^2*Ma %debe ser mayor a cero.

%[1 (Ka/Ma)]
%[1 (B/M) (Ka/M + K/M + Ka/M) (B*Ka)/(M*Ma) (K*Ka)/(M*Ma)]*M