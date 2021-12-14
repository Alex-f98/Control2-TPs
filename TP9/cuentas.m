syms A B C D L K Ka s Ts

A_= [A B*K B*Ka; L*C (A -L*C + B*K) B*Ka; -C*Ts 0 0]
B_= [0; 0; Ts]
C_= [C 0 0]
D_= 0

I=eye(3)

H= C_ * inv([s*I-A_])*B_ + D_
%Lim S*Y= S*H*R= S*H*(1/S)=H(s)= 1
LIM= s*H*(1/s)
LIM=H %(s->0): -(B*C*Ka/-B*C*Ka= 1... no importa el Ka? claro. pero sirve de algo?

%%
clear all; close all; clc
A = [0 1; 0 -5]
B = [0; 50]
C = [1 0]
D = [0]

K= [-1.74213821 -0.17366766]
L= [0.57463229; 6.55795836]

zeta = 0.707
wn = 10
zeta_obs = 0.5
wn_obs = 20
%lambd = np.roots([1 2*zeta*wn wn**2])
%lambd_obs = np.roots([1 2*zeta_obs*wn_obs wn_obs**2])

Ts = 2*pi/(max(wn, wn_obs))/10 % Por como es la respuesta oscilatoria subamortiguada, llega a representar al menos 10 puntos por período
Ts = round(Ts/10e-3)*10e-3 % Redondea a decenas de ms

motorcc= ss(A,B,C,D)
motorcc_d = c2d(motorcc, Ts) 

A_t= [motorcc_d.A  motorcc_d.B*K; L*motorcc_d.C (motorcc_d.A-L*motorcc_d.C+motorcc_d.B*K)]
B_t= [motorcc_d.B; motorcc_d.B]
C_t= [C zeros( 1, length(A) )]
D_t=zeros(1,1)

sys_d= ss(A_t, B_t, C_t, D_t, Ts)

theta0 = pi
omega0 = 10
theta0_o = 0
omega0_o = 0
x0 = [theta0; omega0; theta0_o; omega0_o]

t_ini = 0
t_final = 1
t_step = Ts
%t = linspace(t_ini, t_final, t_step)
t= t_ini:t_step:t_final
[y,tiempo,x]= initial(sys_d, x0, t)

figure()
plot( tiempo, x(:,1),'.-')
hold on
plot( tiempo, x(:,3),'.-')
legend({'theta', 'theta_o'})
grid
title('Comparación de evolución de x_1 y X*_1, K=0')

figure()
plot( tiempo, x(:,2),'.-')
hold on
plot( tiempo, x(:,4),'.-')
legend({'omega', 'omega_o'})
grid
title('Comparación de evolución de x_2 y X*_2, K=0')
%%
%VEAMOS AHORA LA ACCION INTEGRAL
Ka= 1

Aa = [  [ sys_d.A  [motorcc_d.B*Ka; motorcc_d.B*Ka]  ] ;  [ -motorcc_d.C zeros(1,  length(B*K)+ 1 ) ] ]
Ba = [ [zeros(length(sys_d.A), 1 )]; [ones(1, 1)] ]   
Ca = [sys_d.C zeros(1,1)]
Da = zeros(1,1)

sys_rbt = ss(Aa, Ba, Ca, Da, Ts)

x0 = [theta0, omega0, theta0_o, omega0_o, 10]
t= 0:Ts:4
escalon= ones(length(t),1);

%[salida,tiempo,x]= initial(sys_rbt, x0, t)
[Y,tiempo,xout]= lsim(sys_rbt, escalon, t, x0)

figure()
plot( tiempo, escalon,'.-')
hold on
plot( tiempo, Y,'.-')
legend({'escalon', 'rta-escalon'})
grid
title('respuesta al escalon unitario')

%%
Ka=1
lambd = roots([1 2*zeta*wn wn^2])
lambd_obs = roots([1 2*zeta_obs*wn_obs wn_obs^2])

K= -place(A, B,[lambd])
L= place(A', C', [lambd_obs])'

Aa= [ A+B*K -B*K B*Ka; L*C A-L*C+B*K B*Ka; -C zeros(1,  length(B*K)+ 1 )]
Ba = [ [zeros(2*length(A), 1 )]; [ones(1, 1)] ]   
Ca = [C zeros(1, 3)]
Da = zeros(1,1)

sys_rbt = ss(Aa, Ba, Ca, Da)
sys_rbt_d = c2d(sys_rbt, Ts)

x0 = [theta0, omega0, theta0_o, omega0_o, 10]
t= 0:Ts:4
escalon= pi*ones(length(t),1);

%[salida,tiempo,x]= initial(sys_rbt, x0, t)
[Y,tiempo,xout]= lsim(sys_rbt_d, escalon, t, x0)

figure()
plot( tiempo, escalon,'.-')
hold on
plot( tiempo, Y,'.-')
legend({'escalon', 'rta-escalon'})
grid
title('respuesta al escalon unitario')


%%
%Accion de control

X_casa= xout(:,[2, 3])
Xa= xout(:, 4)

U_accion_control= K(1)*X_casa(1) + K(2)*X_casa(2) + Xa*Ka

figure()
plot( tiempo, U_accion_control,'.-')
legend({'accion de conttrol'})
grid
title('U(t)')

