%double integrator lqr
%x_doubledot = u

close all
clear all
%initialise
t_f = 5;
dt = 0.001;
P_f= 100*eye(2);
P0 =P_f;

%plant dynmaics
A = [0 1;0 0];
B = [0 ;1];

%initial state
X0=[10;0];
X_eval(:,1) = X0;

%cost function
Q = eye(2);
R = 1;

%controller
t_d = 0:dt:t_f;
t_r = t_f:-dt:0;

%evaluate all p
P_all_history(:,length(t_r))= P_f(:);
for i = length(t_r):-1:2
    P =reshape(P_all_history(:,i),size(A));
    dPdt = -(A.'*P + P*A - P*B*B.'*P + Q); 
    P = P - dt*(dPdt);
    P_all_history(:,i-1)= P(:);
end
P_all_history = P_all_history';

%calculate state at net step using p
for i = 2:length(t_d)
    P_eval = reshape(P_all_history(i-1,:),size(A));
    U_eval(:,i-1) = -inv(R)*B'*P_eval*X_eval(:,i-1);
    X_eval(:,i) = X_eval(:,i-1) + dt* (A*X_eval(:,i-1) + B*U_eval(:,i-1) );  
end

figure
plot(t_d(1:end-1),U_eval)
ylabel('control input')
xlabel('time')

figure;
plot(t_d,X_eval)
ylabel('states')
xlabel('time')
legend({'position','velocity'})