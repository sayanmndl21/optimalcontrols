%double integrator lqr
%x_doubledot = u

close all
clear all
%initialise
t_f = 10;
dt = 0.001;
P_f= 100*eye(2);
P0 =P_f;

%plant dynmaics
A = [0 1;0 0];
B = [0 ;1];

%initial state
x0=[10;0];
X_eval(:,1) = x0;

%cost function
Q = eye(2);
R = 1;

%gain threshold
K_high = [200 40];
K_low = [30 .4];

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

%check controllability
M_a = rank(ctrb(A,B));
t = 0:0.01:10;

%care to calculate P in ricatti eqn
%can use dare
P = care(A,B,Q,R);
K = inv(R)*B'*P;

%calculate for optimum K from above
K_opt = K
sys_fun = @(t,x)[x(2); -K_opt*x];
[t X_opt] = ode45(sys_fun,t,x0);

%calculate for high gain
K = K_high;
sys_fun = @(t,x)[x(2); -K_high*x];
[t X_high] = ode45(sys_fun,t,x0);

%calculate for low gain
K = K_low;
sys_fun = @(t,x)[x(2); -K_low*x];
[t X_low] = ode45(sys_fun,t,x0);

figure;
plot(t,-K_high*X_high',t,-K_low*X_low',t,-K_opt*X_opt',t_d(1:end-1),U_eval)
legend('Control @ High gain','Control @ Low gain','Control @ Optimum gain','control @ Fixed time')
figure;
plot(t,X_high(:,1),t,X_low(:,1),t,X_opt(:,1),t_d,X_eval(1,:))
legend('Position @ High gain','Position @ Low gain','Position @ Optimum gain','Position @ Fixed time')
