clc
close all
clear all 

%define time and steps
t_f = 20;
dt = 0.0001;
t_d = 0:dt:t_f;
N = length(t_d);

%define cost function
S{N} = 100*eye(2);
R = dt;
Q = .0001*eye(2)*dt;
K{N} = [0 0];

%define plant dynamics
A = [1 0.0001;0 1];
B = [0 ;0.0001];

for i = N-1:-1:1
    S{i} =  A'*inv(eye(2) + S{i+1}*B*inv(R)*B')*S{i+1}*A + Q; % Third form
    K_norm(i) = norm(K{i});
end

X(:,1) = [1;0];
X_dlqr(:,1) = [1;0];
P_dlqr = dare(A,B,Q,R);
K_dlqr = inv(R)*B'*P_dlqr;

for i = 2:N
    u(i-1) = -inv(R)*B'*S{i-1}*X(:,i-1);
    X(:,i) = A * X(:,i-1) + B*u(i-1);
    X_dlqr(:,i) = A * X_dlqr(:,i-1) - B*K_dlqr*X_dlqr(:,i-1) ;

end


figure;
subplot(2,1,1)
plot(t_d,X(1,:),t_d,X_dlqr(1,:))
legend('Dynamic Programming',' Discrete LQR')
ylabel('position')
subplot(2,1,2)
plot(t_d,X(2,:),t_d,X_dlqr(2,:))
legend('Dynamic Programming','Discrete LQR')
ylabel('velocity')
xlabel('time')