clc
close all
clear all 

%define time and steps
t_f = 4;
dt = 0.02;
t_dis = 0:dt:t_f;
N = length(t_dis);

%define cost function
S{N} = 100*eye(2);
R = dt;
Q = .0001*eye(2)*dt;
K{N} = [0 0];
K_norm(N)=0;

%define plant dynamics
A = [1 0.02;0 1];
B = [0 ;0.02];

for i = N-1:-1:1
    K{i} = inv(R + B'*S{i+1}*B)*B'*S{i+1}*A;
    S{i} = Q + K{i}'*R*K{i} + (A-B*K{i})'*S{i+1}*(A-B*K{i});
    K_norm(i) = norm(K{i});
end

X(:,1) = [1;0];
X_dlqr(:,1) = [1;0];
P_dlqr = dare(A,B,Q,R);
K_dlqr = inv(R)*B'*P_dlqr;

for i = 1:N-1
    u(i) = -K{i}*X(:,i);
    X(:,i+1) = A * X(:,i) + B*u(i);
    X_dlqr(:,i+1) = A * X_dlqr(:,i) - B*K_dlqr*X_dlqr(:,i) ;

end

figure;
subplot(2,1,1)
plot(t_dis,X(1,:),t_dis,X_dlqr(1,:))
legend('Dynamic Programming','LQR')
ylabel('position')
subplot(2,1,2)
plot(t_dis,X(2,:),t_dis,X_dlqr(2,:))
legend('Dynamic Programming','LQR')
ylabel('velocity')
xlabel('time')