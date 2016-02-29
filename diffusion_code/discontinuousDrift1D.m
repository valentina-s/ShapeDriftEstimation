T = 100000;
dt = 0.1;
N = round(T/dt);
theta = 1;

% generating the diffusion
x = zeros(N,1);
x(1) =5;
for k = 1:N-1
    x(k+1,1) = x(k,1)-theta*sign(x(k))*dt + sqrt(dt)*randn;
    %hist(x(1:k+1,1),100)
    %pause(0.001)
end

% calculating the likelihood ratio
dx = circshift(x,[-1]) - x;

theta = -(sign(x)'*dx)/T






