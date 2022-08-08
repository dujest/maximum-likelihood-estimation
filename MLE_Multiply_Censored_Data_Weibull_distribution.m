clc;
clear all;

% matrices or samples of the failure and repair times - multiply censored data
t = [353 287 189 136 34 154 286 334] % the failure times
c = [145 200 380 500 500 500 500]; % the reapir times
k = [t, c];
r = length(t);

% Sorts the failure and reapir times from the smallest to the largest one
for i=1:(length(k));
    for j=i+1:(length(k));
        if (k(j)<k(i));
            temp=k(i);
            k(i)=k(j);
            k(j)=temp;
        end
    end
end

disp('Sorted data:');
disp(k);

beta=2.1; 


y=funmcd(beta, k, t, r); % the nonlinear equation with the shape parameter (beta)

% the Newton Raphson method obtains the optimal beta
TOL=0.001;
step=0.001;
while abs(y)>TOL
    gradient = (funmcd(beta+step, k, t, r)-y)/step;
    beta=beta-y/gradient;
    y=funmcd(beta, k, t, r);
end

disp('The optimal shape parameter (beta):');
disp(beta);

theta=(sum((k.^beta)./r))^(1/beta);

disp('The optimal scale parameter (theta):');
disp(theta);

[beta1,theta1]=meshgrid((beta-beta):(beta/50):(beta+beta) , (theta-theta):(theta/50):(theta+theta));

grid_size=size(beta1);
MLE=zeros(size(beta1));

for i=1:grid_size(1);
    for j=1:grid_size(2);
        f=(((t./theta1(i,j)).^(beta1(i,j)-1)).*(beta1(i,j)/theta1(i,j))).*exp(-((t./theta1(i,j)).^(beta1(i,j))));
        R=exp(-((c./theta1(i,j)).^beta1(i,j)));
        MLE(i,j)=(prod(f))*prod(R);
    end
end

surf(beta1, theta1, MLE);
xlabel('beta');
ylabel('theta');
zlabel('MLE');
hold on;

beta_opt=round(grid_size(1)/2);
teta_opt=round(grid_size(2)/2);
MLEmax=MLE(beta_opt, teta_opt);
plot3(beta, theta, MLEmax, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'r');
hold on;

figure;
mesh(beta1, theta1, MLE);
xlabel('beta');
ylabel('theta');
zlabel('MLE');
hold on;

plot3(beta, theta, MLEmax, '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'r');
hold on;

% plots the confidence intervals for the optimal distribution parameters
figure;
[k1, h1]=contourf(beta1, theta1, (100-((MLE/MLEmax)*100)), 9);
clabel(k1, h1);
colorbar;
xlabel('beta');
ylabel('theta');
hold on;

function y = funmcd( beta, k, t, r )

y=1/beta+sum(log(t).*(1/r))-(sum((k.^beta).*log(k)))/(sum(k.^(beta)));

end