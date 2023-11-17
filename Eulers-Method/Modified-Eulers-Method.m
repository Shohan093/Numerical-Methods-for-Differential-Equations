clear all
% using clear all at the beginning for avoiding
% conflict among variables

% initializing symbols
syms exact_sol t;

% defining the ODE function y' = y - t^2 + 1
dydt = @(t, y) y - t^2 + 1;

% defining the initial conditions
t0 = 0;
y0 = 0.5;

% defining the final value of x
t_final = 2;

% step size
step = 0.1;

% initialinzing the arrays to store the variables
t_val = t0 : step : t_final;
y_val = zeros(size(t_val));
y_val(1) = y0;

% applying Modified Euler's method
for i = 1 : length(t_val) - 1
    % estimating slope at the next point using Euler's Method
    slope = dydt(t_val(i), y_val(i));
    % using the estimated slop to predict the next point
    y_estimate = y_val(i) + step * slope;
    % calculate the slope at the next point based on estimate
    slope_estimate = dydt(t_val(i) + step, y_estimate);
    % updating y_val using average slope
    y_val(i + 1) = y_val(i) + 0.5 * step * (slope + slope_estimate);
end

% exact solution
exact_sol = (t + 1)^2 - 0.5 * exp(t);

% displaying the solutions
for i = 1 : length(t_val)
    exact = subs(exact_sol, t_val(i));
    fprintf("t: %0.2f \ty: %0.7f\t exact: %0.7f\n", t_val(i), y_val(i), exact);
end

% plotting both the solutions
% numeric solution
plot(t_val, y_val, 'b-o');
hold on;
% symbolic solution
fplot(exact_sol, [0 2], 'r');
title('Numerical solution using Euler''s method');
xlabel('t');
ylabel('y');
legend('numeric', 'analytic');
