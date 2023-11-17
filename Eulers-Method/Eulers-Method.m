clear all
% using clear all at the beginning to avoid conflict among variables

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

% initializing the arrays to store the variables
t_val = t0 : step : t_final;
y_val = zeros(size(t_val));
y_val(1) = y0;

% applying Euler's method
for i = 1 : length(t_val) - 1
    y_val(i + 1) = y_val(i) + step * dydt(t_val(i), y_val(i));
end

% exact solution
exact_sol = (t + 1)^2 - 0.5 * exp(t);

% displaying the solutions
for i = 1 : length(t_val)
    exact = subs(exact_sol, t_val(i));
    fprintf("t: %0.2f \ty: %0.7f\t exact: %0.7f\n", t_val(i), y_val(i), exact);
end

% plotting both the solutions
% numerical plot
% plotting both the solutions
plot(t_val, y_val, 'b-o');
hold on;
% symbolic plot
fplot(exact_sol, [0 2], 'r');
title('Numerical solution using Euler''s method');
xlabel('t');
ylabel('y');
legend('numeric', 'analytic');
