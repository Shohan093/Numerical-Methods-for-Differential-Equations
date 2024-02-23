clear all

% exact solution
syms x;
sol = (x + 1)^2 - exp(x) * 0.5;

% defining the ODE y' = y - t^2 + 1
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

% for AMM Four-Step method: 
%   - find the necessary approximation by using 
%     other methods (i.e RK-4, Euler's)
% apply Adams-Moulton Four-step to update the values
for i = 1 : length(t_val) - 1
    % getting the next value by using RK-4
    k1 = step * dydt(t_val(i), y_val(i));
    k2 = step * dydt(t_val(i) + step * 0.5, y_val(i) + k1 * 0.5);
    k3 = step * dydt(t_val(i) + step * 0.5, y_val(i) + k2 * 0.5);
    k4 = step * dydt(t_val(i + 1), y_val(i) + k3);
    y_val(i + 1) = y_val(i) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    
    % applying the Adams-Moulton of Four-Step to update
    if(i > 4)
        one_step_ahead = dydt(t_val(i + 1), y_val(i + 1));
        current = dydt(t_val(i), y_val(i));
        one_step_back = dydt(t_val(i - 1), y_val(i - 1));
        two_step_back = dydt(t_val(i - 2), y_val(i - 2));
        three_step_back = dydt(t_val(i - 3), y_val(i - 3));
        y_val(i + 1) = y_val(i) + step / 720 * (251 * one_step_ahead + 646 * current - 264 * one_step_back + 106 * two_step_back - 19 * three_step_back);
    end
end

% printing the approximations
for i = 1 : length(t_val)
    exact = subs(sol, x, t_val(i));
    err = abs(exact - y_val(i));
    fprintf("t: %0.2f\ty: %0.07f\texact: %0.7f\terror: %0.7f\n", t_val(i), y_val(i), exact, err);
end

% plotting both the solutions
plot(t_val, y_val, 'b-o');
xlabel('t');
ylabel('y');
hold on;
fplot(sol, [0, 2], 'r'); % for exact solution
title('Numerical solution using Adams-Moulton Four-step');
legend('numeric', 'analytic');
hold off;
