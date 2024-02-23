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
y_predator = zeros(size(t_val));
y_corrector = zeros(size(t_val));
y_predator(1) = y0;
y_corrector(1) = y0;

for i = 1 : length(t_val) - 1
    % getting the next value by using RK-4
    k1 = step * dydt(t_val(i), y_predator(i));
    k2 = step * dydt(t_val(i) + step * 0.5, y_predator(i) + k1 * 0.5);
    k3 = step * dydt(t_val(i) + step * 0.5, y_predator(i) + k2 * 0.5);
    k4 = step * dydt(t_val(i + 1), y_predator(i) + k3);
    y_predator(i + 1) = y_predator(i) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    y_corrector(i + 1) = y_predator(i + 1);
    
    % predicting the values
    if (i > 3)
        current = dydt(t_val(i), y_predator(i));
        one_step_back = dydt(t_val(i - 1), y_predator(i - 1));
        two_step_back = dydt(t_val(i - 2), y_predator(i - 2));
        three_step_back = dydt(t_val(i - 3), y_predator(i - 3));
        y_predator(i + 1) = y_predator(i) + step / 24 * (55 * current - 59 * one_step_back + 37 * two_step_back - 9 * three_step_back);
    end
    
    % updating the values
    if(i > 3)
        one_step_ahead = y_predator(i + 1);
        current = dydt(t_val(i), y_predator(i));
        one_step_back = dydt(t_val(i - 1), y_predator(i - 1));
        two_step_back = dydt(t_val(i - 2), y_predator(i - 2));
        y_corrector(i + 1) = y_predator(i) + step / 24 * (9 * one_step_ahead + 19 * current - 5 * one_step_back + two_step_back);
    end
end

% printing the approximations
for i = 1 : length(t_val)
    exact = subs(sol, x, t_val(i));
    err = abs(exact - y_corrector(i));
    fprintf('Predicted y: %0.7f\tCorrected y: %0.7f\texact: %0.7f\tError: %0.7f\n', y_predator(i), y_corrector(i), exact, err);
end

% plotting the approximations
plot(t_val, y_predator, 'b-x');
hold on;
plot(t_val, y_corrector, 'g-o');
fplot(sol, [0, 2], 'r'); % for exact solution
title('Predator-Corrector using ABM & AMM');
xlabel('t');
ylabel('y');
legend('predator', 'corrector', 'exact', 'Location', 'Northwest');
hold off;
