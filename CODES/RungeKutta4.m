function [t, Y, Extra] = RungeKutta4(odeFun, tSpan, dt, Y0)

    % This Function is the Implementation of the 4th-Order Runge-Kutta Method
    % for the Numeric Solve of Differetial Equations...
    % Input:
    %     odeFun ==> ODE Function
    %     tSpan  ==> Time Span for the ODE Solve
    %     dt     ==> Time Resolution
    %     Y0     ==> Initial Condition
    %
    % Returns:
    %     t ==> Simulation Time
    %     Y ==> ODE Numerical Solution


    StepNum  = floor((max(tSpan) - min(tSpan)) / dt) + 1;    % Number of Simulation Steps
    StateNum = numel(Y0);                                    % Number of Sys States

    % Initialize Solution and Time Vars
    Y = zeros(StateNum, StepNum);
    t = (min(tSpan) : dt : max(tSpan))';

    Y(:, 1) = Y0;

    % Main Loop
    for k = 1:StepNum - 1

        [K1, Tmp] = feval(odeFun, t(k), Y(:, k));
        K1 = dt * K1;

        Extra.CtrlSignal(:, k) = Tmp.u;
        Extra.SliSurf(:, k) = Tmp.s;
        Extra.e(:, k) = Tmp.e;
        Extra.de(:, k) = Tmp.de;

        K2 = dt * feval(odeFun, t(k) + dt/2 , Y(:, k) + K1 / 2);
        K3 = dt * feval(odeFun, t(k) + dt/2 , Y(:, k) + K2 / 2);
        K4 = dt * feval(odeFun, t(k) + dt   , Y(:, k) + K3);

        Y(:, k+1) = Y(:, k) + (K1 + K2 * 2 + K3 * 2 + K4) / 6;
    end

    % Transpose to Make this and ODE45 Look Alike
    Y = transpose(Y);
end