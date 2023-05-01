clc
clear
close all

%% Define Robot Parameters

Mass = [1,1,1]; % Links Masses
Length = [2,2,2]; % length of Links
LengthCOM = 1/2 * Length; % Center of Msses lengths
InertiaMoment = 1/12 * Mass .* Length.^2; % Inertia momemnts
GravConst = 9.81; % gravity constant
NumLinks = 3; % Number of Links

%% Define symbolic Parameters (joint angels and velocities)

q = sym('q',[NumLinks 1]); % joint angels
q = sym(q,'real');
dq = sym('dq',[NumLinks 1]); % joint velocities
dq = sym(dq,'real');

%% Robot Difination

Robot = ThreeDOF(Mass,Length,LengthCOM,InertiaMoment,GravConst,q,dq); % create object

%% ode Function

odeFunc = @(t,x)odeSolver(t,x,Robot); % ode Function

%% Solve Equation of Motion

dt = 0.01; % Time Step
SimTime = 5; % Simulation Time (s)
tspan = 0:dt:SimTime; % Time Span for solve eq.
PosCond = [pi/5, pi/10, -pi/7]; % Initial Positions
VelCond = zeros(1,NumLinks); % Initial Velocities
InitCond = [PosCond, VelCond]; % Initial Conditions
[T,Q,Info] = RungeKutta4(odeFunc,tspan,dt,InitCond); % Solve odeFunc
MaxIter = size(T,1); % Maximum Iteration (for plot)

%% Robot Motion Plot

fig = figure('Name','Robot','Units','normalized','Position',[0 0 1 1]);
title('Robot Behaviour')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis(5*[-1 1 -1 1 -2.5/5 1])
grid minor
View = '3D';

for i = 1:MaxIter
       P = Robot.Motion(Q(i,1:3),View,'LineWidth',4,'Color',[0, 0.4470, 0.7410],...
        'Marker','.','MarkerSize',30);
    drawnow expose
    if i ~= MaxIter
        delete(P{1});
        delete(P{2});
        delete(P{3});
    end 
end

%% Plot Results

Robot.PlotControlSignal(T(1:end-1),Info.CtrlSignal)
Robot.PlotSlidSurface(T(1:end-1),Info.SliSurf)
Robot.PlotError(T(1:end-1),Info.e,dt,'error')
Robot.PlotError(T(1:end-1),Info.de,dt,'derror')
Robot.PlotStates(T,Q)