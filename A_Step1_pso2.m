%
%OUTPUTS:Res_M4_i.mat- mat file with 
%
%INPUTS:
%
%
%

% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com

clc;
clear;
close all;

%% Problem Definiton
problem.CostFunction = @(x) vec_fit_alt(x);  % Cost Function (DEFINE MAGNITUDE IN VEC_FIT_ALT.m)
problem.VarMin = [ 0.1,    0, -1.3*ones(1,8),    -1, 3.5,0.25]; % Lower Bound of Decision Variables
problem.VarMax = [ 6.0, 3.00, -0.6*ones(1,8),   0.5, 7.5,0.8];   % Upper Bound of Decision Variables
problem.nVar = length(problem.VarMin(1,:));      % Number of Unknown (Decision) Variables

%% Parameters of PSO
% Constriction Coefficients: strategy from Zhang et al. (2014). See
% inversion paper for more info and citation.
    popsss =120;
    mu = 0.675;
    tp = pi*mu;
    damp = sqrt((log(mu)^2)/((log(mu)^2)+(pi^2)));
    chi = exp((((log(mu^2)))/tp));
    term = (tan(((pi/tp)-pi))^2);
    phi1 = -((chi+chi*term-2*sqrt(chi*(term+1))+term+1)/(-term-1));
    phi2 = phi1;
    kappa = 1;
 
params.MaxIt = 300;    % Maximum Number of Iterations
params.nPop = popsss;    % Population Size (Swarm Size)
params.w = chi;        % Inertia Coefficient
params.wdamp = damp;     % Damping Ratio of Inertia Coefficient
params.c1 = phi1;      % CP, Personal Acceleration Coefficient
params.c2 = phi2;        % CG, Social Acceleration Coefficient
params.ShowIterInfo = true; % Flag for Showing Iteration Informatin

%% Calling PSO
for i = 1:10
dat = PSO(problem, params);
filename = sprintf('Res_M4_i%d.mat',i+50);

save(filename,'dat');
%% Results
end

%%%%OLD/ORIGINAL from Yarpiz
% % params.w = chi;             % Inertia Coefficient
% % %params.wdamp = 1;           % Damping Ratio of Inertia Coefficient
% % params.c1 = chi*phi1;       % CP, Personal Acceleration Coefficient
% % params.c2 = chi*phi2;       % CG, Social Acceleration Coefficient

