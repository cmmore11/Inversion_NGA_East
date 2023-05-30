%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function dat = PSO(problem, params)

    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function

    nVar = problem.nVar;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables


    %% Parameters of PSO

    MaxIt = params.MaxIt;   % Maximum Number of Iterations

    nPop = params.nPop;     % Population Size (Swarm Size)

    w = params.w;           % Intertia Coefficient
    wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
    c1 = params.c1;         % Personal Acceleration Coefficient
    c2 = params.c2;         % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = inf;

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        [particle(i).Cost,particle(i).PSA,particle(i).Exceed] = CostFunction(particle(i).Position);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;
        particle(i).Best.PSA = particle(i).PSA;
        particle(i).Best.Exceed = particle(i).Exceed; 

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, 1);
    BestPops = zeros(MaxIt,nVar);
    BestExceeds = zeros(MaxIt,1);


    %% Main Loop of PSO
%% Main Loop of PSO
    it = 1;
    same = 0;
    while (it<MaxIt) && (same <50)

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position % turncate two decimals.
            particle(i).Position = fix((particle(i).Position + particle(i).Velocity)*10^2)/(10^2);
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = fix((max(particle(i).Position, VarMin))*10^2)/(10^2);
            particle(i).Position = fix((min(particle(i).Position, VarMax))*10^2)/(10^2);

            % Evaluation
            [particle(i).Cost,particle(i).PSA,particle(i).Exceed] = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;
                particle(i).Best.PSA = particle(i).PSA;
                particle(i).Best.Exceed = particle(i).Exceed; 

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                    
                else
                    
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(it) = GlobalBest.Cost;
        BestPops(it,:) = GlobalBest.Position;
        BestExceeds(it) = GlobalBest.Exceed;
        if it>1
            if round(GlobalBest.Cost,4) == round(BestCosts(it-1),4)
                same = same+1;
            else
                same = 0;
            end
        end

        % Display Iteration Information
        if ShowIterInfo
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
            disp(['   xmin =' mat2str(BestPops(it,:))]);
            disp(['   exceedances =' num2str(BestExceeds(it))]);
        end

        % Damping Inertia Coefficient
        w = w * (1-wdamp);
        it = it+1;

    end

    dat.PSA = GlobalBest.PSA;
    dat.xmin = GlobalBest.Position;
    dat.fxmin = GlobalBest.Cost;
    dat.xmingen = BestPops;
    dat.fxmingen = BestCosts;
        
end