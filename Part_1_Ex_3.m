%% EnKF equation
% We will now analyze how Data Assimilation (DA) through the EnKF combines model 
% and observation information. See lecture "Introduction to Data Assimilation" 
% for more information on the EnKF equation. 
%% Setting up the DA run
% The DA run is built upon an ensemble run, as the ensemble information will 
% be used to compute the uncertainty of the model. *Please copy-paste the missing 
% lines from Exercise 2 (omit last line, |clear mP P|).*

clc, clear all, close all
addpath('data'); addpath(genpath('functions'));

% Define parameter K
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% Define timestep
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% precipitation data
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% initialize storage vector
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% Define ensemble size
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% Load gaussian noise files
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% Set noise amplitude
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% Perturb variables
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%
% clear mP P
%% 
% We will now continue setting up the Data Assimilation run. First, we need 
% to get some observations. In this case, as this bucket model is a ficticious 
% scenario, we will generate synthetic water storage observations. For that purpose, 
% we start by generating a synthetic ground truth:

% Generate synthetic ground truth
sim = [1; 24]; % simulation phase of model
S0_true = 8;   % true initial state
S_true = zeros(sim(2)-sim(1)+1,1); S_true(1)=S0_true; % initialize storage vector
[Q_true, S_true] = linearModel(S_true, P(sim(1):sim(2)), K, Dt); % run linear reservoir model
clear mP P
%% 
% We will now generate the observations, which represent a noisy measurement 
% of this synthetic truth. *[Marker for extra exercises].*_

% Generate observations
% Load observation noise level
load('data/noise_S_obs');
% Determine noise amplitude and offset
sigma_sObs = 0.3; % 30% of value (multiplicative error model)
offset = 0;
% Compute noisy observations
S_obs = S_true + offset + S_true .* sigma_sObs .* noise_S_obs;
%% 
% Now that the observations are generated, we will perturb them to then assimilate 
% them throught the EnKF.

% Matrix of Observation ensemble
Y = repmat(S_obs,1,Ne);

% Matrix of Observation error
load('data/noise_dY.mat'); % Load observation uncertainty
dY = repmat(S_true,1,Ne) .* sigma_sObs .* noise_dY;
Y_ens = Y+dY;
IN = (Y_ens<0); Y_ens(IN) = 0; clear IN; % Remove negative water storage values that might have generated during perturbation.
%% 
% Another variable we need for the DA run is the design matrix, which relates 
% the model state to the observed variable. In this case, the relation is 1-to-1.

A = [1]; % design matrix
%% Initialization run
% We will now intialize the model over a period of 12 timesteps. *Can you copy-paste 
% the missing lines from Exercise 2?*

% initialization phase of model
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% initialize storage vector for initialization phase
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% run linear reservoir model for each ensemble member
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

%% Open Loop run (Ensemble run)
% We will also run the ensemble run as done in the previous exercise, so that 
% we can plot it and copare it with the DA results. In the field of land DA, this 
% is often referred to as "Open Loop (OL)" run.
% 
% *Please copy-paste the missing lines from Exercise 2.*

% simulation phase of model
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% initialize storage vector for simulation phase
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

% run linear reservoir model for each ensemble member
%%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%

%% Data Assimilation run
% We finally launch the DA run. First, let us initialize the variables that 
% we will be using.

% Initialization
S_t = S0_ens; % Water storage states
xMinus_3d = zeros(1,Ne,sim(2)-sim(1)+1); xPlus_3d = xMinus_3d; % Model states before and after DA
Cxx_3d = zeros(1,1,sim(2)-sim(1)+1); Cxpxp_3d = Cxx_3d; Cxpxp_a_3d = Cxx_3d; % Model error structure
Sll_2d = zeros(sim(2)-sim(1)+1,1); % Observation error structure
%% 
% During DA, we will iterate over the timesteps combining two steps:

% loop over time
for t=sim(1):sim(2)
    
% Step 1: forward model run
% in this step, we will run the model forward of one timestep for all of the 
% ensemble members:

    % step 1:
    for ii=1:Ne % loop over ensemble
        % model forward integration
        [Q_tp1(1,ii), S_tp1(1,ii)] = linearModel(S_t(1,ii), P_ens(t,ii), K_ens(1,ii), Dt);
    end
    
% Step 2: DA update
% In this step, we will update the model state using the observations as well 
% as model and observation uncertainties. First, we will add the model state in 
% the vector xMinus. We can then use this information to compute the empirical 
% error covariance matrix of the model state.

    % step 2: update (ensemble Kalman filter)
    % model prediction vector
    xMinus = [S_tp1];
    % empirical model error covariance matrix
    Cxx = 1/(Ne-1)*(xMinus-repmat(mean(xMinus,2),1,Ne))*(xMinus-repmat(mean(xMinus,2),1,Ne))';
%% 
% Similarly, we will use the observation perturbations to compute the error 
% covariance matrix of the observations:

    % empirical observation error covariance matrix
    Sll = 1/(Ne-1)*dY(t,:)*dY(t,:)';
%% 
% We now have sufficient information to fill the Kalman gain. *Please add the 
% Kalman gain equation.* You will have to use the matlab function |inv()| to compute 
% the inverse of matrices.

    % Kalman gain
    %%%%%%%%% TO BE FILLED %%%%%%%%%%%%%%%%%
    
%% 
% We can now apply the EnKF update to the model state stored in xMinus, thus 
% generating xPlus.

    % Kalman update
    xPlus = xMinus + K * (Y_ens(t,:)-A*xMinus);
%% 
% We can now initialize the next run with the updated information

    % re-initialization
    S_t = xPlus(1,:);
%% 
% Before we go forward, we will save all the update information in the preset 
% matrices, so that we can plot them after the process is done.

    % empirical covariance matrices of filter update
    Cxpxp = 1/(Ne-1)*(xPlus - repmat(mean(xPlus,2),1,Ne))*(xPlus - repmat(mean(xPlus,2),1,Ne))';
    
    % save time series
    xMinus_3d(:,:,t-sim(1)+1) = xMinus;
    xPlus_3d(:,:,t-sim(1)+1) = xPlus;
    Cxx_3d(:,:,t-sim(1)+1) = Cxx;
    Cxpxp_3d(:,:,t-sim(1)+1) = Cxpxp;
    Sll_2d(t-sim(1)+1,1) = Sll;
    
end
%% Visualize
% Let us visualize the results.

F_Visualize_Ex3(Y_ens,S_obs,S_true,S_ens,xPlus_3d,Sll_2d,Cxx_3d,Cxpxp_3d)
%% Extra exercises
%% 
% # Go to the _[Marker for extra exercises]_ and triple the observation uncertainties 
% (that is, triple the variable |sigma_sObs|). Then run the whole script again, 
% and see what happens.
% # Now, make the observation uncertainty 3 times smaller. Run the whole script 
% again and see what happens. 
%% Reflection questions
%% 
% # Do you think the observation and their uncertainty is realistic, given the 
% ground truth? 
%% 
%% 
% # What happens to the model estimates in the DA process? 
%% 
%% 
% # What happens to the model uncertainties in the DA process? 
%% 
%% 
% # What do you think is the benefit of DA for models? 
%% 