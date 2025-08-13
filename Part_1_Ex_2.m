%% Ensemble run
% Here, we will investigate what the ensemble represents in a Monte-Carlo approach. 
% See lecture "Introduction to Data Assimilation" for more information on the 
% ensemble approach.
%% Setting up the ensemble run
% As a basis for the ensemble run, we will use a similar setting to the single 
% run in the previous exercise (i.e. definition of model parameter, timestep, 
% net precipitation data and initial storage value). *Please copy-paste the missing 
% lines from Exercise 1.*

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
%% 
% To generate a model ensemble, first, we have to define the ensemble size.

Ne = 30; % ensemble size
%% 
% We also need to load some Gaussian noise files, which we will use to perturb 
% the input data, parameters and initial water storage state. *[Marker for extra 
% exercises].*_

% Load gaussian noise files
load('data/noise_K');
load('data/noise_S0_ens');
load('data/noise_mP');
%% 
% Here, the properties for the Gaussian noise are defined, and then use to perturb 
% the different variables:

% Set noise amplitude
Kmin = 0.01; Kmax = 0.99; % min and max value of model parameter
S0min = 2; S0max = 8;       % min and max value of initial state
Pmin = 0.8; Pmax = 1.2;   % multiplicator

% Perturb variables
K_ens = Kmin + (Kmax-Kmin).*noise_K;     % initial parameter values
S0_ens = S0min + (S0max-S0min).*noise_S0_ens; % initial state values
mP =  Pmin + (Pmax-Pmin).*noise_mP;       % precipitation multiplicator
P_ens = repmat(P,1,size(mP,2)) .* repmat(mP,length(P),1); % input ensemble
clear mP P
%% Initialization of ensemble spread
% Before starting with the ensemble run, we need to launch a "warm-up" run to 
% make sure that the model storage state has an appropriate ensemble spread that 
% reflects the uncertainties of (1) initial storage; (2) precipitation and (3) 
% model parameters. Let us run a warm-up over a period of 12 timesteps:

ini = [1; 12]; % initialization phase of model

% initialize storage vector for initialization phase
S_ini = zeros(ini(2),Ne);
S_ini(1,:) = S0_ens;

% run linear reservoir model for each ensemble member
for ii = 1:Ne
    [~, S_ini(ini(1):ini(2),ii)] = linearModel(S_ini(ini(1):ini(2),ii), P_ens(ini(1):ini(2),ii), K_ens(1,ii), Dt);
end
%% Ensemble run
% We will now initialize the final ensemble run with the ouput of the "warm-up" 
% run. The final ensemble run will be performed over a period of 24 timesteps.

sim = [1; 24]; % simulation phase of model

% initialize storage vector for simulation phase
S_ens = zeros(sim(2),Ne);
S_ens(1,:) = S_ini(end,:);
%% 
% We are now ready to launch the ensemble run. We will iterate and run the linearModel 
% function for each ensemble member.

% run linear reservoir model for each ensemble member
for ii = 1:Ne
    [Q_ens(sim(1):sim(2),ii), S_ens(sim(1):sim(2),ii)] = linearModel(S_ens(sim(1):sim(2),ii), P_ens(sim(1):sim(2),ii), K_ens(1,ii), Dt);
end
%% Visualization
% Let us visualize the ensemble run results.

F_Visualize_Ex2(P_ens,S_ens,Q_ens,sim)
%% Extra exercises
%% 
% * Go back to the _"Marker for extra exercises"._ Triple the precipitation 
% uncertainty (variable |noise_mP|) and run the whole script to see what happens.
% * At the same place, reduce the ensemble size from 30 to 5. You can do this 
% by modifying the size of the loaded noise variables (|noise_K|, |noise_mP| and  
% |noise_S0_ens|). Run the whole script to see what happens.
%% Reflection questions
%% 
% # What does the ensemble spread represent? 
%% 
%% 
% # How does the uncertainty of S and R depend on the uncertainty of P – ET?      
%% 
%% 
% # What could be the benefits and challenges of increasing the ensemble? 
%% 
% 