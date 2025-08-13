%% Simple “bucket” hydrological model
% We will start working with a simple “bucket” hydrological model. See lecture 
% "Introduction to Data Assimilation" for more information on the one-bucket model.
%% Setting up model variables
% First, we will set up the model variables.

clc, clear all, close all
addpath('data'); addpath(genpath('functions'));
%% 
% This simple bucket model only has one model parameter, K, which relates the 
% relation between the water stored and runoff.

% Define parameter K
K  = 0.2;  % model parameter, element of (0,1), for practical implementation [0.01 0.99]
%% 
% Let us define the timestep, load the forcing (precipitation minus evapotranspiration, 
% P-ET) data and initialize our water storage state. 

% Define timestep
Dt =   1;  % time step

% precipitation data
fnP = 'Precip.dat'; % filename of precipitation input
P=load(fnP); clear fnP

% initialize storage vector
S0 =   4;  % initial state
S = zeros(length(P),1); S(1)=S0;
%% Model run
% We are now ready to run the linear reservoir model. For that, we will call 
% the linearModel function.

% run linear reservoir model
[Q, S] = linearModel(S, P, K, Dt);
%% Visualization
% We can now visualize the outputs of the model run.

F_Visualize_Ex1(P,S,S0,Q,K);
%% Reflection questions
% *Reflection questions:*
%% 
% # How does the net precipitation P-ET relate to the storage S? 
%% 
% 
% 2. How does the runoff R relate to the storage S? 
% 
% 
% Note: real hydrological models have a more complex representation of the hydrology.