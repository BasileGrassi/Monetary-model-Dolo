%The parameters of the IRF
%paramirf.range= the range of time that will be plot
%paramirf.style= the style of the line (markers, dashed, etc...) e.g: '-bs'
%paramirf.width= the width of the line
%paramirf.paramirf.istatediscreate=the index of a discreate state variable

%%%Scenario configuration

close all
clear all

%addpath('..')
addpath('../solver_lib')



%% Configure the exogenous fundalmentals
exo.n=2;
N=50;
exo.T=N;

%Define exogeneous path

model_name='baseline';

addpath(['../models/' model_name]);

load( [ model_name '_sol']);

ieta=strmatch('eta', model.states,'exact');
irho=strmatch('rho',model.parameters,'exact');
irho=15;

eta10y=linspace(1,0.925,10);
eta10y= eta10y';
eta_last=model.s_ss(ieta);

epsi_constr=bargain_constraint(eta10y,eta_last,irho,ieta, [1,10],model);
epsi=[epsi_constr;zeros(N-10,1)];

crisis=[zeros(30,1);1;zeros(N-30-1,1)];
exo.e=[epsi,crisis];


%% Benchmark

%% Load the model and the decision rule
%this was compute using main
%contains model, grid and rule

model_name='baseline';

addpath(['../models/' model_name]);

load( [ model_name '_sol']);

paramirf.istatediscreate=strmatch('crisis',model.states,'exact'); %give the index of a discreate state variable

%The shape of irf's benchmarks
paramirf.range=[0 N];
paramirf.style='-k';
paramirf.width=2;

%Compute and Plot IRF
RES_benchmark=irf(grid, rule, exo, paramirf, model);

%% No uncertainty and endogenous crisis proba
model_name='baseline_no_uncertainty';

addpath(['../models/' model_name]);

load( [ model_name '_sol']);

paramirf.istatediscreate=strmatch('crisis',model.states,'exact'); %give the index of a discreate state variable



%The shape of irf
paramirf.range=[0 N];
paramirf.style='--r';
paramirf.width=2;

%Compute and Plot IRF
RES=irf(grid, rule, exo, paramirf, model);

%% No uncertainty and exogenous crisis proba
model_name='baseline_no_uncertainty_probaexo';

addpath(['../models/' model_name]);

load( [ model_name '_sol']);

paramirf.istatediscreate=strmatch('crisis',model.states,'exact'); %give the index of a discreate state variable


%The shape of irf
paramirf.range=[0 N];
paramirf.style=':b';
paramirf.width=2;

%Compute and Plot IRF
RES=irf(grid, rule, exo, paramirf, model);

