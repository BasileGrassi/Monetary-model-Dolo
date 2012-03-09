%The parameters of the IRF
%paramirf.range= the range of time that will be plot
%paramirf.style= the style of the line (markers, dashed, etc...) e.g: '-bs'
%paramirf.width= the width of the line
%paramirf.paramirf.istatediscreate=the index of a discreate state variable

%%%Scenario configuration

close all

addpath('..\')
addpath('..\..\..\solver_lib')

%% Load the model and the decision rule
%this was compute using iterative_main
%contains model, grid and rule

model_name='baseline';

load( [ model_name '_sol']);

%% Configure the exogenous fundalmentals
exo.n=2;
N=50;
exo.T=N;
paramirf.istatediscreate=4; %give the index of a discreate state variable

%% Benchmark
%Define exogeneous path
ieta=3;
irho=15;


eta10y=linspace(1,0.925,10);
eta10y= eta10y';
eta_last=model.s_ss(ieta);

epsi_constr=bargain_constraint(eta10y,eta_last,irho,ieta, [1,10],model);
epsi=[epsi_constr;zeros(N-10,1)];

crisis=[zeros(30,1);0;zeros(N-30-1,1)];

exo.e=[epsi,crisis];

%The shape of irf's benchmarks
paramirf.range=[0 N];
paramirf.style='-';
paramirf.width=2;

%Compute and Plot IRF
RES_benchmark=irf(grid, rule, exo, paramirf, model);

%% An other exogeneous path
%Define exogeneous path
ieta=3;
irho=15;

eta10y=linspace(1,0.925,10);
eta10y= eta10y';
eta_last=model.s_ss(ieta);

epsi_constr=bargain_constraint(eta10y,eta_last,irho,ieta, [1,10],model);
epsi=[epsi_constr;zeros(N-10,1)];

crisis=[zeros(30,1);1;zeros(N-30-1,1)];
exo.e=[epsi,crisis];

%The shape of irf
paramirf.range=[0 N];
paramirf.style='--r';
paramirf.width=2;

%Compute and Plot IRF
RES=irf(grid, rule, exo, paramirf, model);


