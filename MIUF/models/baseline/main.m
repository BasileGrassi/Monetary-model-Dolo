%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need a model written as:
%
% snext=G(s,x,e)
% E F(s,x,e,snext,xnext)=0
%
% Do not forget to add the compecon libraries
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Addpath libraries

addpath('../../solver_lib')

%% set model name
model_name = 'baseline';

%% Parameters
none=[];
model = feval( [model_name '_model']);

%% Define regular shocks (non discrete)
N_shocks = 5;
[epsi,w]=hernodes(N_shocks);
isigma=strmatch('sigma',model.parameters, 'exact');
%isigma=23;
sigma = [[model.params(isigma)]];
e = sigma*epsi;


%% Defined the grid
nstate=length(model.states);
k_pts       = 10;
l_pts       = 10;
z_pts       = 4;
d_pts       = 2;


%Size of the State Space
faktor = 1;

% K is beginning of period capital stock (level)
ik=strmatch('k', model.states,'exact');
kmin    = model.s_ss(ik) - 0.5*faktor;
kmax    = model.s_ss(ik) + 6.0*faktor;
kstep   = (kmax-kmin)/k_pts;
k_grid  = kmin:kstep:kmax;
nk      = length(k_grid);

% L is beginning period of debt (level)
il=strmatch('l', model.states,'exact');
lmin    = model.s_ss(il) - 0.05*faktor;
lmax    = model.s_ss(il) + 0.60*faktor;
lstep   = (lmax-lmin)/l_pts;
l_grid  = lmin:lstep:lmax;
nl      = length(l_grid);

% Z is current bargaining power
zmin    = 1 - 0.08*faktor;
zmax    = 1 + 0.0065*faktor;
zstep   = (zmax-zmin)/z_pts;
z_grid  = zmin:zstep:zmax;
nz      = length(z_grid);

% D is current default status
d_grid  = [0,1];
nd=size(d_grid,2);
   
%% Create the grid
grid=gridmake(k_grid',l_grid', z_grid', d_grid');
ns = size(grid,1);
    
% Define interpolator
order=[nk nl nz nd];
gridmin=[kmin lmin zmin 0];
gridmax=[kmax lmax zmax 1];
cdef=fundefn('lin',order,gridmin,gridmax);


% Convergence criteria
tol=1e-10;
maxiteration=1000;

% Initialization using first order d.r.
x_ss = model.x_ss;
s_ss = model.s_ss;
X_s = initial_guess(model, model.s_ss, model.x_ss);    %X_s = model.X{2};
xinit=x_ss*ones(1,ns)+X_s*(grid'-s_ss*ones(1,ns));
x=xinit';


iteration=1;
converge=0;

hom_n = 4;
homvec = linspace(0,1,hom_n);
hom_i = 1;
hom = 0;
err0 = 1e6;

disp('Starting policy rule iteration.');
disp('_________________________________________________________');
disp('iter       error        gain    hom   inner    elap.(s)  ');
disp('_________________________________________________________');

tic;
t0 = tic;

while converge==0 && iteration < maxiteration
    
    [coeff,B]=funfitxy(cdef, grid, x);
    
    %fobj = @(xt) step_residuals_endo(grid, xt, e, w, model.params, model, coeff, cdef, hom);
    %[x_up, nit] = newton_solver(fobj, x, 50);

    fobj = @(xt) step_residuals_endo_nodiff(grid, xt, e, w, model.params, model, coeff, cdef, hom);
    [x_up, nit] = newton_solver_diff(fobj, x, 50);
    
    err=sum(sum(abs(x-x_up)));
    if (err < tol);
        converge=1;
    end;
    
    t1 = tic;
    elapsed = double(t1 - t0)/1e6;
    t0 = t1;
    
    gain=err/err0;
    fprintf('%d\t%e\t%.2f\t%.2f\t%d\t%.2f\n', iteration, err, gain, hom, nit, elapsed)

    
    if (err < 1) && (hom_i < hom_n);
        hom_i = hom_i + 1;
        hom = homvec(hom_i);
    end;
    
    err0 = err;  
    
    x=x_up;
    iteration = iteration+1;
    
    
end;
disp('___________ ________________________________');
toc;

if iteration > maxiteration
    disp('The model could not be solved');
end

%% Save the grid
grille.nstate= nstate;
grille.pts= [k_pts l_pts z_pts d_pts];
grille.max= gridmax;
grille.min=gridmin;
grille.step= [kstep lstep zstep 1];
grille.order=[nk nl nz nd];
grille.grid=grid;
grille.npts=ns;

grid=grille;


%% Save the rule
rule.x=x;
rule.cdef=cdef;
rule.coeff=coeff;


save([model_name '_sol'],'model','grid','rule');
