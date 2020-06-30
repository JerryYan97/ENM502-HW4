% Kaushik Shankar
% Solves the system of equations defined in
% Eq.1, G.Q. Sun, Nonlinear Dynamics, 69, 1097-1104 (2012)

clear;
clc;

% Set number of PDEs that you want to solve for
% Ref: https://www.mathworks.com/help/pde/ug/createpde.html
numberofpde = 2;
model = createpde(numberofpde);

% Generate and plot the domain
% Set the width and height of the domain
% In the paper, this is 200 by 200 system
% Here, we solve the same for a 50 by 50 system
width = 50; 
height = 50;
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.geometryfromedges.html
gdm = [3 4 0 width width 0 0 0 height height]';
g = decsg(gdm, 'S1', ('S1')');
geometryFromEdges(model,g);


% Set the coefficients of the system of PDEs
% The general form of the equation is:
% m*(?/?t)(?u/?t) + d*?u/?t ? ?·(c*?u)+ a*u = f
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.specifycoefficients.html

% Model parameter from the paper
% Ref: https://www.mathworks.com/help/pde/ug/f-coefficient-for-specifycoefficients.html
A = 1;
f1 = A;
f2 = 0;
f = [f1; f2];

% Ref: https://www.mathworks.com/help/pde/ug/m-d-or-a-coefficient-for-systems.html
d1 = 1;
d2 = 1;
d = [d1; d2];

% Ref: https://www.mathworks.com/help/pde/ug/m-d-or-a-coefficient-for-systems.html
% Position-dependent non-linear coefficients can be assigned using
% function-handles if needed
% The format for specifying the function is given below
a = @acoeffunction;

% Model parameter from the paper 
D1=6; 
D2=1;
beta = 32;
mu = 1.8;

% Ref: https://www.mathworks.com/help/pde/ug/c-coefficient-for-systems-for-specifycoefficients.html
c = [D1; 0; 0; D1; D2; 0; 0; D2];


specifyCoefficients(model,'m',0, ...
    'd',d, ...
    'c',c, ...
    'a',a, ...
    'f',f);

% Set position-dependent initial conditions
% For this model, the initial condition is a random perturbation about the
% stable steady state of the system
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.setinitialconditions.html
setInitialConditions(model, @randomperturb);

% Set boundary conditions
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.applyboundarycondition.html
% Dirichlet (or value) BC example
% applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',[0; 0]);

% Neumann (or flux) BC example
applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,'g',[0; 0],'q',[0;0;0;0]);

% Generate the finite element meshing for the domain
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.generatemesh.html
generateMesh(model,'Hmax',2)

% Ref: https://www.mathworks.com/help/pde/ug/pdemesh.html
figure;
pdemesh(model);

% List of times you need the solution to be evaluated at
tlist = 0:5:500;

% Solve the system of PDEs using conditions specified earlier
% The solution is stored only in the list of time points specified in tlist
% Ref: https://www.mathworks.com/help/pde/ug/pde.pdemodel.solvepde.html
results = solvepde(model, tlist);

% Plot the solution to the PDE at times that are of interest to you
% Ref: https://www.mathworks.com/help/pde/ug/pdeplot.html
figure;
pdeplot(model, 'XYData', results.NodalSolution(:,1,end));
colormap('jet');

save('case1Model.mat', 'model');
save('case1_50x50Result.mat', 'results');

function ic = randomperturb(location)
nr = numel(location.x);
ic = zeros(2,nr);
beta = 32;
mu = 1.8;
A = 1;
d = 1;
%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
randomVec = 0.001 * rand(2, nr);
S_SS = (A * beta - sqrt(A.^2 * beta.^2 - 4 * d.^3 * beta - 8 * d.^2 * beta * mu - 4 * d * beta * mu.^2)) / (2 * d * beta);
I_SS = (2 * d * (d + mu)) / (A * beta - sqrt(A.^2 * beta.^2 - 4 * d.^3 * beta - 8 * d.^2 * beta * mu - 4 * d * beta * mu.^2));
for idx = 1:nr
    ic(1,idx) = randomVec(1,idx) + S_SS;
    ic(2,idx) = randomVec(2,idx) + I_SS;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function amatrix = acoeffunction(location, state)
nr = numel(location.x);
amatrix = zeros(2,nr);
u = state.u;
beta = 32;
mu = 1.8;
d = 1;
%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
for idx = 1:nr
    amatrix(1,idx) = d + beta * u(2,idx).^2;
    amatrix(2,idx) = d + mu - beta * u(1,idx) * u(2,idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end