function res = example_stl_vanderpol
% example of signal temporal logic
% checking of the vanderpol oscillator
%
% Syntax:
%    res = example_stl_vanderpol()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean

% Authors:       Besset Antoine
% Written:       22-Jul-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------
f = @(x,u)[
    x(5)*(x(3)-x(1));
    x(5)*(x(4)-x(2));
    x(4);
    0.5*(1 - x(3)^2)*x(4) - x(3);
    0   % K is constant
];
sys = nonlinearSys(f);
% Parameter ---------------------------------------------------------------

params.tFinal = 5.01; %%% Minimum simulation time to have a conclusive evaluation at 2.9s
%params.startLoc = 1;
params.R0 = zonotope([4;4;4;4;1.15], diag([0.01; 0.01; 0.01; 0.01; 0.15]));

% Reachability Settings ---------------------------------------------------

options.alg = 'poly-adaptive';

% settings for continuous reachability
options.taylorTerms = 4;
options.zonotopeOrder = 10;
options.timeStep = 0.01;

% settings for hybrid systems
options.enclose = {'box','pca'};
options.guardIntersect = 'zonoGirard';

% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

A1 = [
   1  0  0  0  0;   % x1 <= 0.75
  -1  0  0  0  0;   % x1 >= -0.75
   0  1  0  0  0;   % x2 <= 1.1
   0 -1  0  0  0    % x2 >= -1
];
b1 = [0.75; 0.75; 1.1; 1];
p1 = stl('p1', atomicProposition(polytope(A1, b1)));

A2 = [
   1  0  0  0  0;
  -1  0  0  0  0;
   0  1  0  0  0;
   0 -1  0  0  0
];
b2 = [1.9; -1.5; -0.8; 1.5];
p2 = stl('p2', atomicProposition(polytope(A2, b2)));

A3 = [
   1  0  0  0  0;
  -1  0  0  0  0;
   0  1  0  0  0;
   0 -1  0  0  0
];
b3 = [4; -3; 0.5; 1.5];
p3 = stl('p3', atomicProposition(polytope(A3, b3)));

phi = {and(globally(~p1, interval(3,5)), globally(or( finally(p2, interval(1,2)), ~p3 ), interval(3, 4)))};

res = true;

% 1. Incremental with four-valued signals fastest verification: 5.775s or 4.6625s

tic

for i = 1:length(phi)
    assertLoop(modelChecking(R,phi{i},'incremental','propFreq',100,'verbose',true),i);
end

tComp = toc;
disp(['computation time with incremental: ',num2str(tComp)]);

% ------------------------------ END OF CODE ------------------------------
