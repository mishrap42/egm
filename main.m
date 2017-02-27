% This program graphs the consumption function for unemployment model
close all;
clear all;
run EGM

%% Parameters
s.Beta = 0.7;     % hyperbolic discounting parameter
s.Delta=0.96;     % time preference -- discount rate of future actions.
s.Rho=2;          % CRRA -- higher it is, more risk averse.
s.n=20;           % number of grid points
s.G=1.03;         % permanent income growth rate; period-to-period growth
s.p=0.005;        % probability of losing all income
s.W= 1.00;        % wage function (constant)
s.R=1.04;         % interest rate function (constant)
s.numPds = 99;    % number of periods over which gridpoints are run
s.Sigma = 0.357;

% set up triple exponential grid to have more points in higher alpha range
% when spline needs to be more sensitive to slope changes
s.AlphaVec=exp(exp(exp(linspace(0.00,log(log(log(10+1)+1)+1),s.n))-1)-1)-1; 

%% Exogenous Variables

permval=[0.90,1.00,1.10];   
% permanent shock values; when permanent shocks to probability occur, can
% either setback, push forth, or stay same with probabilities on next line
permprob=[0.25,0.50,0.25];  % permanent shock probabilities

%transitory shock, if it occurs, is the permanent values times the
%probability of zero income. note that thanks to liquidity constraints
%inducing a saving structure if p is 0 this makes no difference
tranval=[permval/(1-s.p),0.0];   % transitory  shock values
tranprob=[permprob*(1-s.p),s.p];   % transitory shock probabilities

s.zval = s.W*tranval;
s.zprob = s.W*tranprob;

% integrate into a structure as a single structural argument
% s = struct('Beta',0.7,'Delta',0.96, 'Rho',2,'n',20','G',1.03,'p',0.005,...
%             'W',1.00,'R',1.04,'numPds',99,'Sigma',0.357,...
%             'AlphaVec',AlphaVec,'zval',zval,'zprob',zprob

%% Call Functions
% beta-delta / present biased
[M,C,L] = EGM.gridpoints(s);

pbias = C(:,end);
pbias_m = M(:,end);
pbias_m_expec = s.R*s.AlphaVec' - s.W .* L(:,end) + s.zval*s.zprob';

% rational
s.Beta = 1;

[M,C,L] = EGM.gridpoints(s);

rational = C(:,end);
rational_m = M(:,end);
rational_m_expec = s.R*s.AlphaVec' - s.W .* L(:,end) + s.zval*s.zprob';

% Consumption Functions
figure()
plot(pbias_m,pbias);
title('Consumption Function');
xlabel('beginning-of-period resources (m)');
ylabel('c');
xlim([0 10]);
hold on;
plot(rational_m,rational);
hold off;
legend('Present Biased, \beta = 0.7','Rational');

% Consumption as Proportion of Starting Capital
figure();
plot(pbias_m,pbias./pbias_m);
title('Consumption Function');
xlabel('beginning-of-period resources (m)');
ylabel('share of resources consumed (c/m)');
xlim([0 10]);
hold on;
plot(rational_m,rational./rational_m);
hold off;
legend('Present Biased, \beta = 0.7','Rational');

%Market goods at end of period vs. beginning
figure()
plot( pbias_m, pbias_m_expec );
title('Market Goods Dynamics');
xlabel('market resources');
ylabel('expected market resources');
xlim([0 10]);
hold on;
plot( rational_m, rational_m_expec );
plot(0:length(rational), 0:length(rational),':k');
hold off;
legend('Present Biased, \beta = 0.7','Rational', '45 Degree Line', 'location','northwest');

% add steady state to plot
rational_m_ss = fsolve(@(x) interp1(rational_m,rational_m_expec,x) - x, 2);
pbias_m_ss = fsolve(@(x) interp1(pbias_m,pbias_m_expec,x) - x, 2);
hold on; 
plot([rational_m_ss rational_m_ss],[0 rational_m_ss],'--k');
plot([pbias_m_ss pbias_m_ss],[0 pbias_m_ss],'--k');
hold off;
