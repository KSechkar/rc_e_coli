%% figS2cd.m

% MCMC FITTING
% Figure S2: c,d

% From the MCMC fit, find the eigenvalues of the Fisher Information Matrix
% (displayed in Figure c, not plotted using Matlab) and the parameter
% sensitivities (Figure d)

% Warning: requires the DREAM package

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% LOAD the MCMC chains
load DREAM.mat
ParSet = genparset(chain); DREAMPar.N = size(chain,3);
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : DREAMPar.d ); % take the last 25% of the posterior samples
N_Pars = size(Pars,1); % get number of posterior samples

%% Figure S2 c - FIM AND EIGENVALUES
% estimate Fisher Info Matrix as pseudo-inverse of the variance-covariance matrix for the posterior samples
FIM=pinv(cov(ParSet(:,1:4))) % print FIM
eigenvals=eig(FIM); % find FIM eigenvalues
eigenval_magnitudes=abs(eigenvals); % FIM eigenvalue magnitudes
eigenval_magnitudes_norm=eigenval_magnitudes./max(eigenval_magnitudes); % normalised eigenvalue magnitudes
eigenvals_log=log10(eigenval_magnitudes_norm) % print normalised eigenvalue logarithms

%% Figure S2 d - PARAMETER SENSITIVITIES
% get the matrix C that diagonalises FIM
[C,DIAG]=eig(FIM);

% get parameter sensitivities
sens2=zeros([1 4]);
for j=1:size(sens2,2)
    for i=1:size(sens2,2)
        sens2(j)=sens2(j)+eigenvals(i)*(C(i,j)^2);
    end
end
disp('Parameter sensitivities (for a_a, a_r, nu_max, K_t and kcm respectively) ')
sens_log=log10(sens2./sum(sens2)) % print parameter sensitivities

% make a bar chart
param_names={'K_\epsilon=K_\nu'; 'k_{cm}'; '\alpha_r : \alpha_a'; '\nu_{max}';};


Fbar=figure('Position',[0 0 385 255]);
set(Fbar, 'defaultAxesFontSize', 9)
set(Fbar, 'defaultLineLineWidth', 1)
bar(sens_log([3,4,1,2]))
set(gca,'xticklabel',param_names)
ylabel('log(normalised sensitivity)')
axis square
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on
