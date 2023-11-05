%% Global estimation for data anaysis in thermometry
% 
% This code simulates the outcomes of mu energy measurements on a gas of n
% non-interacting spin-1/2 particles in thermal equilibrium. These are then
% processed into a temperature estimate by by using global quantum thermometry
% as put forward in:
%
%   J. Rubio, J. Anders, and L. A. Correa, Phys. Rev. Lett. 127, 190402 (2021).
%
% Running the code generates a plot with the result. 
% 
% Notes:
%
%   - data_opt selects either (1) a fresh simulation, or (2) the simulated
%   data set employed in Fig. 1b of Phys. Rev. Lett. 127, 190402 (2021).
%
%   - In order to compare with local estimation, this code also calculates
%   the temperature one would get from an optimal local estimator with initial
%   'hint' at the temperature T_true.
%
%   - The units have been chosen such that T is dimensionless.
%
% Jesús Rubio, PhD
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
%
% Created: May 2021
% Last updated: November 2023
clear all

% Data set option
data_opt=2;

if data_opt == 1
    % New data set
    n = 250; % number of particles
    mu = 5000; % number of trials
    
elseif data_opt == 2
    % Data set in manuscript
    n = 150;
    mu = 500;
end

% Outcomes space
r = 0:n;

%% Temperature estimation from simulated measurements: Bayesian approach

% Prior information
Tmin = 0.1; % lower limit
Tmax = 10; % upper limit
dT = 10^(-3);
dimT = round((Tmax - Tmin)/dT); % differential
T = linspace(Tmin, Tmax, dimT); % dimensionless temperature space
measure = sparse(1./T); % measure representing complete ignorance for scale parameters
prior = measure/trapz(T,measure); 

% Likelihood function
likelihood = zeros(n+1, length(T));
for xAux = 1:n+1
    likelihood(xAux,:) = sparse(exp(-r(xAux)./T-n*log(1+exp(-1./T))+ncklog(n,r(xAux))));
end

%% Simulation

% 'True' temperature 
if data_opt == 1
    index_real = 3900;
elseif data_opt == 2
    index_real = 3900;
end

% Energy measurements (total energy of all n particles in each trial)
if data_opt == 1
    prob_sim = likelihood(:, index_real);
    rng('shuffle') % seed for the random generator
    for runs = 1:mu
        
        auxiliar = cumsum(prob_sim) - rand;
        
        for x = 1:n+1
            if auxiliar(x) > 0
                outcome_index = x;
                break
            end
        end
        
        outcomes(runs) = r(outcome_index); %#ok<SAGROW>
        outcomes_index(runs) = outcome_index; %#ok<SAGROW>
        
    end
    
elseif data_opt == 2
    outcomes = load('data_sample.txt');
    outcomes_index = zeros(1, mu);
    for runs = 1:mu 
        indexAux = find(outcomes(runs) == r);
       outcomes_index(runs) = indexAux;
    end
end
        
% Inference
prob_temp = prior;
opt_est = zeros(1, mu); opt_err = zeros(1, mu);
locEst = zeros(1, mu); loc_err = zeros(1, mu);
for runs = 1:mu
        
    % Likelihood, joint, evidence and posterior functions
    joint = prob_temp.*likelihood(outcomes_index(runs), :); % joint probability
    evidence = trapz(T,joint); % normalisation of Bayes theorem
    
    if evidence > 1e-16
        posterior = joint/evidence; % posterior probability
    else
        posterior = 0;
    end
    
    prob_temp = posterior; % this updates the posterior with the info of each new trial (Bayes theorem)
            
    % Optimal estimator
    aux = sparse(posterior.*log(T));
    opt_log_est = trapz(T, aux);
    opt_est(runs) = exp(opt_log_est); % estimator in Eq.(6)
    
    % Optimal uncertainty
    opt_err(runs) = trapz(T, aux.*log(T)) - opt_log_est^2;
    
    % Local theory
    local_point = T(index_real) - 1;
    loc_est(runs) = local_point+4*local_point^2*mean(outcomes(1:runs))*cosh(1/(2*local_point))^2/n - local_point^2*(1+exp(-1/local_point));   
    loc_err(runs) = 2*local_point^2*cosh(1/(2*local_point))/sqrt(n*runs);
    
end
opt_err_bar = opt_est.*sqrt(opt_err);

% Plots
shaded_error_bar(1:mu, opt_est, opt_err_bar, 'lineProps', 'b');
hold on
shaded_error_bar(1:mu, loc_est, loc_err, 'lineProps', 'k');
plot(1:mu, T(index_real)*ones(1,mu), 'r-', 'LineWidth', 1.5)
hold off
fontsize=25;
xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$k_B\tilde{\theta}(\textbf{\emph{r}})/(\hbar \omega)$', 'Interpreter', 'latex', 'FontSize', fontsize);
legend('Opt. estimator, $\tilde{\vartheta}(\textbf{\emph{r}})$', 'Loc. estimator, $\tilde{\theta}_L(\textbf{\emph{r}})$', 'True temperature, $T$', 'Interpreter', 'latex', 'Location', 'southwest')
xlim([1 mu])
ylim([T(index_real)-1 T(index_real)+1])
set(gca, 'FontSize', fontsize, 'FontName', 'Times')
box on
grid