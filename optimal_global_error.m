function [e_opt] = optimal_global_error(n, Tmin, Tmax)
% Given
%
%   - n non-interacting spin-1/2 particles in thermal equilibrium, and
%   - a prior temperature range [Tmin, Tmax],
% 
% the function
%
%     OptimalGlobalError(n,Tmin,Tmax)
%
% returns the minimium mean logarithmic uncertainty, as indicated by Eq.(7) 
% of Phys. Rev. Lett. 127, 190402 (2021). Note that the chosen prior 
% represents complete ignorance within the interval [Tmin, Tmax]. 
%
% Remarks:
%
%   - The units have been chosen such that T is dimensionless.
%
%   - The discretisation of the continuous variable T has been chosen such 
%     that the calculation is accurate for the interval [0.1, 10].
%
%   - This code relies on the auxiliary function: ncklog(n, k).
%
% Dr Jesús Rubio
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
%
% Created: September 2020
% Last updated: November 2023

%% Initialisation
if n > 10^5
    warning('This code may not be accurate for systems with more than 10^5 particles.')
end

% Number of spin-1/2 particles
r = 0:n;

% Parameter space
if Tmax>10
    warning('This code may be very slow when Tmax>>10.')
end
dT = 10^(-3)/2;
dimT = round((Tmax - Tmin)/dT);
T = linspace(Tmin, Tmax, dimT);

%% Inference

% Prior
measure = sparse(1./T); % complete ignorance for a scale parameter   
prior = measure/trapz(T, measure); 

log_info = 0;
for index = 1:length(r)
       
    % Likelihood, joint, evidence and posterior functions
    fermion_likelihood = sparse(exp(-r(index)./T-n*log(1 + exp(-1./T)) + ncklog(n,r(index)))); 
    joint = prior.*fermion_likelihood;    
    evidence = trapz(T, joint);

    if evidence > 1e-16
        posterior = joint/evidence;
    else
        posterior = 0;
    end
        
    % Optimal logarithmic estimator
    aux1 = sparse(posterior.*log(T));
    opt_log_est = trapz(T, aux1);
    
    % Auxiliary term to calculate the mean logarithmic error 
    log_info = log_info + evidence*opt_log_est^2;
end

% Optimal global mean logarithmic error
aux2 = prior.*log(T).*log(T);
e_opt = trapz(T,aux2) - log_info;

end