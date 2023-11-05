function [C] = ncklog(n,k)
% Given the integers n and k, with k in [n, 0), 
%
%   ncklog(n, k) 
%
% returns the natural logarithm of the binomial coefficient C(n, k), 
% i.e., log[C(n, k)].
%
% Dr Jesús Rubio
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
%
% Created: June 2020
% Last updated: November 2023

aux = 0;
for x = k + 1:n
    aux = aux + log(x) - log(x-k);
end
C = aux;

end