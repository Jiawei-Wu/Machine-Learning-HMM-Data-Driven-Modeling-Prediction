function dir = newdirichlet_sample(alpha)
% Samples a dirichlet distributed random vector.

dir = gamrnd(alpha, 1);%Generate Gamma Distribution
dir = dir ./ sum(dir);