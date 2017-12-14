function [hidden_seq,tr_mat,emit_mat] = hmmtrain_Gibbs(Y,symbol,ifplot)


%The number of hidden states
num_init_symbol=symbol;
%Number of iterations
num_iterations=300;

num_sample_out=1;
hypers.alpha0_a = 4;
hypers.alpha0_b = 2;
hypers.gamma_a = 3;
hypers.gamma_b = 6;
hypers.H = ones(1,num_init_symbol) * 0.3;
length_hmm=size(Y,2);
init_stat_seq=ceil(rand(1,length_hmm)*20);

[S ,stats] = iHmmSampleGibbs(Y, hypers, num_iterations, num_sample_out, 1, init_stat_seq);

hidden_seq=S{1}.S;
[tr_mat, emit_mat] = hmmestimate(Y, S{1}.S);
%#############
%The Returned Hidden State Sequence

% S{1}.S
% fprintf('The Transmission Matrix:')
% tr_mat
% fprintf('The Emission Matrix:')
% emit_mat

% Plot JLL
if ifplot
    figure(5)
    plot(stats.jll)
    title('Joint Log Likelihood')
end