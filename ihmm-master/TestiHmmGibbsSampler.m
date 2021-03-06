% Test the iHMM Gibbs sampler on the expansive HMM experiment from the
% original iHMM paper [Beal et al. 2002]. The plotted transition and
% emission matrices should be equivalent up to permutation.(排列？？)

T = 500;                       % Length of HMM
K = 4;                         % # of states
L = 8;                         % # of symbols in emission alphabet

stream = RandStream('mt19937ar','seed',21);
RandStream.setGlobalStream(stream);

% Parameters for the HMM that generates the data.
A = [ 0.0 0.5 0.5 0.0;
      0.0 0.0 0.5 0.5;
      0.5 0.0 0.0 0.5;
      0.5 0.5 0.0 0.0 ];
E = [ 1 0 0 0 0 0 1 1;
      1 1 1 0 0 0 0 0;
      0 0 1 1 1 0 0 0;
      0 0 0 0 1 1 1 0 ] / 3;
pi = [1.0; zeros(K-1,1)];

% Generate data.
[Y, STrue] = HmmGenerateData(1,T,pi,A,E);%###这里生成一个样本，返回的数据是观测序列与隐状态序列



STrue
%fprintf('%d  ',Y(1:10))；



% Sample states using the iHmm Gibbs sampler.
hypers.alpha0_a = 4;
hypers.alpha0_b = 2;
hypers.gamma_a = 3;
hypers.gamma_b = 6;
hypers.H = ones(1,L) * 0.3;
tic
[S ,stats] = iHmmSampleGibbs(Y, hypers, 500, 1, 1, ceil(rand(1,T) * 10));
toc


%fprintf(S{1}.S);
%fprintf(S{1}.K);

S{1}.S

% Plot some stats
figure(1)
subplot(3,2,1)
plot(stats.K)
title('K')
subplot(3,2,2)
plot(stats.jll)
title('Joint Log Likelihood')
subplot(3,2,3)
plot(stats.alpha0)
title('alpha0')
subplot(3,2,4)
plot(stats.gamma)
title('gamma')
subplot(3,2,5)
imagesc(SampleTransitionMatrix(S{1}.S, zeros(1,S{1}.K))); colormap('Gray');
title('Transition Matrix')
subplot(3,2,6)
imagesc(SampleEmissionMatrix(S{1}.S, Y, S{1}.K, hypers.H)); colormap('Gray');
title('Emission Matrix')