function [RMSE,MAPE,MABE,r] = xxy(states,symbols,ifplot,ifpause,every_day_analysis,temper_file,radia_file)

%%%%params%%%%
if nargin < 1
    states = 20;
    symbols = 25;
    
    ifplot=true;
    ifpause=true;
    
    temper_file = 'temp_LA_SUS.csv';
    radia_file = 'radiation.csv';
    every_day_analysis = true;
elseif nargin<3
    
    ifplot=false;
    ifpause=false;
    
    temper_file = 'temp_LA_SUS.csv';
    radia_file = 'radiation.csv';
    every_day_analysis = true;
elseif nargin<4
    
    ifpause=false;
    every_day_analysis = true;
    
    temper_file = 'temp_LA_SUS.csv';
    radia_file = 'radiation.csv';
elseif nargin<6
    
    temper_file = 'temp_LA_SUS.csv';
    radia_file = 'radiation.csv';
end
    
    
if every_day_analysis == true
    fprintf 'DAY ANALYSIS!\n'
else
    fprintf 'HOUR ANALYSIS!\n'
end

    
    
temper_24 = csvread(temper_file,1,1);
radia_96 = csvread(radia_file,1,1);

a = size(temper_24);
days = a(1);

% radia_24
radia_24=[];
for i=1:24
    a = radia_96(:,4*i-4+1:4*i);
    radia_24 = [radia_24 mean(a, 2)];
end


if every_day_analysis == true
    csvwrite('flag',1)
    radia_24 = mean(radia_24,2);
    temper_24 = mean(temper_24,2);
else
    csvwrite('flag',0)
end
    
% plot(

% size(radia_24)

max_temper = max(max(temper_24));
min_temper = min(min(temper_24));

max_radia = max(max(radia_24));
min_radia = min(min(radia_24));

temper_region_len = (max_temper - min_temper) / symbols;
radia_region_len = (max_radia - min_radia) / states;

temp = [];
for i = 1:symbols
    temp = [temp min_temper + i * temper_region_len];
end
temper_regions = temp;

temp = [];
for i = 1:states
    temp = [temp min_radia + i * radia_region_len];
end
radia_regions = temp;


temper_vec = reshape(temper_24,1,[]);
radia_vec = reshape(radia_24,1,[]);

temper_int_vec = [];
for i = 1:length(temper_vec)
    for j = 1:symbols
        if temper_vec(i) < temper_regions(j)
            break
        end
    end
    temper_int_vec = [temper_int_vec j];
end
% 
% test_int_vec = [];
% for i = 1:length(test_raw)
%     for j = 1:symbols
%         if test_raw(i) < temper_regions(j)
%             break
%         end
%     end
%     test_int_vec = [test_int_vec j];
% end
% % test_int_vec = reshape(test_int_vec,length(test_int_vec)/2,2)


radia_int_vec = [];
for i = 1:length(radia_vec)
    for j = 1:states
        if radia_vec(i) < radia_regions(j)
            break
        end
    end
    radia_int_vec = [radia_int_vec j];
    
end

% length(radia_int_vec)
% length(temper_int_vec)

temper_float_vec = temper_int_vec.*temper_region_len-temper_region_len/2+min_temper;
radia_float_vec = radia_int_vec.*radia_region_len-radia_region_len/2+min_radia;

fprintf('Temperature range: %f ¡ãC - %f ¡ãC. \n',min_temper, max_temper);
fprintf('Radiation range: %f W/m2 - %f W/m2. \n ',min_radia, max_radia);

index = 1:length(temper_float_vec);

if ifplot
    figure(1)
    plot(index,temper_float_vec)
    title('Temperature Changing Curve')
    if every_day_analysis
        xlabel('Day index')
    else
        xlabel('Hour index')
    end
    ylabel('Temperature (¡ãC)')

    figure(2)
    plot(index,radia_float_vec)
    title('Radiation Changing Curve')
    if every_day_analysis
        xlabel('Day index')
    else
        xlabel('Hour index')
    end
    ylabel('Radiation (W/m2)')
    
end

if ifpause
    fprintf('Presss any key to continue...\n')
    pause;
end


%%%%%%%%%%%%%Likelyhood Function estimation%%%%%%%%%
fprintf('Parameting emtimating now...\n')

[tr_mat, emit_mat] = hmmestimate(temper_int_vec, radia_int_vec);

if ifplot
    fprintf('HMM Parameters Estimated using Likelyhood Function:')
    tr_mat
    emit_mat
end

% decode = hmmdecode(temper_int_vec, tr_mat, emit_mat);
% [a,decode_index]=max(decode,[],1);
decode_index = hmmviterbi(temper_int_vec, tr_mat, emit_mat);

decode_float_vec = decode_index.*radia_region_len-radia_region_len/2+min_radia;

fprintf 'The Hidden States Sequence Estimated by Likelyhood Function£º'
decode_float_vec
csvwrite('like_predictedSequence.csv',decode_float_vec)


%%%%%%%%%EM Estimation%%%%%%%%%%

[em_tr_mat,em_emit_mat] = hmmtrain(temper_int_vec,tr_mat,emit_mat);

decode_index = hmmviterbi(temper_int_vec, em_tr_mat, em_emit_mat);

% decode = hmmdecode(temper_int_vec, em_tr_mat, em_emit_mat)
% [a,decode_index]=max(decode,[],1);
% 

if ifplot
    figure(6)
    plot(index,decode_index)
    
    title('Viterbi Original Output')
    if every_day_analysis
        xlabel('Day index')
    else
        xlabel('Hour index')
    end
    ylabel('Prob')
end

    
em_decode_float_vec = decode_index.*radia_region_len-radia_region_len/2+min_radia;

fprintf "Parameter Estimation Finished !!!¡£\n"


if ifplot
    fprintf 'HMM Transmission and Emission Matrix Estimated by EM£º\n'
    em_tr_mat
    em_emit_mat
end

fprintf 'HMM Hidden State Sequence Estimated by EM£º'
em_decode_float_vec

csvwrite('em_predictedSequence.csv',em_decode_float_vec)


%%%%%%%%%Gibbs Estimation%%%%%%%%%%
fprintf 'Gibbs Sampling is on the way...\n'
[hidden_seq,gibbs_tr_mat,gibbs_emit_mat] = hmmtrain_Gibbs(temper_int_vec,symbols, ifplot);
fprintf 'Gibbs Sampling Finished!!!¡£\n'

if ifplot
    fprintf 'HMM Transmission and Emission Matrix Estimated by Gibbs Sampling£º\n'
    gibbs_tr_mat
    gibbs_emit_mat
    
end

final_states = length(gibbs_tr_mat);

fprintf('The Number of Hidden States after Convergence: %d ¡£ \n',final_states);

gibbs_region_len = (max_radia - min_radia)/final_states;
gibbs_float_vec = hidden_seq.*gibbs_region_len-gibbs_region_len/2+min_radia;

fprintf 'HMM Hidden State Sequence Estimated by Gibbs Sampling£º'
gibbs_float_vec

csvwrite('gibbs_predictedSequence.csv',gibbs_float_vec)




%%%%%%%%%%%%%%%Estimation Curve PLOT%%%%%%%%%%
if ifplot
    figure(3)
    plot(index,radia_float_vec(index),index,decode_float_vec(index),index,em_decode_float_vec(index),index,gibbs_float_vec(index))
    legend('Raw Data','Likelihood Estimation','EM Estimation','Gibbs Estimation')

    title('Estimated Radiation Tensity Changing Curve')
    if every_day_analysis
        xlabel('Day index')
    else
        xlabel('Hour index')
    end
    ylabel('Radiation (W/m2)')
    hold on;
end




%%%%%%%%%%%ERROR PLOT%%%%%%%%%%%%%%
decode_float_vec = em_decode_float_vec;
% decode_float_vec = gibbs_float_vec

error = decode_float_vec-radia_vec;

if ifplot
    figure(4)
    plot(index,error)
    title('Error Curve')
    if every_day_analysis
        xlabel('Day index')
    else
        xlabel('Hour index')
    end
    ylabel('Radiation (W/m2)')
    ylim([-200 200])
    
end


RMSE = norm(decode_float_vec-radia_float_vec,2)/sqrt(length(decode_float_vec))
    

MAPE = mean(abs((decode_float_vec-radia_float_vec)./radia_float_vec)*100)
    
    
MABE = norm(decode_float_vec-radia_float_vec,1)/length(decode_float_vec)
    
r = corr2(decode_float_vec,radia_float_vec)




    
    
    
    
    
    
    
    
    
    



