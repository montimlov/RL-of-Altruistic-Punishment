%% Data simulation, PR & MI script %% Adpated from script written by Pat Lockwood & Marco Wittmann 2019 by Jo Cutler 2020
%%% simulation script for parameter recovery and  model identifiability for different
%%% reinforcement learning models for prosocial learning task
%%%


clear all;
close all;
clearvars

addpath('models');
addpath('tools');

beep off;

models{1} = 'model_RWPL'; % enter model to run here **
models{2} = 'model_RWPL_temp_bias'; % enter model to run here **
models{3} = 'model_RWPL_temp_4beta'; % enter model to run here **

% Load in schedule
% -------------------------------------------- %
load trialorderPLZ0302.mat % specify trial order file here **
stim = unnamed;
nTrls = size(stim,1);
nBlocks = 1; % specify number of blocks here **
% Jialu's note: 一些全局参数
nagent = 2;
mod_ID = 2; %当前检验的模型编号

alphabounds = [0, 1]; % enter bounds on alpha values here **
betabounds = [0, 1.803]; % enter bounds on beta values here **
% Jialu's note: 不知道为啥betabounds的边界是0.3
biasbounds = [0,0.5];

alphamin = min(alphabounds);
alphamax = max(alphabounds);
betamin = min(betabounds);
betamax = max(betabounds);
biasmin = min(biasbounds);
biasmax = max(biasbounds);
rng default % resets the randomisation seed to ensure results are reproducible (MATLAB 2019b)
% Jialu's note: 在这个地方改要检验什么模型(对应上边28-30行models123那个地方)
modelsTR = mod_ID; % enter the model number to run PR on here - numerical index in models variable


% Jialu's note: 模拟的人数
nSubj = 10000; % how many subjects to simulate if not defined by grid **

for m = modelsTR % loop over model number(s) specified above
    
    clearvars -except models* *bounds *min *max stim nTrls nBlocks nagent nSubj tr nRounds mle* r m all_*
    
    modelID = models{m};
    s.PL.expname = 'ProsocialLearn';
    % Jialu's note: 根据模型选择参数，这里需要改
    switch modelID % if adding new model functions also add the parameters here **
        case {'model_RWPL'}
            params  = {'alpha','beta'};
        case {'model_RWPL_temp'}
            params  = {'alpha_self_front', 'alpha_other_front','alpha_self_later', 'alpha_other_later','beta_self','beta_other'};
        case {'model_RWPL_temp_bias'}
            params  = {'alpha_self_front', 'alpha_other_front','alpha_self_later', 'alpha_other_later','beta_self','beta_other','bias'};
        case {'model_RWPL_temp_4beta'}
            params  = {'alpha_self_front', 'alpha_other_front','alpha_self_later', 'alpha_other_later',...
                'beta_self_front', 'beta_other_front','beta_self_later', 'beta_other_later'};
        otherwise
            error(['No parameters defined for model ', modelID, '. Check modelID parameter'])
    end
    
    modelIDs = [modelID, '_simulate']; % function name to use for a model that simulates the dta
    modelIDr = [modelID, '_real']; % function name to use for a model to fit the simulated data
    nParam  = length(params);
    
    msg = ['Estimating the optimal alpha for: ', modelID, ', calculating ', num2str(nParam), ' parameters: ', char(params{1})];
    for n = 2:nParam
        msg = [msg, ', ', char(params{n})];
    end
    disp(msg) % show details in command window
    
    % Set parameters to simulate
    % -------------------------------------------- %
    
   
     % Jialu's note:在这里生成所有的参数组合，相当于随机参数赋值给不同的模拟被试
    for sub = 1:nSubj
        for ip = 1:length(params)
            thisp=params{ip};
            if contains(thisp, 'alpha') == 1
                allCombs(sub,ip) = betarnd(1.1,1.1); % define distribution **
                while norm2alpha(alpha2norm(allCombs(sub,ip))) > (alphamax) || norm2alpha(alpha2norm(allCombs(sub,ip))) < (alphamin)
                    allCombs(sub,ip) = betarnd(1.1,1.1);
                end
            elseif contains(thisp, 'beta') == 1
                allCombs(sub,ip) = gamrnd(1.2,5); % define distribution **
                while norm2beta(beta2norm(allCombs(sub,ip))) > (betamax) || norm2beta(beta2norm(allCombs(sub,ip))) < (betamin)
                    allCombs(sub,ip) = gamrnd(1.2,5);
                end
            elseif contains(thisp, 'bias') == 1
                allCombs(sub,ip) = betarnd(1.1,1.1); % define distribution **
                while norm2bias(bias2norm(allCombs(sub,ip))) > (biasmax) || norm2bias(bias2norm(allCombs(sub,ip))) < (biasmin)
                    allCombs(sub,ip) = betarnd(1.1,1.1);
                end
            else
                error('Define param as one of above cases');
            end
        end
    end
    noise = 0;
    
    % Transform parameters
    % -------------------------------------------- %
    
    allCombsNorm = NaN(size(allCombs,1),size(allCombs,2));
    for ip=1:length(params)
        thisp=params{ip};
        if contains(thisp, 'alpha') == 1
            allCombsNorm(:,ip)=alpha2norm(allCombs(:,ip));
        elseif contains(thisp, 'beta') == 1
            allCombsNorm(:,ip)=beta2norm(allCombs(:,ip));
        elseif contains(thisp, 'bias') == 1
            allCombsNorm(:,ip)=bias2norm(allCombs(:,ip));
        else
            error(['Can`t detect whether parameter ', thisp, ' is alpha or beta']);
        end
    end
    
    % Plot the distribution of parameters we are using to simulate behaviour
    % -------------------------------------------- %

    figure('color','w');
    for param=1:nParam
        subplot(1,nParam,param);
        thisp=params{param};
        if contains(thisp, 'alpha') == 1
            histogram(norm2alpha(allCombsNorm(:,param)+noise*randn(length(allCombsNorm),1)),'FaceColor',[0.5 0.5 0.5]);
        elseif contains(thisp, 'beta') == 1
            histogram(norm2beta(allCombsNorm(:,param)+noise*randn(length(allCombsNorm),1)),'FaceColor',[0.5 0.5 0.5]);
        elseif contains(thisp, 'bias') == 1
            histogram(norm2bias(allCombsNorm(:,param)+noise*randn(length(allCombsNorm),1)),'FaceColor',[0.5 0.5 0.5]);
        end
        hold on;box off;title(params{param});
    end
    
    % Simulate all combinations of parameters
    % -------------------------------------------- %
    
    allChoices = nan(nSubj, nagent);
    allOutcomes = nan(nSubj, nagent);
    
    for simS=1:nSubj
        allbl=[];
        if simS < 0.5*nSubj
            env = ones(60,1); %环境鼓励惩罚
        else
            env = 2* ones(60,1); %环境鼓励不惩罚
        end
        
        % Jialu's add: 随机生成PE序列
        PEs = [];
        for iPE = 1:60
            randitem = rand;
            if randitem < 0.8
                PEt = 1;
            else
                PEt = 0;
            end
            PEs = [PEs;PEt];
        end
        stim(:,2)=PEs;
        Data(simS).ID        = sprintf('Subj %i',simS);
        Data(simS).data      = nan(nTrls,1); % save choices
        Data(simS).agent     = stim(:,1);
        Data(simS).normalPE  = stim(:,2);
        Data(simS).block     = env;
        Data(simS).trueModel = modelIDr;
        
        for param=1:nParam % Add some noise to grid parameters
            thisp=params{param};
            if contains(thisp, 'alpha') == 1
                pmin = alpha2norm(alphamin);
                pmax = alpha2norm(alphamax);
            elseif contains(thisp, 'beta') == 1
                pmin = beta2norm(betamin);
                pmax = beta2norm(betamax);
            elseif contains(thisp, 'bias') == 1
                pmin = bias2norm(biasmin);
                pmax = bias2norm(biasmax);
            end
            truep(param) = allCombsNorm(simS,param) + noise*randn(1);
            while truep(param) < pmin || truep(param) > pmax
                truep(param) = allCombsNorm(simS,param) + noise*randn(1);
            end
        end
        
        Data(simS).trueParam = truep;
        
        simfunc = str2func(modelIDs);
        
        try% Jialu's note: 这里是拟合函数的入口
            [f,allout] = simfunc(Data(simS).data, Data(simS).normalPE, Data(simS).block, Data(simS).agent, Data(simS).trueParam, alphabounds, betabounds); %%%% VARIABLES TO FEED INTO THE MODEL IN ORDER TO SIMULATE CHOICES
        catch
            disp('Error in call to simulate - possibly argument "allout" (and maybe others) not assigned')
        end
        
        s.PL.beh{1,simS}.choice = allout.all_data';
        s.PL.beh{1,simS}.agent = allout.all_agent';
        s.PL.beh{1,simS}.outcome = allout.all_outcome';
        s.PL.beh{1,simS}.block = allout.all_block;
        s.PL.ID{1,simS}.ID = simS;
        
        for a = 1:nagent% Jialu's note: 这里是针对每个agent进行选择
            agentTrials = find(s.PL.beh{1,simS}.agent == a);
            choicePer = sum(s.PL.beh{1,simS}.choice(agentTrials)) / length(agentTrials);
            outcomePer = sum(s.PL.beh{1,simS}.outcome(agentTrials)) / length(agentTrials);
            allChoices(simS, a) = choicePer;
            allOutcomes(simS, a) = outcomePer;
        end
        
    end
    
end

% Jialu's note: 整理数据成画图的格式，现在是合并内外群的形式
allComb = [allCombs(:,5); allCombs(:,6)];
allChoice = [allChoices(:,1); allChoices(:,2)];
allOutcome = [allOutcomes(:,1); allOutcomes(:,2)];

figure(2);
scatter(allComb, allOutcome)
%figure(3);
%scatter(allComb, allOutcome)

optimaldata = [allComb, allChoice, allOutcome];

ODtab = cell2table(num2cell(optimaldata), 'VariableNames', {'Alpha', 'Choices', 'Outcomes'});
ODname = ['../Prosocial_learning_R_code/Optimal_alpha.xlsx'];
writetable(ODtab,ODname,'WriteVariableNames',true)

% for a = 1:3
%     scatter(allCombs(:,a), allChoices(:,a))
%     scatter(allCombs(:,a), allOutcomes(:,a))
% end