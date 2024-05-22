function [fval,fit] = mod_ms_RWPL_temp_bias(behavData,q, doprior,dofit,varargin)
%%
% -------------------------------------------------------------------------------------
% 1 ) Define free parameters
% -------------------------------------------------------------------------------------

if nargin > 4
    prior      = varargin{1};
end

qt = norm2par('ms_RWPL_temp_bias',q); % transform parameters from gaussian space to model space

% Define free parameters and set unused ones to zero
beta_self   = qt(1);
beta_other  = qt(2);
alpha_self_front  = qt(3);
alpha_other_front = qt(4);   
alpha_self_later  = qt(5);
alpha_other_later = qt(6);   
bias = behavData.bias;

allbetas=[beta_self  beta_other];

%num_blocks=1;
num_cond=2;
num_pic_reps=30;

all_prob = [];
%

all_Va  = [];
all_Vb  = [];
all_PE  = [];

% if (beta_self<0.05), fval=10000000; return; end;
% if (beta_other<0.05), fval=10000000; return; end;
% if (beta_noone<0.05), fval=10000000; return; end;


% Change missed trials coded as '9' to NaNs since other parts of the
% pipeline presume missed trials are NaNs

%%% 0.) Load information for that subject:    % load in each subjects variables for the experiment
choice  = behavData.choice; %matrix of choices for self (stimulus in columns, repetitions in rows)
outcome = behavData.outcome; %matrix of outcomes for self(stimulus in columns, repetitions in rows)
agent   = behavData.agent;
group = behavData.block;
startVa   = behavData.startvalueA;
startVb   = behavData.startvalueB;

choice(choice   == 9) = nan;
outcome(outcome == 9) = nan;
agent(agent     == 9) = nan;


for iagent=1:num_cond; %%%% loop through the 3 conditions, self, other, no one
    

    
    % pick values for this sequence
    choicet = choice(agent==iagent);                                                % choice =1 for high prob, 0 for low prob
    outcomet= outcome(agent==iagent);
    agentt =  agent(agent==iagent);% outcome= 1 for correct, 0 for incorrect
    startVat   = startVa;
    startVbt   = startVb;
    groupt = group(1);

    % make empty matrix for values you want to collect:
    probs_choice =nan(num_pic_reps,1);
    Va           =nan(num_pic_reps,1);      %V_self(1) =  start_bias(1);
    Vb           =nan(num_pic_reps,1);      %V_fri (1) =  start_bias(2);
    PE           =nan(num_pic_reps,1);
    
    Va(1)=startVat(1);
    Vb(1)=startVbt(1);

    for t=1:num_pic_reps  % for the 16 reps of the stimuli                                               % loop through every presentation of the same stimulus
        
        %%% 1. DECISION/CUE PHASE
        
        % first calculate parts of softmax:
        if agentt(t) == 1 %ingroup add bias
            if groupt == 1 % encorage punish, Va - punish, Vb add bias
                if choicet(t)  ==1, softmax_num=exp(Va(t)/allbetas(iagent)); end
                if choicet(t)  ==0, softmax_num=exp((Vb(t)+bias)/allbetas(iagent)); end
                softmax_denom = (exp(Va(t)/allbetas(iagent)) +  exp((Vb(t)+bias)/allbetas(iagent)));
            elseif groupt == 2 % encorage unpunish, Va - unpunish, Va add bias
                if choicet(t)  == 1, softmax_num=exp((Va(t)+bias)/allbetas(iagent)); end
                if choicet(t)  == 0, softmax_num=exp(Vb(t)/allbetas(iagent)); end
                softmax_denom = (exp((Va(t)+bias)/allbetas(iagent)) +  exp(Vb(t)/allbetas(iagent)));
            end
        elseif agentt(t) == 2 % outgroup no change
            if choicet(t)  ==1, softmax_num=exp(Va(t)/allbetas(iagent)); end
            if choicet(t)  ==0, softmax_num=exp(Vb(t)/allbetas(iagent)); end
            softmax_denom =    (exp(Va(t)/allbetas(iagent)) +  exp(Vb(t)/allbetas(iagent)));
        end
        % apply softmax to get choice probability for current trial:
        
        if ~isnan(choicet(t)), probs_choice(t)=softmax_num/softmax_denom; end
        
        %%% 2. FB PHASE
        
        % set defaults for (next) trial, then calculate PEs and update:
        Va(t+1)    = Va(t);                                       % keep value of self,other no ones same at t+1
        Vb(t+1)    = Vb(t);
        PE(t)      = nan;
        if t <= 15
            if agentt(t)==1
                
                if choicet(t)== 1, PE(t)  = outcomet(t) - Va(t);          Va(t+1) = Va(t) +  (alpha_self_front*PE(t));      end;
                if choicet(t)== 0, PE(t)  = outcomet(t) - Vb(t);          Vb(t+1) = Vb(t) +  (alpha_self_front*PE(t));       end;
            
            elseif agentt(t) ==2 %|| agentt(t)==3; for combined alpha
                
                if choicet(t)== 1, PE(t)  = outcomet(t) - Va(t);          Va(t+1) = Va(t) +  (alpha_other_front*PE(t));      end;
                if choicet(t)== 0, PE(t)  = outcomet(t) - Vb(t);          Vb(t+1) = Vb(t) +  (alpha_other_front*PE(t));       end;
                   
            end
        else
            if agentt(t)==1
                
                if choicet(t)== 1, PE(t)  = outcomet(t) - Va(t);          Va(t+1) = Va(t) +  (alpha_self_later*PE(t));      end;
                if choicet(t)== 0, PE(t)  = outcomet(t) - Vb(t);          Vb(t+1) = Vb(t) +  (alpha_self_later*PE(t));       end;
            
            elseif agentt(t) ==2 %|| agentt(t)==3; for combined alpha
                
                if choicet(t)== 1, PE(t)  = outcomet(t) - Va(t);          Va(t+1) = Va(t) +  (alpha_other_later*PE(t));      end;
                if choicet(t)== 0, PE(t)  = outcomet(t) - Vb(t);          Vb(t+1) = Vb(t) +  (alpha_other_later*PE(t));       end;
                   
            end
        end
        if isnan(choicet(t)),  probs_choice(t) = NaN;  end;
        
    end % number of pics rep 16 times each stim is presented
    
    %%% 4. now save stuff:
    all_Va     =  [all_Va Va(1:num_pic_reps)];
    all_Vb     =  [all_Vb Vb(1:num_pic_reps)];
    all_PE     =  [all_PE PE];
    all_prob   =  [all_prob; probs_choice];
    

    
end % cond



% all choice probablities
ChoiceProb=all_prob';


% -------------------------------------------------------------------------------------
% 4 ) Calculate model fit:  
% -------------------------------------------------------------------------------------

nll =-nansum(log(ChoiceProb));                                                % the thing to minimize                      

if doprior == 0                                                               % NLL fit
   fval = nll;
elseif doprior == 1                                                           % EM-fit:   P(Choices | h) * P(h | O) should be maximised, therefore same as minimizing it with negative sign   
   fval = -(-nll + prior.logpdf(q));
end

% % make sure f is not just low because of Nans in prob-variable:

sumofnans=sum(sum(isnan(choice)));
if sum(isnan(ChoiceProb))~=sumofnans disp('ERROR NaNs in choice and choice prob dont agree'); keyboard; return; end               

% -------------------------------------------------------------------------------------
% 5) Calculate additional Parameters and save: 
% -------------------------------------------------------------------------------------

if dofit ==1
   %vsum = o1_val + o2_val ;
   
   fit         = struct;
   fit.xnames  = {'beta_self'; 'beta_other';'alpha_self_front'; 'alpha_other_front';'alpha_self_later';'alpha_other_later';'bias'};
   
   fit.choiceprob = [ChoiceProb];% %%% NEW as choice prob 400* 1 and values stored as 25*16
   fit.mat    = [all_Va all_Vb all_PE ];
   fit.names  = {'Va';'Vb';'all_PE';};
end






