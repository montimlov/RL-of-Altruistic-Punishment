% M4 = 'RWPL_SON_LR_SON_beta'
% original model from PNAS with 3 learning rates and 3 betas - beta self, beta other, beta no one, self LR, other LR, no one LR

function [f, allout] = model_RWPL_temp_bias_simulate(data,normalPE,block,agent,p,alphabounds,betabounds)

%%%%% 1. Assign free parameters and other stuff:

% use norm2alpha and norm2beta to get them in a sensible range

alpha_self_front  = norm2alpha(p(1));
alpha_other_front = norm2alpha(p(2));
alpha_self_later = norm2alpha(p(3));
alpha_other_later = norm2alpha(p(4));
beta_self   = norm2beta(p(5));
beta_other   = norm2beta(p(6));
bias = norm2bias(p(7));

all_alphas_front = [alpha_self_front, alpha_other_front];
all_alphas_later = [alpha_self_later; alpha_other_later];
allbetas=[beta_self; beta_other];

alphamin = min(alphabounds);
alphamax = max(alphabounds);
betamin = min(betabounds);
betamax = max(betabounds);
biasmin = 0;
biasmax = 0.5;

num_blocks=1;
num_cond=2;
num_pic_reps=30;

all_prob = [];
all_Va  = [];
all_Vb  = [];
all_PE  = [];

all_data    = [];
all_outcome = [];
all_agent   = [];
all_block = [];

%%% constrain learning rates:
if (alpha_self_front<alphamin || alpha_self_front>alphamax), f=10000000; return; end;
if (alpha_other_front<alphamin || alpha_other_front>alphamax), f=10000000; return; end;
if (alpha_self_later<alphamin || alpha_self_later>alphamax), f=10000000; return; end;
if (alpha_other_later<alphamin || alpha_other_later>alphamax), f=10000000; return; end;
if (beta_self<betamin || beta_self>betamax), f=10000000; return; end;
if (beta_other<betamin || beta_other>betamax), f=10000000; return; end;
if (bias<biasmin || bias>biasmax), f=10000000; return; end;

for u=1:num_blocks
    
    for group=1:num_cond %%%% loop through the 3 conditions, self, other, no one
       
        % pick values for this sequence
        normalPEt= normalPE(agent==group); % normalPE = 1 for normal, 0 for PE trial (correct not rewarded, incorrect rewarded)
        agentt =  agent(agent==group);
        env = block(1);
        
        % make empty matrix for values you want to collect:
        probs_choice =nan(num_pic_reps,1);
        Va           =nan(num_pic_reps,1);
        Vb           =nan(num_pic_reps,1);
        PE           =nan(num_pic_reps,1);
        outcomet     =nan(num_pic_reps,1);
        choicet      =nan(num_pic_reps,1);
        
        [Va0,Vb0] = getV0();
        Va(1)=Va0;
        Vb(1)=Vb0;
        
        for t=1:num_pic_reps                                               % loop through every presentation of the same stimulus
            
            %%% 1. DECISION/CUE PHASE
            if agentt(t) == 1 %ingroup add bias
                if env == 1 % encorage punish, Va - punish, Vb add bias
                    probs_choice(t) = exp(Va(t)/allbetas(group))/(exp(Va(t)/allbetas(group))+exp(bias+Vb(t)/allbetas(group)));
                    choicet(t) = double(rand < (probs_choice(t)));
                elseif env == 2 % encorage unpunish, Va - unpunish, Va add bias
                    probs_choice(t) = exp(bias+Va(t)/allbetas(group))/(exp(bias+Va(t)/allbetas(group))+exp(Vb(t)/allbetas(group)));
                    choicet(t) = double(rand < (probs_choice(t)));
                end
            elseif agentt(t) == 2 % outgroup no change
                probs_choice(t) = exp(Va(t)/allbetas(group))/(exp(Va(t)/allbetas(group))+exp(Vb(t)/allbetas(group)));
                choicet(t) = double(rand < (probs_choice(t)));
            end
            
            %%% 2. FB PHASE
            
            % set defaults for (next) trial, then calculate PEs and update:
            Va(t+1)    = Va(t);                                       % keep value of self,fri,str same at t+1
            Vb(t+1)    = Vb(t);
            PE(t)      = nan;
            
            if choicet(t)== 1 % if choose the high reward option
                if normalPEt(t) == 1 % if the trial is normal so high reward is rewarded
                    outcomet(t) = 1;
                elseif normalPEt(t) == 0
                    outcomet(t) = 0;
                end
                PE(t)  = outcomet(t) - Va(t);
                if t<=15
                    Va(t+1) = Va(t) +  (all_alphas_front(group)*PE(t));
                else
                    Va(t+1) = Va(t) +  (all_alphas_later(group)*PE(t));
                end
            end
            
            if choicet(t)== 0 % if choose the low reward option
                if normalPEt(t) == 1 % if the trial is normal so low reward is not rewarded
                    outcomet(t) = 0;
                elseif normalPEt(t) == 0
                    outcomet(t) = 1;
                end
                PE(t)  = outcomet(t) - Vb(t);
                if t<=15
                    Vb(t+1) = Vb(t) +  (all_alphas_front(group)*PE(t));
                else
                    Vb(t+1) = Vb(t) +  (all_alphas_later(group)*PE(t));
                end
            end
            
        end % number of pics rep 16 times each stim is presented
        
        %%% 4. now save stuff:
        all_Va     =  [all_Va Va(1:num_pic_reps)];
        all_Vb     =  [all_Vb Vb(1:num_pic_reps)];
        all_PE     =  [all_PE PE];
        all_prob   =  [all_prob; probs_choice];
        all_data    =  [all_data choicet'];
        all_outcome =  [all_outcome outcomet'];
        all_agent   =  [all_agent agentt'] ;
        
    end % cond
    
end % block

% all choice probablities
prob=all_prob';

% make sure f is not just low because of Nans in prob-variable:
sumofnans=sum(sum(isnan(all_data)));  % calculate number of "too slow" responses, as they are set to "nan" in the "prob"-variable
if sum(isnan(prob))~=sumofnans, f=10000000; return; end                    % choice prob can only be nan if it was a "too slow" response, not otherwise

% now calculate f:
f=-nansum(log(prob));
    
allout.all_Va       = all_Va;
allout.all_Vb       = all_Vb;
allout.all_PE       = all_PE;
allout.prob         = prob;
allout.nll          = f;
allout.all_data     = all_data ;
allout.all_outcome  = all_outcome;
allout.all_agent    = all_agent;
allout.all_block    = block';
f=allout;

end