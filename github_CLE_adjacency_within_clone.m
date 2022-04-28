%% calculate the adjacency matrix for sequences with a clone 


function A=github_CLE_adjacency_within_clone(seqs,dates)


%% set the parameters and shift for the gamma distribution from McAloon for probability of transmission between individuals
gk=-10.5;
gshape=21.96;
gscale=1.57; 
pd=makedist('Gamma','a',gshape,'b',1/gscale);

minday=-5; % change the min days preceeding from 

%%
% the adjacency matrix 
A=sparse(length(dates),length(dates)); 

% date differences
A0=zeros(length(dates),length(dates));

for i=1:length(dates)
    A0(i,:)=dates'-dates(i);
end

ii1=(A0(:)>=minday & A0(:)<=17); % within the feasible time frame

PP=pdf(pd,A0(ii1)-gk); % the probability of the feasible times

A(ii1)=-log(PP); 
A=A-diag(diag(A)); % remove self connections 


