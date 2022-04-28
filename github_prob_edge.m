%% calculate the probability a clone transmitted to another based on dates of sequence submission
function [dp,dpj]=github_prob_edge(Cclone,Dates)
    
nclones=length(Cclone);
% dp=sparse(nclones,nclones); % the prob clone j transmitted to clone k in dp(j,k). 
%                                            % the node in k is always
%                                            % the root.
% dpj=sparse(nclones,nclones); % the index j twithin the Cclone{j} sequences that gives the max prob

Cclone_dates=cellfun(@(x) (Dates(x)), Cclone,'UniformOutput',false);
[Cclone_root,Cclone_iroot]=cellfun(@(x) min(x), Cclone_dates); % the root date and the index within the clone


%% set the parameters and shift for the gamma distribution from McAloon et al., BMC Public Health, 2021, for probability of transmission between individuals
gk=-10.5;
gshape=21.96;
gscale=1.57; 
pd=makedist('Gamma','a',gshape,'b',1/gscale);

minday=-5; % change the min days preceeding another. Otherwise the penalty is very high and creates too many additional edges. 

% vector of index and dp and dpj values
jval=zeros(round(nclones^1.8),1); % a guess as to number of nonzero values
kval=jval;
dpval=jval;
dpjval=jval;
icount=0;

for j=1:nclones
    for k=1:nclones
        if k~=j
            fred=Cclone_root(k)-Cclone_dates{j};
            if any(fred>= minday & fred<=17)
                [dpval0,dpjval0]=max(pdf(pd,fred-gk));
                 icount=icount+1;
                jval(icount)=j;
                kval(icount)=k;
                dpval(icount)=dpval0;
                dpjval(icount)=dpjval0;
               
            end
        end
    end
end

dp=sparse(jval(1:icount),kval(1:icount),dpval(1:icount),nclones,nclones);
dpj=sparse(jval(1:icount),kval(1:icount),dpjval(1:icount),nclones,nclones);
