%% calculate the adjacency matrix for clone groups

function [A, bingrps]=github_CLE_adjacency_clone(dp,dn,maxnt_dist,min_bin_size)

% the adjacency matrix 
A=sparse(size(dp,1),size(dp,1)); 

ii1=(dn(:)-1000)<=maxnt_dist & dp(:)>0; % only allow connections if nt dist <= maxnt_dist and time feasbile

PP=-log(dp(ii1)); % converts prob to higher pos numbers for less likely connection

A(ii1)=PP; 
A=A-diag(diag(A)); % remove self connections and will have edges between same sequences
% have to prioritise connections within the same set of mutations. Make
% connections between groups separated by k nt diffs, penalised by
% k*mut_pen. That way preferences will be for connections via sequential
% single muts.
mut_pen=10;

for k=1:maxnt_dist
    ii3=(mod(dn(:),10)==k & A(:)>0); % the nt dist calc has 1000 added to it for feasible connections so remove these
    if sum(ii3)>0
        A(ii3)=A(ii3)+k*mut_pen;
    end
end

G=digraph(A);

% get the connected components
[bins,binsizes]=conncomp(G,'Type','weak');

ii2=find(binsizes>=min_bin_size);
if ~isempty(ii2)
    ii3=ismember(bins,ii2);
    bingrps=bins.*ii3; % this will be nonzero only for the ones in a weakly connected group of size >=min_bin_size
                    % the number will represent the group that seq belongs
                    % to.
else
    bingrps=zeros(size(bins));
end

end
