%% for a given set of sequences, calculated pairwise length vector dd, and separate into
% input:    1) seq, the set of sequences ordered in increasing nonACTG values
%           2) george, the nonACTG submatrix for these sequences
% output:
% Cgood = cell of clones (zero distance to the root or singletons) - no
% further separation is required for these. These are listed as the set of
% seq indices within seq. Have to work back to get the original indices

% Cdist = cell of 1) indices and 2) the sequences and 3)
% the nonACTG submatrix of those the same
% distance from the root. This will be empty if this path is finished

function [Cgood, Cdist]=github_clone_sep(seqs,george)

Cgood={};
Cdist={};

% specify the sequence indices
ii00=(1:length(seqs));

% pairwise distance from the first, with the first value 0
dd=pair_distance(seqs,george);

% igrps the sequence numbers in each grp
igrps=unique(dd); % thedifferent distance grps to the first sequence

igood=1; % the counter of finalised groups
idist=1; % the counter for nonfinalised groups

for idd=1:length(igrps)
    % get the sequences with these distances from the sequence with fewest
    % nonACTG posns
    kd=igrps(idd);
    ii0=(dd==kd);
    if kd==0
            Cgood(igood)={ii00(ii0)}; % the indices within seq. The first seq will alway have 0 dist to itself
            igood=igood+1;            
    elseif kd>0 && sum(ii0)==1
        Cgood(igood)={ii00(ii0)};
        igood=igood+1;    
    else
        Cdist(idist,1)={ii00(ii0)};
        Cdist(idist,2)={seqs(ii0)};
        Cdist(idist,3)={george(ii00(ii0),:)};
        idist=idist+1;
    end
end


end
