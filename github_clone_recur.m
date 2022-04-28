%% a recursive function to separate into clones
function [Cclone,iclone,Cdist,iid1]=github_clone_recur(iid,Cdist0,iclone,Cclone)

% input
%   1) iid, the sequence indices relative to the original sequence set
%   2) Cdist0, the groups of sequences equidistant from the first
%   3) iclone, the current clone count,
%   4) george, the nonACTG matrix for these sequences

% output
%   1) Cclone, the new clone groups found
%   2) iclone, the updated clone counter
%   3) Cdist, the new groups of equidistant seqs


if ~isempty(Cdist0)    
    % 2nd level
    for ii1=1:size(Cdist0,1) % the number of groups still needing separation
        iid01=iid(Cdist0{ii1,1});
        seqs=Cdist0{ii1,2};
        [Cgood,Cdist1]=github_clone_sep(seqs,Cdist0{ii1,3});
        
        for ii2=1:length(Cgood)
            Cclone(iclone)={iid01(Cgood{ii2})};
            iclone=iclone+1;
        end
        [Cclone,iclone,Cdist,iid1]=github_clone_recur(iid01,Cdist1,iclone,Cclone);
    end
else
    Cdist={};
    iid1=[];
end
