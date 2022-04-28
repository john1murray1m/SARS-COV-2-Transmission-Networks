%% determine the nt difference between groups 
function dn=github_nt_delta_edge1(seqs1,Cclone,dp)
% dn=sparse(size(dp,1),size(dp,1)); % the number of nt differences between the sequences in each group with the fewest nonACTG

% the first sequence in each clone will be the one with the fewest nonACTG
% sites
ii1=cellfun(@(x) x(1,1), Cclone); % the seq in each group used for the comparison
seqs=seqs1(ii1);

% only use the domain nt for comparison
ndomain=size(seqs(1).Domains,1);

for i=1:length(seqs)
    fred=seqs(i).Domains;
    george=[];
    for k=1:ndomain
        george=union(george,fred(k,1):fred(k,2));
    end
    len_seq=length(seqs(i).Sequence);
    if george(end)>len_seq % some orf did not have the final stop but were given the nt after the sequence end
        ii1=george<=len_seq;
        george=george(ii1);
    end
    seqs(i).Sequence=seqs(i).Sequence(george);
end
AAA=arrayfun(@(x) nt2int(x.Sequence), seqs,'UniformOutput',false);
AAA=cell2mat(AAA');
iiAAA=any(diff(AAA,1)~=0,1); % the posns that are not all the same

AAA1=AAA(:,iiAAA);
fred=sparse(AAA1>4 & AAA1<16); % the nonACTG posns

% only need the comparisons for those within a feasible day connection
     [jindex,kindex]=find(triu(dp+dp')>0);

% vector of index and dp and dpj values
jval=zeros(length(jindex),1); % a guess as to number of nonzero values
kval=jval;
dnval=jval;
icount=0;

for i=1:length(jindex)
    j=jindex(i);
    k=kindex(i);
    iicommon=(fred(j,:)==0 & fred(k,:)==0);
    dnval0=sum(AAA1(j,iicommon) ~= AAA1(k,iicommon))+1000; % have to add 1000 to differentiate those that
                                    % are the same distance (should not
                                    % happen with the clones) and those
                                    % where connection not feasible
    icount=icount+1;
    jval(icount)=j;
    kval(icount)=k;
    dnval(icount)=dnval0;

end

dn=sparse(jval(1:icount),kval(1:icount),dnval(1:icount),size(dp,1),size(dp,1));

% make the matrix symmetric
dn=(dn+dn');

% now get rid of the entries that are not possible from a date persepctive
ii3=dp(:)==0;
dn(ii3)=0;

