%% determine pairwise distances and group into clone sets for delta

function [tab_s,mm1,mm2,mm3,iikeep,Cclone1,tot_clones,missing_seq]=github_ACTG_delta_clones_recur(seqs1)


% remove dates that are NaN
iikeep=find(~isnan(datenum([seqs1.Date])));
% only keep the nonnan date sequences since this will be used to connect -
seqs=seqs1(iikeep);

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


    
% the bad positions for the group of sequences. Produces a 0-1 matrix
% of (no. seqs)*(length of each sequence) with 1 where there is a
% nonACTG nt value
numnonnactg=arrayfun(@(x) nt2int(x.Sequence)>4 & nt2int(x.Sequence)<16, seqs,'UniformOutput',false);
fred=cell2mat(numnonnactg');
fred=sparse(fred);
len_seq=size(fred,2);

% reorder sequences in increasing order of problem positions so pairwise
% distances will be done relative to the most complete sequences
[~,iiorder]=sort(sum(fred,2));
fred=fred(iiorder,:);
seqs=seqs(iiorder);

% the sequence with the largest number of nonACTG
[mm1,mm2]=max(sum(fred,2));
 mm3=min(sum(fred,2));


% how big are the problems in each sequence
tab_s=tabulate(sum(fred,2));
ii1=tab_s(:,2)>0;
tab_s=tab_s(ii1,:);
%    if you want to plot the number sequences with values of nonACTG positions  
%     figure(201)
%     clf
%     subplot(2,2,1)
%     plot(tab_s(:,1),tab_s(:,2))
%     xlabel('Sequences with x nonACTG')
%     ylabel('Number of sequences')
%     
%     subplot(2,2,2)
%     survnum=sum(tab_s(:,2))-[0; cumsum(tab_s(1:(end-1),2))];
%     semilogy(tab_s(:,1),survnum)
%     xlabel('Sequences with >=x nonACTG')
%     ylabel('Number of sequences')

%% calculate the separate groups for this set of sequences

iclone=1;
% specify the sequence indices
ii00=(1:length(seqs));

% 1st level
[Cgood, Cdist]=github_clone_sep(seqs,fred);
for ii0=1:length(Cgood)
    Cclone(iclone)={ii00(Cgood{ii0})};
    iclone=iclone+1;
end
% the index of sequences relative to the original indices
iid=ii00;

[Cclone,iclone,Cdist1,iid1]=clone_recur(iid,Cdist,iclone,Cclone);


% put back into correct order according to the whole seq set pre-order
Cclone1=cellfun(@(x) iikeep(iiorder(x)), Cclone, 'UniformOutput', false);

% how many seq included in the groups
num_clones=cellfun(@(x) length(x), Cclone);
tot_clones=sum(num_clones);

missing_seq=length(seqs)-tot_clones;