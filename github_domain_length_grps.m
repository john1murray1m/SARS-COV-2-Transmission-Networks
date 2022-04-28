%% split seqs up into different groups determined by domain length

function ngrp_feas=github_domain_length_grps(seqs,grp_min)

ndomain=12; 
iframes=[2 1 400 10000 16000 13218 266;... % 1a
    1 12000 15000 4000 9000 8088 13768; ...% 1b
     2 17000 24000 3000 5000 3822 21563; ... % S % 
     1 24000 27000 600 1200 828 25393; ...% 3a
     1 25000 27000 130 350 228 26245; ...% E
     3 25500 27500 400 900 669 26523; ...% M
     1 26500 27700 100 280 186 27202;... % 6
     1 26500 28000 200 500 366 27394; ... % 7a
     3 26400 28500 90 200 132 27756; ... % 7b
     3 26500 28500 300 450 366 27894; ... % 8 
     2 27000 28500 750 1400 1260 28274; ...% N
     2 28000 30000 90 300 117 29558]; % 10

% don't want to include any of the sequences that have NaNs in Domains or
% have some start less than a previous one
iia=arrayfun(@(x) any(isnan(x.Domains(:))) || any(diff(x.Domains(:,1))<0), seqs);

% tabulate the different lengths
tab_gene_lengths=cell(ndomain,1);
for i=1:ndomain
    fred=arrayfun(@(x) x.Domains(i,2)-x.Domains(i,1),seqs);
     tb=tabulate(fred);
    ii1=tb(:,2)>0;
    tb=tb(ii1,:);
    tab_gene_lengths(i)={tb};
end

% for each of the ORFs determine the lengths with at least grp_min sequences and
% where the length is within the iframes bounds
tab_gene_lengths_grps=cell(ndomain,1);
for i=1:ndomain
    tb=tab_gene_lengths{i};
    ii1=(tb(:,2)>=grp_min & tb(:,1)>=iframes(i,4) & tb(:,1)<=iframes(i,5));
    tab_gene_lengths_grps(i)={tb(ii1,:)};
end

%% classify each sequence by the lengths of each domain - add in the length from ORF1a to ORF10 end
nseq=length(seqs);
domain_lengths_ind=NaN(nseq,ndomain+1);
for i=1:nseq
     domain_lengths_ind(i,:)=[seqs(i).Domains(:,2)-seqs(i).Domains(:,1); ...
                                seqs(i).Domains(ndomain,2)-seqs(i).Domains(1,1)];
end

% don't include the problem sequences
domain_lengths_ind(iia,:)=NaN;

[domain_grps,~,iidgrp]=unique(domain_lengths_ind,'rows');

% can separate into the different length groups
tab_grps=tabulate(iidgrp);

% determine groups of at least min size and within length bounds
ii1=(tab_grps(:,2)>=grp_min);
ii2=(sum(domain_grps(:,1:ndomain)>=repmat(iframes(:,4)',size(domain_grps,1),1),2)==ndomain);
ii3=(sum(domain_grps(:,1:ndomain)<=repmat(iframes(:,5)',size(domain_grps,1),1),2)==ndomain);
iigrp_feas=find(ii1&ii2&ii3);

ngrp_feas=cell(length(iigrp_feas),1); % the  sequence indices (within seqs subset in each group
for i=1:length(iigrp_feas)
    ngrp_feas(i)={find(iidgrp==iigrp_feas(i))};
end

