%% determine the delta lineage and domain length groups of size at least grp_min
 
% % save('Delta_seq_domains','dom_align_calcs','domain_refs',"ndomain",'seqs_delta')
% load('Delta_seq_domains')

% %     save('Delta_seqs_domain_orfs',"seqs_delta",'seq_group_probs')
%     load('Delta_seqs_domain_orfs')

lineages=[seqs_delta.Lineage];

tab_lin=tabulate(lineages);

grp_min=20;
tab_lin_num=cell2mat(tab_lin(:,2));
iij=find(tab_lin_num>=grp_min);

delta_lineages=tab_lin(iij,1);
ngrps_delta_lin=cell(length(delta_lineages),3); % this will be the outer grouping of lineage grps
                % within these will be the separate groups determined by
                % domain lengths. % the 2nd value is the set of sequence
                % indices within seqs_delta
                % and the 1st is the lineage name
                % the 3rd will contain cells with the seq indices split into different groups determined by domain length 
 for ilin=1:length(delta_lineages)
    ngrps_delta_lin(ilin,1)=delta_lineages(ilin);
    ii1=strcmp([seqs_delta.Lineage],delta_lineages(ilin));
    seqs=seqs_delta(ii1);
    ngrps_delta_lin(ilin,2)={find(ii1)};

    % now split into equivalent domain lengths
    ngrp_feas=github_domain_length_grps(seqs,grp_min);
    ngrps_delta_lin(ilin,3)={ngrp_feas};
end
% remove the empty ones
len_grps=cellfun(@(x) length(x),ngrps_delta_lin(:,3));
ii5=len_grps==0;
ngrps_delta_lin(ii5,:)=[];

save('Delta_lineage_domain_groups',"ngrps_delta_lin",'delta_lineages','grp_min')