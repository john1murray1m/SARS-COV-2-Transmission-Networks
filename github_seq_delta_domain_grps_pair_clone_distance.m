%% calculate the distance and date distance matrices for the lineage/length groups

% 
% % save('Delta_seq_data_all','cov_delta_meta','no_seq_delta_all',...
% %     'seqs_delta','seq_ref','save_codes')
% load('Delta_seq_data_all')
% 
% % save('Delta_seq_domains','dom_align_calcs','domain_refs',"ndomain",'seqs_delta')
% load('Delta_seq_domains')

% %     save('Delta_seqs_domain_orfs',"seqs_delta",'seq_group_probs')
%     load('Delta_seqs_domain_orfs')

% % save('Delta_lineage_domain_groups',"ngrps_delta_lin",'delta_lineages','grp_min')
% load('Delta_lineage_domain_groups')

 Dist_delta_clone_grps=cell(size(ngrps_delta_lin,1),1); % within each remaining lineage group 
        % 1: tab_s tabulation of the number nonACTG posns ('N' etc) in the sequences
        % 2: [mm1, mm2, mm3] mm1 is the max number of nonACTG posns in a (mm2) sequence
        %  mm3 is min
        % 3: iinonACTG the sequences includes in the dist calcs. The others have
        % been omitted with problems with the domain
        % calcs
        % 4: clone groups
        % 5: the number of clones grouped and the number not grouped so far 
  for i=1:size(ngrps_delta_lin,1)
    seqs0=seqs_delta(ngrps_delta_lin{i,2}); 
    % now get any subgroups specified by domain lengths
    sally=ngrps_delta_lin{i,3};
    Dist_grps=cell(length(sally),5);
    for j=1:length(sally)
        seqs=seqs0(sally{j});

        [tab_s,mm1,mm2,mm3,iikeep,Cclone,tot_clones,missing_seq]=github_ACTG_delta_clones_recur(seqs);
        Dist_grps(j,1)={tab_s};
        Dist_grps(j,2)={[mm1,mm2,mm3]};
        Dist_grps(j,3)={iikeep}; % the sequences with nonNAN dates
        Dist_grps(j,4)={Cclone}; % the indexing is relative to the original set of sequences (seqs) pre limiting by iikeep and reordering
        Dist_grps(j,5)={[tot_clones,missing_seq]};
    end
    Dist_delta_clone_grps(i)={Dist_grps};
  end

   save('seq_delta_pair_clone_distance','Dist_delta_clone_grps')
