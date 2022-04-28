%% determine the differences between the root nodes of each clonal network for a all delta lineage

% load('Delta_seqs_domain_orfs')
% load('seq_delta_pair_clone_distance')
% load('save_delta_pair_clone_Gopt')
% choose the lineage
% choose the larger ngrps_delta_lin 1:8
icount=0;
clear idx_seq_all date_seq_all

for jj0=[2 1 3:8]
    CLE_adj_grps=CLE_adj_delta_clone_grps{jj0};
    Dist_grps=Dist_delta_clone_grps{jj0};
    
    % the seq indices for the whole group
    ii1=ngrps_delta_lin{jj0,2};
    
    % the number of separate opt graphs for each different length
    % there is a row for each different length
    grp_opt=cellfun(@(x) length(x), CLE_adj_grps(:,2));
    
     % get the root sequence indices for each and the number of sequences
     % within that network
    iiy=find(grp_opt>0); % where there are some networks
    for i0=1:length(iiy) % the lengths where there is at least one optimal network
        i=iiy(i0);
        ii2=ngrps_delta_lin{jj0,3}{i}; % the indices within ii1
        Cclone=Dist_grps{i,4}; % the clone groups
        GG=CLE_adj_grps{i,4}; % the cell containing opt networks
        DD=CLE_adj_grps{i,2}; % the dates of the root nodes
        for j=1:grp_opt(i) % how many connected optimal graphs
            icount=icount+1;
            Gopt=GG{j};
            iroot=find(indegree(Gopt)==0 & outdegree(Gopt)>0); % the index within Gopt of the root clone
            iroot_name=Gopt.Nodes.Name(iroot);
            root_clone_index=str2num(iroot_name{1}); % the index within Cclone of the root node
            idx_seq_root_node=Cclone{root_clone_index}(1); % the first seq in this root Clone
            idx_seq_all(icount)=ii1(ii2(idx_seq_root_node)); % the sequence index within seqs_delta
            date_seq_all(icount)=DD(j);
            num_seqs_all(icount)=sum(Gopt.Nodes.Size);
        end
    end
end

 lengths_roots=zeros(length(idx_seq_all),13);
 for i=1:length(idx_seq_all)
     lengths_roots(i,:)=[(seqs_delta(idx_seq_all(i)).Domains(:,2)-seqs_delta(idx_seq_all(i)).Domains(:,1))' ...
         seqs_delta(idx_seq_all(i)).Domains(12,2)-seqs_delta(idx_seq_all(i)).Domains(1,1)];
 end
 seqs=seqs_delta(idx_seq_all);
 seqs(length(seqs)+1)=multi2(16); % the ref seq with Domains included

 multi2all=multialign(seqs,'UseParallel',true)       

%  multi2lim=multi2; % restrict sequence to common start (616)and end (29618)
%  for i=1:length(multi2)
%      multi2lim(i).Sequence=multi2(i).Sequence(616:29618);
%  end
% multi2alll was same as multi2all but Header changed to lineage and grp
% no.

 save('Save_outer_clone_network_all', "idx_seq_all",'multi2all','date_seq_all','num_seqs_all','multi2alll')