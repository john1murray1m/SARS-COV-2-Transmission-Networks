%% calculate adjacency connections between delta clones within each group
 
% % save('Delta_seq_domains','dom_align_calcs','domain_refs',"ndomain",'seqs_delta')
% load('Delta_seq_domains')
% 
% %     save('Delta_seqs_domain_orfs',"seqs_delta",'seq_group_probs')
%     load('Delta_seqs_domain_orfs')
% 
% % save('Delta_lineage_domain_groups',"ngrps_delta_lin",'delta_lineages','grp_min')
% load('Delta_lineage_domain_groups')
% 
% %    save('seq_delta_pair_clone_distance','Dist_delta_clone_grps')
%    load('seq_delta_pair_clone_distance')


% set the max nt distance between connected sequences - same as ntdiffmax
% in nonACTG_calculate
maxnt_dist=6; 

% set the minimum group size - this is  for clones but leave it at 10 to
% be consistent with omicron
min_grp_size=10;

CLE_adj_delta_clone_grps=cell(size(ngrps_delta_lin,1),1);
% the lineage, the indices rel to seqs_new, the cell of subgroup 
% indices within the larger index set,
% the cell of optimal graphs for the
% subgroups, the date representing day 0

for jk=1:size(ngrps_delta_lin,1)
    lineage=ngrps_delta_lin{jk,1};
    seqs0=seqs_delta(ngrps_delta_lin{jk,2}); 
    % now get any subgroups specified by domain lengths
    sally=ngrps_delta_lin{jk,3};
    Dist_grps=Dist_delta_clone_grps{jk};
    CLE_adj_grps=cell(length(sally),4);
    for jjk=1:length(sally)
        seqs=seqs0(sally{jjk});
        Dates=datenum([seqs.Date]);
    
        % first calculate possible date probabilities of group j transmitting
        % to group k using the prob
        Cclone=Dist_grps{jjk,4};
    
        nclones=length(Cclone);
        dp=sparse(nclones,nclones); % the prob clone j transmitted to clone k in dp(j,k). 
                                                   % the node in k is always
                                                   % the root.
        dpj=sparse(nclones,nclones); % the index within the Cclone{j} sequences that gives the max prob
    
        [dp,dpj]=github_prob_edge(Cclone,Dates); % calculate the probability and the index within the j clone that trhere is
                                            % a transmission from j to k, based
                                            % on date differences.
    
         %% now determine the nt difference between the groups that have a possible connection
        dn=github_nt_delta_edge1(seqs,Cclone,dp);
    
        % make the node size relative to the number of sequences
        node_size=cellfun(@(x) length(x),Cclone);
        node_date=cellfun(@(x) min(Dates(x)), Cclone); % the earliest date for each clone
    
         [A0,bingrps]=github_CLE_adjacency_clone(dp,dn,maxnt_dist,min_grp_size);
        
         if sum(bingrps>0) % there will be some grps of suff size
            iib0=bingrps>0;
            bin_nums=unique(bingrps(iib0));
            ind_grp_grps=cell(1,length(bin_nums)); % the indices rel to ii11
            opt_grp_grps=cell(1,length(bin_nums)); % the optimal digraph for this group
            node_grp_grps=cell(1,length(bin_nums)); % the NodeTable for this group
            datenum0_grp_grps=zeros(1,length(bin_nums)); % the starting datenum for this group
    
            for i=1:length(bin_nums)
                ii13=find(bingrps==bin_nums(i));
                ind_grp_grps(i)={ii13};
                [Date0,iid0]=min(node_date(ii13));
                datenum0_grp_grps(i)=Date0;
                A=A0(:,ii13);
                A=A(ii13,:);
                % also if necessary remove any self-loops
                A=A-diag(diag(A));
                % specify the root of the tree - this does not always work if
                % mutations between the earliest sequences don't allow a strong
                % connection from this root
                ind_root=iid0;
    
                % remove any edges coming into the root
                A(:,ind_root)=0;
    
                dday=node_date(ii13)-Date0; % the days relative to earliest
    
                %% set up the initial digraph
                Nodes=string(ii13); % label the nodes to be their indices within the whole group
                 NodeTable00=table(Nodes',dday',node_size(ii13)','VariableNames',{'Name','Day','Size'});
                NodeTable0=table(Nodes','VariableNames',{'Name'});
                
                G=digraph(A,NodeTable0);
    
                % number the edges
                G.Edges.Number=(1:numedges(G))';
    
                % the weights will be the values in A which are the prob
                % calculations modified by nt distance
                % modify the weights so they are unique (otherwise the
                % arborescence calculations can have problems)
                G.Edges.Weight=G.Edges.Weight+(1:numedges(G))'/(numedges(G)*1e4); % so no edges have same weight
    
                [Gopt,errorflag]=github_CLE_optimal(G,ind_root);
    
                % get the endNodes and convert them to indices 
                fred=Gopt.Edges.EndNodes;
                harry=cellfun(@(x) str2num(x), fred);
                pre_idx=zeros(size(harry,1),1); 
                mut_dist=zeros(size(harry,1),1); 
                for jj=1:length(pre_idx)
                    pre_idx(jj)=dpj(harry(jj,1),harry(jj,2));
                    mut_dist(jj)=dn(harry(jj,1),harry(jj,2));
                end
                Gopt.Edges.Pre_Idx=pre_idx;
                Gopt.Edges.Mut_Dist=mut_dist;
    
    %             % add back in the Node attributes
                Gopt.Nodes.Day=dday';
                Gopt.Nodes.Size=node_size(ii13)';
    %             figure(1)
    %             clf
    %             plot(Gopt,'MarkerSize',Gopt.Nodes.Size)
    
                if errorflag==0 % CLE calculation successful
                    opt_grp_grps(i)={Gopt};
                end
            end
            CLE_adj_grps(jjk,1)={lineage};
            CLE_adj_grps(jjk,3)={ind_grp_grps};
            CLE_adj_grps(jjk,4)={opt_grp_grps};
            CLE_adj_grps(jjk,2)={datenum0_grp_grps};
         end
    end
    CLE_adj_delta_clone_grps(jk)={CLE_adj_grps};
end

save('save_delta_pair_clone_Gopt_mutpen_10','CLE_adj_delta_clone_grps','maxnt_dist','min_grp_size')