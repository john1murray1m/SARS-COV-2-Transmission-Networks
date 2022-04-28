%% plot the optimal networks, mutational paths, and total numbers over time for each lineage

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
% 
% %  save('save_delta_pair_clone_Gopt_mutpen_10','CLE_adj_delta_clone_grps','maxnt_dist','min_grp_size')
%  load('save_delta_pair_clone_Gopt_mutpen_10')

% % save('Save_delta_lengths','T_delta_lengths')
% load('Save_delta_lengths')

T_delta_lengths=readtable('Delta_lengths.xlsx'); % read in the lengths of
% the orfs for each group

domain_name(1)={'1a'};
domain_name(2)={'1b'};
domain_name(3)={'S'};
domain_name(4)={'3a'};
domain_name(5)={'E'};
domain_name(6)={'M'};
domain_name(7)={'6'};
domain_name(8)={'7a'};
domain_name(9)={'7b'};
domain_name(10)={'8'};
domain_name(11)={'N'};
domain_name(12)={'10'}; 

%specify the start and end of the 1a shifts for each j0. Refers to the rows
%in T_delta_lengths for AANoToAddFor1aMut. Nothing for j0=3
AAshift_ind=[16 17; 1 15; 0 0; 18 18; 19 20; 21 24; 25 30; 31 31];

% choose the lineage group from 1:8 (nothing for 3)
jj0=2;

CLE_adj_grps=CLE_adj_delta_clone_grps{jj0};

% get the vector of AA modifications to 1a for this jj0
oneamodvec=T_delta_lengths{AAshift_ind(jj0,1):AAshift_ind(jj0,2) ,'AANoToAddFor1aMut'};

% choose the row within this that picks out the group of same domain
% lengths
jrow=find(cellfun(@(x) ~isempty(x), CLE_adj_grps(:,1)));

% keep track of the figure counter
figcount=0;
max_clone_size=[];

for jk0=1:length(jrow)
    jk=jrow(jk0);
    
    % choose the particular connected component within this
    for jjk=1:length(CLE_adj_grps{jk,2})
    
        % update the figure counter
        figcount=figcount+1;
        
        Gopt=CLE_adj_grps{jk,4}{jjk};
        lin_name=CLE_adj_grps{jk,1};
        
        Dist_grps=Dist_delta_clone_grps{jj0};
        
        Date_val=datestr(CLE_adj_grps{jk,2}(jjk));
        figure(figcount)
        clf
        popt=plot(Gopt,'NodeLabel',Gopt.Nodes.Day,'MarkerSize',sqrt(Gopt.Nodes.Size));
        title([lin_name,', grp = ',num2str([jj0 jk jjk]),', no. clones = ',num2str(height(Gopt.Nodes)),...
            ', no. seqs = ',num2str(sum(Gopt.Nodes.Size)),', ',Date_val])
        if height(Gopt.Nodes)>50
                     layout(popt,'force','UseGravity',true)
        else
                     layout(popt,'layered')
        end        
        
        % generate a power law plot of clone size
        sally=tabulate(Gopt.Nodes.Size);
        ii1=sally(:,2)>0;
        sally=sally(ii1,:);
        s1=sally(:,1);
        s2=cumsum(sally(:,2),'reverse')/sum(sally(:,2));
        % generate the linear best fit on a log-log basis
        [bp,~,~,~,stats]=regress(log(s2(2:end)),[ones(size(s1(2:end))) log(s1(2:end))]);
        alpha(figcount)=-bp(2); %the power law exponent x^(-alpha)
        alphacnum(figcount)=height(Gopt.Nodes); % the number of clones
        R2(figcount)=stats(1);
        reg_length(figcount)=length(s1);
        ss=(1:max(s1));
        figure(100+figcount)
        clf
        loglog(s1,s2,'bo','MarkerFaceColor','b');
        ylabel('P(Clone size>=x)','Fontweight','bold')
        xlabel('x','Fontweight','bold')
        
        hold on
        loglog(s1,exp(bp(1)+bp(2)*log(s1)),'b--')
        text(max([2,0.1*max(s1)]), 0.5,['\alpha = ',num2str(-bp(2),2)])
        title([lin_name,', grp = ',num2str([jj0 jk jjk]),', no. clones = ',num2str(height(Gopt.Nodes)),...
            ', no. seqs = ',num2str(sum(Gopt.Nodes.Size)),', ',Date_val])
        
        max_clone_size(figcount,1:3)=[max(Gopt.Nodes.Size), sum(Gopt.Nodes.Size), max(Gopt.Nodes.Size)/sum(Gopt.Nodes.Size)*100];
        
        %% determine the longest path so can produce a plot of mutations
        ind_root=find(indegree(Gopt)==0 & outdegree(Gopt)>0); % the root node (ignore singletons)
        
        % get the leaf nodes - the ones with outdegree 0
        ileaf=find(outdegree(Gopt)==0 & indegree(Gopt)>0);
        
        % all paths will be from the root to these nodes
        nleafs=length(ileaf);
        pathi=cell(nleafs,1); % the paths
        for i=1:nleafs
            pathi(i)=allpaths(Gopt,ind_root,ileaf(i));
        end
        
        length_pathi=cellfun(@(x) length(x), pathi); % their lengths
        % the longest path
        [milong,ilong]=max(length_pathi);
        
        figure(figcount+200)
        clf
         popt1=plot(Gopt,'NodeLabel',[],'MarkerSize',sqrt(Gopt.Nodes.Size));
        title([lin_name,', grp = ',num2str([jj0 jk jjk]),', no. clones = ',num2str(height(Gopt.Nodes)),...
            ', no. seqs = ',num2str(sum(Gopt.Nodes.Size)),', ',Date_val])
        if height(Gopt.Nodes)>50
                     layout(popt1,'force','UseGravity',true)
        else
                     layout(popt1,'layered')
        end
        highlight(popt1,pathi{ilong},'EdgeColor','r','LineWidth',1.5,'NodeColor','r')
        highlight(popt1,pathi{ilong},'EdgeColor','r','LineWidth',1.5,'NodeColor','r')
        highlight(popt1,ind_root,'NodeColor','g','MarkerSize',sqrt(Gopt.Nodes.Size(ind_root)))
        
        %% plot the mut pathway
        seqs0=seqs_delta(ngrps_delta_lin{jj0,2});
        % now the subset for this domain length
        iss=ngrps_delta_lin{jj0,3};
        seqs=seqs0(iss{jk});
            
        Cclone=Dist_grps{jk,4};
        oneamod=oneamodvec(figcount);
        path_diff=github_plot_delta_mutpath1(seqs,Cclone,pathi{ilong},Gopt.Nodes,domain_name,300+figcount,oneamod);
        
        %% plot the dynamics with state info
        
        github_cov_state_plot(seqs, 400+figcount)
        
        %% determine the probability of mutation per transmission
        % the number of mutations over the whole network
        num_mut(figcount)=sum(Gopt.Edges.Mut_Dist-1000);
        num_clone_edges(figcount)=height(Gopt.Edges);
        num_within_edges(figcount)=sum(Gopt.Nodes.Size-1);
        prob_mut(figcount)=num_mut(figcount)/(num_clone_edges(figcount)+num_within_edges(figcount));
    end
end

AAA=[num_mut', num_clone_edges', num_within_edges', prob_mut', alpha' R2', reg_length'];
AAA=AAA(1:figcount,:)

max_clone_size