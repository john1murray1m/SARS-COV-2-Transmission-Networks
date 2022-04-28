%% determine Poisson and NBin fits to the networks for a Delta lineage
% extract just the largest group and analyse the mutational probability

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
% %  save('save_delta_pair_clone_Gopt','CLE_adj_delta_clone_grps','maxnt_dist','min_grp_size')
%   load('save_delta_pair_clone_Gopt')

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

% choose a lineage
jj0=8;
CLE_adj_grps=CLE_adj_delta_clone_grps{jj0};
% choose the row within this that picks out the group of same domain
% lengths
jrow=find(cellfun(@(x) ~isempty(x), CLE_adj_grps(:,1)));

% keep track of the figure counter
figcount=0;

clear pd_params
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
            
        
        %% determine the distribution of mutations
        nmut=[Gopt.Edges.Mut_Dist-1000; zeros(sum(Gopt.Nodes.Size)-height(Gopt.Nodes),1)];
        
        if mean(nmut)>var(nmut) % use Poisson
            pd_mut=fitdist(nmut,'Poisson');
            pd_params(figcount,1:2)=[pd_mut.ParameterValues, NaN];
            xx=0:6;
            y=pdf(pd_mut,xx);
            lambda=paramci(pd_mut);
            ydown = pdf('Poisson',xx,lambda(2));
            yup = pdf('Poisson',xx,lambda(1));
            
            figure(figcount+700)
            clf
            histogram(nmut)
            hold on
            plot(xx,y*length(nmut),'ro-','MarkerFaceColor','r')
            plot(xx,yup*length(nmut),'r--')
            plot(xx,ydown*length(nmut),'r--')
            xlabel('Mutations per transmission')
            set(gca,'YScale','log')
            ylabel('Number')
            title('Poisson')
        
        else
            % now fit with neg bin
            pd_mut_nbin=fitdist(nmut,'nbin');
            pd_params(figcount,1:2)=pd_mut_nbin.ParameterValues;
            y=pdf(pd_mut_nbin,xx);
            lambdan=paramci(pd_mut_nbin);
            ydown = nbinpdf(xx,lambdan(1,1),lambdan(1,2));
            yup = nbinpdf(xx,lambdan(2,1),lambdan(2,2));
            
            figure(figcount+700)
            clf
            histogram(nmut)
            hold on
            plot(xx,y*length(nmut),'ro-','MarkerFaceColor','r')
            plot(xx,yup*length(nmut),'r--')
            plot(xx,ydown*length(nmut),'r--')
            xlabel('Mutations per transmission')
            ylabel('Number')
            title('Neg Binomial')
            set(gca,'YScale','log')
        end
    end
end