%% determine the digraph within a specified clone and permute days to generate more likely graph structures


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
%  load('save_delta_pair_clone_Gopt_mutpen_10')

% specify markers and colour for mutations
mymarkers={'s','d','p','^','>','v','<'}; % 0, 1:6
cmap=colormap(hsv(7));

flag_out_permute=1; % if this is set to 1 the within nodes connecting to 
% outside clones are randomly chosen from nodes with the same day value

flag_in_permute=1; % if this is set to 1 then the date of submission is 
% replaced by a value modified by a lognormal distribution - then only use
% a single run


% choose one of the lineage groups in ngrps_delta_lin 1:8 (not 3)
jj0=2;
CLE_adj_grps=CLE_adj_delta_clone_grps{jj0};

jk=11; % choose one of the groups in that lineage
jjk=1; % choose the particular graph within this row

Gopt=CLE_adj_grps{jk,4}{jjk};

lin_name=CLE_adj_grps{jk,1};

Dist_grps=Dist_delta_clone_grps{jj0};

seqs00=seqs_delta(ngrps_delta_lin{jj0,2});
% now the subset for this domain length
iss=ngrps_delta_lin{jj0,3};
seqs0=seqs00(iss{jk});

Cclone=Dist_grps{jk,4};
clone_size=cellfun(@(x) length(x), Cclone);
% choose one of the clones
% ii1=find(clone_size==552); % for a particular sized clone or the 
% largest as calculated below. Cannot choose
% the root clone as the calculations below assume there is a preceding
% clone.
[~,ii1]=sort(clone_size,'descend');
iroot=find(indegree(Gopt)==0 & outdegree(Gopt)>0);
if iroot~=ii1(1) % the largest root is not the first one
    ii1=ii1(1);
else
    ii1=ii1(2);
end

Clone_node_name=string(ii1); % the name of the node for this clone
iij=find(strcmp(Gopt.Nodes.Name,Clone_node_name));
Date0_within=Gopt.Nodes.Day(iij); % the  date of the within seq root from the clone graph

seqs=seqs0(Cclone{ii1});
Dates=datenum([seqs.Date]);
[Date0,iid0]=min(Dates);

if flag_in_permute
    nperm=20;
    % lognormal
    logm=1;
    logs=1;
    
    logmu=log(logm^2/sqrt(logs+logm^2));
    logsig=sqrt(log(logs/logm^2+1));
    Dates00=repmat(Dates,nperm,1)-lognrnd(logmu,logsig,nperm,length(Dates)); % log normal dist sd 1 around the original dates
else
    nperm=1;
    Dates00=Dates;
end

dates_orig=Dates-min(Dates);
% collect the mean degree and table for these networks
tab_outdegree=cell(nperm,1);
for iperm=1:nperm
    dates0=Dates00(iperm,:);
    [dmin,iid0]=min(dates0);
    dates=dates0-dmin;
    ind_root=iid0;
    A=github_CLE_adjacency_within_clone(seqs,dates);
    
    % remove any edges coming into the root
    A(:,ind_root)=0;

    %% set up the initial digraph
    Nodes=string(Cclone{ii1}); % label the nodes to be their indices within the whole group
     NodeTable00=table(Nodes',dates','VariableNames',{'Name','Day'});
    NodeTable0=table(Nodes','VariableNames',{'Name'});
    
    G=digraph(A,NodeTable0);

    % number the edges
    G.Edges.Number=(1:numedges(G))';

    % the weights will be the values in A which are the prob
    % calculations modified by nt distance
    % modify the weights so they are unique
   G.Edges.Weight=G.Edges.Weight+(1:numedges(G))'/(numedges(G)*1e4); % so no edges have same weight

    [Gopt_within,errorflag]=CLE_optimal(G,ind_root);
    if errorflag==0 % CLE calculation successful
        Gopt_within.Nodes.Day=dates';
        opt_within={Gopt_within};
    end

    % get the Clone node of the clone graph coming in to the root
    fred=Gopt.Edges.EndNodes;
    ii2in=find(cellfun(@(x) strcmp(x,Clone_node_name), fred(:,2)));
    in_node=fred(ii2in,1); % the name of the node connecting into this clone at the root
    in_node_mut=Gopt.Edges.Mut_Dist(ii2in); % the mut dist between these nodes +1000
    in_node_idx=Gopt.Edges.Pre_Idx(ii2in); % the mut dist between these nodes +1000
    in_clone_date=datenum(seqs0(Cclone{str2num(in_node{1})}(in_node_idx)).Date)-Date0;
    in_clone_name=['Clone ',in_node{:}];
    % determine the name of the root node within this clone
    ii5=find(indegree(Gopt_within)==0 & outdegree(Gopt_within)>0);
    if length(ii5)>0
        [~,iia]=min(Gopt_within.Nodes.Day(ii5));
        ii5=ii5(iia); % choose the one with the lowest day value
    end
    root_within=Gopt_within.Nodes{ii5,1};

    % get the max of edge numbers
    in_max=max(Gopt_within.Edges.Number);
    
    NewNode=table({in_clone_name},in_clone_date,'VariableNames', {'Name' 'Day'});
    H=addnode(Gopt_within,NewNode);
    NewEdge=table([in_clone_name root_within],in_node_mut,in_max+1,...
        'VariableNames',{'EndNodes','Weight','Number'});
    H=addedge(H,NewEdge);

    % get the out edges going to the roots of other clones
    ii2out=find(cellfun(@(x) strcmp(x,Clone_node_name), fred(:,1)));
    out_node=fred(ii2out,2); % the name of the node connected to this clone 
    out_node_mut=Gopt.Edges.Mut_Dist(ii2out); % the mut dist between these nodes +1000
    out_node_idx=Gopt.Edges.Pre_Idx(ii2out); % the index within this clone of the startfor each edge
    out_node_idx_name=Gopt_within.Nodes{out_node_idx,1};
    % get the unique idx and distribute the number of times this is
    out_node_idx_nameu=out_node_idx_name; % a randomly drawn set of indices with the same day value
    if flag_out_permute % this is set to 1 if dont want all out edges originating from a single node
                        % for each connection from this clone
        [out_nameu,ia,ic]=unique(out_node_idx_name);
        for k=1:length(out_nameu)
            iia=(ic==k); % the number of times this index is used
            dayu=dates_orig(out_node_idx(ia(k))); % the day value
            iij=dates_orig==dayu; % the nodes with the same day
            iijj=randi(sum(iij),sum(iia),1); % indices drawn randomly amoung 
            % the same day indices for the number of times used in the
            % graph for out edges
            iij0=find(iij);
            out_node_idx_nameu(iia)=Gopt_within.Nodes.Name(iij0(iijj));
        end
    end
    out_clones_dates=zeros(length(out_node),1);
    out_clone_name=cell(length(out_node),1);
    for j=1:length(out_node)
        iij=find(strcmp(Gopt.Nodes.Name,out_node(j)));
        out_clones_dates(j)=Gopt.Nodes.Day(iij)-Date0_within;
        out_clone_name(j)={['Clone ',out_node{j}]};
    end
    
    NewNodes=table(out_clone_name,out_clones_dates,'VariableNames', {'Name' 'Day'});
    H=addnode(H,NewNodes);
    NewEdges=table([out_node_idx_nameu out_clone_name],out_node_mut,in_max+1+(1:length(ii2out))',...
        'VariableNames',{'EndNodes','Weight','Number'});
    H=addedge(H,NewEdges);
    
    iihin=find(H.Edges.Number==(in_max+1)); % the incoming edge
    iihout=find(H.Edges.Number>(in_max+1)); % the incoming edge
    figure(500+iperm)
    clf
    hplot=plot(H,'NodeLabel',[],'NodeColor','k','EdgeColor','k');
    highlight(hplot,'Edges',iihout,'EdgeColor','r')
    highlight(hplot,'Edges',iihin,'EdgeColor','g')

    highlight(hplot,in_clone_name,'Marker',mymarkers{in_node_mut-999},...
       'NodeColor',cmap(in_node_mut-998,:));
    for j=1:length(out_node)
       highlight(hplot,out_clone_name(j),'Marker',mymarkers(out_node_mut(j)-999),...
           'NodeColor',cmap(out_node_mut(j)-999,:));
    end
    layout(hplot,'force','UseGravity',true)
    title([lin_name,', grp = ',num2str([jj0 jk jjk]),', Clone ',char(Clone_node_name)])

    tab_outdegree(iperm)={tabulate(outdegree(Gopt_within))};
end

%% plot the histograms of the outdegrees
cmap1=colormap(hsv(nperm));
figure(300)
clf
for iperm=1:nperm
    sally=tab_outdegree{iperm};
    semilogy(sally(:,1),sally(:,2)/sum(sally(:,2)),'o-','Color',cmap1(iperm,:),...
        'MarkerFaceColor',cmap1(iperm,:));
    if iperm==1
        hold on
    end
end
pd=fitdist(outdegree(Gopt_within),'Poisson');
pd1=fitdist(outdegree(Gopt_within),'nbin');
xx=0:max(sally(:,1));
semilogy(xx,pdf(pd1,xx),'kd-','LineWidth',2)
    ylabel('P(Out degree=x)','Fontweight','bold')
    xlabel('x','Fontweight','bold')
    nci95=paramci(pd1);
% semilogy(xx,nbinpdf(xx,nci95(1,1),nci95(1,2)),'k--')
% semilogy(xx,nbinpdf(xx,nci95(2,1),nci95(2,2)),'k--')
set(gca,'YLim',[1e-3,1])

nmuts=sum(out_node_mut-1000);
mean_muts=nmuts/(numedges(Gopt_within)+length(out_node));
disp(['Mean number of mutations per transmission is ',num2str(mean_muts)])