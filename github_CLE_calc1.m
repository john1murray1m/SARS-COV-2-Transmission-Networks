function [G1,G2,KickOut,flag_opt]=github_CLE_calc1(G20,KickOut,ncount)
% perform one iteration of the Chu-Liu-Edmonds algorithm
% Input:
% G20 is the output from the previous iteration of the modified graph with
% one contracted node
% KickOut is a cell which lists the edges kicked out by each other edge,
% cell(nedges,1)
% ncount is the number of nodes in the graph to date

% Output:
% G1 is the graph produced at the hopeful optimal step. It contains the
% edges for each iteration that will be used in the expansion stage
% G2 is the updated graph,
% KickOut is the updated cells for edges kicked out by others
% flag_opt returns 1 if optimal and stops further iterations

tOut=G20.Edges.EndNodes(:,2);
sOut=G20.Edges.EndNodes(:,1);
% determine the incoming edge with min weight and
% form the subgraph - incoming edges have been previously removed from root
nnodes=unique(tOut);
nedges=zeros(size(nnodes)); % the edge index for each terminal node
nweights=zeros(size(nnodes)); % the weights of the lowest incoming edge
edgenum=zeros(size(nnodes)); % the original edge numbers
wts=G20.Edges.Weight;
nums=G20.Edges.Number;
for i=1:length(nnodes)
    ii1=find(strcmp(G20.Edges.EndNodes(:,2),nnodes(i)));
    [nweights(i),ii2]=min(wts(ii1));
    nedges(i)=ii1(ii2);
    edgenum(i)=nums(ii1(ii2));
end

EdgeTable=table([sOut(nedges),tOut(nedges)],nweights,edgenum,...
        'VariableNames',{'EndNodes','Weight','Number'});
NodeTable=table(G20.Nodes.Name,'VariableNames',{'Name'});
G1=digraph(EdgeTable,NodeTable);
% figure(20) % if want to plot progression 
% clf
% subplot(1,2,1)
% plot(G1,'EdgeLabel',G1.Edges.Number)
% subplot(1,2,2)
% plot(G1,'EdgeLabel',G1.Edges.Weight)

if hascycles(G1)
    % this is not optimal
    [cycles,edgecycles] = allcycles(G1);
%     highlight(p,'Edges',edgecycles{1},'EdgeColor','r','LineWidth',1.5,'NodeColor','r','MarkerSize',6)


    ncount=ncount+1; % add a new node
    ii1=find(~ismember(G20.Nodes.Name,cycles{1})); % the nodes not in the cycle
    G2=subgraph(G20,G20.Nodes.Name(ii1));    
    G2=addnode(G2,string(ncount));
    
    % adjust the edges from the cycle nodes (the new composite node)
    % to the rest
    end_nodes=G20.Edges.EndNodes(:,2); % edge end nodes 
    start_nodes=G20.Edges.EndNodes(:,1);  % edge start nodes
    ind_outedges=find(ismember(start_nodes,cycles{1}) & ~ismember(end_nodes,cycles{1}) );
    % sort out which of these are different
    if ~isempty(ind_outedges)
        end_nodes1=unique(end_nodes(ind_outedges));
        len_out=length(end_nodes1);
        wt_out=NaN(len_out,1); % the weight of the new edge going out from the cycle
        edge_out=NaN(len_out,1); % the edge number of the modified edge going out from the cycle
        for j=1:len_out
            ii1=find(strcmp(end_nodes(ind_outedges),end_nodes1(j)));
            [wt_out(j),iimin]=min(G20.Edges.Weight(ind_outedges(ii1)));
            edge_out(j)=G20.Edges.Number(ind_outedges(ii1(iimin)));
        end
        EdgeTable=table([string(ncount*ones(len_out,1)),end_nodes1],wt_out,edge_out,...
                'VariableNames',{'EndNodes','Weight','Number'});
        G2=addedge(G2,EdgeTable);
%         figure(30)
%         clf
%         subplot(1,2,1)
%          plot(G2,'EdgeLabel',G2.Edges.Number);
%         subplot(1,2,2)
%          plot(G2,'EdgeLabel',G2.Edges.Weight);
    end
      
    % now the edges coming into the cycle 
    ind_inedges=find(ismember(end_nodes,cycles{1}) & ~ismember(start_nodes,cycles{1}) );
    if ~isempty(ind_inedges)>0
        cy_endnodes=cycles{1}; % the nodes in the cycle
        len_in=length(cy_endnodes);
        edge_cyend=edgecycles{1}; % these edges in the cycle are listed according to start node
        edge_cyend=[edge_cyend(end), edge_cyend(1:(end-1))]; % now edges are listed according to end node
        cy_wts=G1.Edges.Weight(edge_cyend);
        cy_num=G1.Edges.Number(edge_cyend);
        for j=1:len_in
            ii1=find(strcmp(end_nodes(ind_inedges),cy_endnodes(j)));
            if ~isempty(ii1)
                wt1=cy_wts(j); % its weight
                ikick0=cy_num(j); % its number that will be kicked out if the incoming edge is used
                iienter=ind_inedges(ii1); % the edges entering
                wt_in=G20.Edges.Weight(iienter)-wt1;
                num_in=G20.Edges.Number(iienter);
                %for k=1:length(ii1)
                for k=1:length(iienter)
                    KickOut(num_in(k))={[KickOut{num_in(k)}; ikick0]}; % the edge this one will kick out added to previous list
                end
                EdgeTable=table([start_nodes(iienter),string(ncount*ones(length(ii1),1))],wt_in,...
                        num_in, 'VariableNames',{'EndNodes','Weight','Number'});
                G2=addedge(G2,EdgeTable);
            end
        end
    end
%         figure(40)
%         clf
%         subplot(1,2,1)
%         plot(G2,'EdgeLabel',G2.Edges.Number)
%         subplot(1,2,2)
%         plot(G2,'EdgeLabel',G2.Edges.Weight)

    flag_opt=0;

else
    G2=G1;
    flag_opt=1;

end

end
