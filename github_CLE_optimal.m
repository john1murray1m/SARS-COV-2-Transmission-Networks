%% run the CLE calculation on a given graph and root
function [Gopt,error_flag]=github_CLE_optimal(G,ind_root)

ncount=str2num(G.Nodes.Name{end}); % the last node number so can add new nodes with different names

% include a cell that tracks the edges kicked out by each
KickOut=cell(numedges(G),1);

% figure(2) % if want to plot the progression (slows calculations down
% though)
% clf
% plot(G,'EdgeLabel',G.Edges.Weight)

clear EdgeOpt
flag_opt=0;
G2=G;
iter=0;
while ~flag_opt
    iter=iter+1;
    [G1,G2,KickOut,flag_opt]=github_CLE_calc1(G2,KickOut,ncount);
    EdgeOpt(iter)={G1.Edges.Number};
%     figure(3) % if want to plot the progression (slows calculations down
%     though)
%     clf
%     subplot(1,2,1)
%     plot(G1,'EdgeLabel',G1.Edges.Number)
%     title('G1')
%      subplot(1,2,2)
%     plot(G2,'EdgeLabel',G2.Edges.Number)
%     title('G2')
   
    ncount=ncount+1;
end

    try % should work if ind_root is strongly connected to other nodes -otherwise will need to modify that choice
        j=0;
        OptEdges=EdgeOpt{iter-j};
        for j=0:(iter-2)
            % get the latest optedges
            KO=cell2mat(KickOut(OptEdges));
            prevopt=EdgeOpt{iter-j-1};
            prevopt=setdiff(prevopt,KO);
            OptEdges=union(OptEdges,prevopt);
        end
             
        EdgeTable=table([G.Edges.EndNodes(OptEdges,1),G.Edges.EndNodes(OptEdges,2)],G.Edges.Weight(OptEdges),G.Edges.Number(OptEdges),...
                'VariableNames',{'EndNodes','Weight','Number'});
        Gopt=digraph(EdgeTable,G.Nodes);
    
        error_flag=0;
    catch
        error_flag=1;
        Gopt=[];
    end
end