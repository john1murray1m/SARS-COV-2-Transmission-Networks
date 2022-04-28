%% plot the mutations at all leaf nodes compared to the root node for a group in the clone network


function [all_non,leaf_mut]=github_plot_delta_mut_leaf(seqs0,Cclone,ind_root,ind_leaf,nodes,domain_name,fignum)

% place all indices into path with the root node first
path=[ind_root; ind_leaf];

lin_name=seqs0(1).Lineage;

% nodes will be the digraph node names such as Gopt.Nodes that has a field
% called 'Name' that will correspond to the original clone number rather
% than the number index of path.
% the original clone numbers
iclone=str2num(char(nodes{path,'Name'}));

% the first sequences from each of these clones
ind_seq=cellfun(@(x) x(1,1),Cclone(iclone));

n_clone=cellfun(@(x) length(x),Cclone(iclone));

seqs1=seqs0(ind_seq);

domain_refs=seqs1(1).Domains;

% check to see if any of these are NaN and if so whether any other seqs
% besides the root have a nonNaN
ii1=find(isnan(arrayfun(@(x) sum(x.Domains(:)), seqs1))); % the clones with a problem
if ~isempty(ii1)
    seqs1=seqs1(~ii1); % drop off any where they don't have all domains defined
end

nseqs=length(seqs1);

% get the first sequences of the clones within the path
ndomain=size(seqs1(1).Domains,1);

    domcol=[0 0.9 0.9; 1 1 0; 1 0 0; 0.5 0.5 1;  0.5 1 0.5; 1 0.5 0; 0 1 0.7; 0 0.7 1; 1 0 0.7;...
    0.7 1 0; 1 0.5 0.5; 0 0.5 0.5];

% the positions and frequencies of N
all_non=cell(ndomain,1);

% the positions and frequency of mutations excluding N
leaf_mut=cell(ndomain,1);

for k=1:ndomain
    harry=arrayfun(@(x) nt2int(x.Sequence(x.Domains(k,1):x.Domains(k,2))), seqs1,'UniformOutput',false);
    harry=cell2mat(harry'); % replace with a matrix
    fred=(harry>4 & harry<16); % the mostly N positions
    iirootn=fred(1,:)==1;
    all_non(k)={[1:size(fred,2); sum(fred)]}; % the first row is the position and the second the frequency

    % mutations compared to the first row (root node)
    george=harry-repmat(harry(1,:),size(harry,1),1);

    % don't want to count N so replace thes with 0
    ii1=fred(:)==1;
    george(ii1)=0;

    % also don't want to count where the root node is N
    george(:,iirootn)=0;
    george1=george~=0; % where there are mutation

    leaf_mut(k)={[1:size(george1,2); sum(george1)]};
end

figure(fignum)
clf
ht1=80;
ht2=10;
for k=1:ndomain
    fred=all_non{k};
    plot(domain_refs(k,1)-1+fred(1,:),fred(2,:)/nseqs*100)
    if k==1
        hold on
    end
    pos=domain_refs(k,:);
    pos1=pos(1);
    wid=pos(2)-pos(1);
    col=domcol(k,:);
    rectangle('Position',[pos1,ht1,wid,ht2],'FaceColor',col,'EdgeColor',col); 
    if ismember(k, [5 7 9])
        text(sum(pos)/2,10+ht1+0.75*ht2,domain_name{k},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontWeight','bold','FontSize',11)
    else
        text(sum(pos)/2,10+ht1+ht2/4,domain_name{k},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontWeight','bold','FontSize',11)
    end
end
xlabel('Nt position')
ylabel('nonACTG%')
title(lin_name)
set(gca,'YLim',[0 100])

figure(fignum+1)
clf
for k=1:ndomain
    fred=leaf_mut{k};
    plot(domain_refs(k,1)-1+fred(1,:),fred(2,:)/nseqs*100)
    if k==1
        hold on
    end
    pos=domain_refs(k,:);
    pos1=pos(1);
    wid=pos(2)-pos(1);
    col=domcol(k,:);
    rectangle('Position',[pos1,ht1,wid,ht2],'FaceColor',col,'EdgeColor',col); 
    if ismember(k, [5 7 9])
        text(sum(pos)/2,10+ht1+0.75*ht2,domain_name{k},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontWeight','bold','FontSize',11)
    else
        text(sum(pos)/2,10+ht1+ht2/4,domain_name{k},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontWeight','bold','FontSize',11)
    end
end
xlabel('Nt position')
ylabel('Mutation %')
title(lin_name)
set(gca,'YLim',[0 100])
