%% plot the mutations along a pathway of a group of sequences for omicron
% add in oneamod to be able to modify the 1a posns since they were
% dependent on where the sequencing started. This is the shift to be added
% to the aa posn for 1a muts.

function path_diff=github_plot_delta_mutpath1(seqs0,Cclone,path,nodes,domain_name,fignum,oneamod)

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
domain_refs(1,1)=266; % this is the seq_ref start of 1a

% check to see if any of these are NaN and if so whether any other seqs
% besides the root have a nonNaN
ii1=find(isnan(arrayfun(@(x) sum(x.Domains(:)), seqs1))); % the clones with a problem
if isempty(ii1)
    iflag_dom_prob=0;
    dom_probs_id=[];
else
    iflag_dom_prob=1;    
    dom_probs=cell(length(ii1),1); % the domains that are problems for these sequences
    for i=1:length(ii1)
        ii2=find(isnan(seqs1(ii1(i)).Domains(:,1)));
        dom_probs(i)={ii2};
    end
    dom_probs_id=unique(cell2mat(dom_probs)); % the domains that have issues. for Omicron it is 3 or 11 or both
end


% get the first sequences of the clones within the path
ndomain=size(seqs1(1).Domains,1);

    domcol=[0 0.9 0.9; 1 1 0; 1 0 0; 0.5 0.5 1;  0.5 1 0.5; 1 0.5 0; 0 1 0.7; 0 0.7 1; 1 0 0.7;...
    0.7 1 0; 1 0.5 0.5; 0 0.5 0.5];


%% mutations along the sequence of sequences

% for each domain determine the nt that differ from the root excluding
% nonACTG
path_diff=cell(length(seqs1)-1,ndomain,5);

for k=1:ndomain
    if iflag_dom_prob==0 || (iflag_dom_prob==1 && ~ismember(k,dom_probs_id))
        flag_aant=1;
        path_diff(:,k,:)=github_mut_seqs_dom(seqs1,k,flag_aant);

    else 
        iirootd=seqs1(1).Domains(k,:); % some seqs have problems but at least first OK
        if ~isnan(sum(iirootd))
            flag_aant=1;
            path_diff(:,k,:)=github_mut_seqs_dom(seqs1,k,flag_aant);
        elseif k>1 && k<ndomain
            dlims=arrayfun(@(x) x.Domains(k+1,1)-x.Domains(k-1,2), seqs1);
            iij=~isnan(dlims);
           
            if ~isempty(iij)
                dlims=dlims(iij);
                dlims=unique(dlims); % if there is only one value can use this to define the nt position
                                    % from the 
                if length(dlims)==1 % one of the problem domains without proper 
                    % definition of domain boundaries so use outside ends
                    domain_refs(k,:)=[seqs1(1).Domains(k-1,2), seqs1(1).Domains(k+1,1)]; % replace the NaNs for 
                    
                    if ~isnan(sum(domain_refs(k,:)))
                        % this domain for plotting and reference
                        flag_aant=0;
                        path_diff(:,k,:)=github_mut_seqs_dom(seqs1,k,flag_aant);
                    end
                end
            end
        end
    end
end
                

figure(fignum)
clf
nseqs=length(seqs1);
for j=1:nseqs
    box_ht=0.85/nseqs;
    box_diff=0.9/nseqs;
    yht=0.05+(nseqs-j)*box_diff;
    axes('Position',[0.05,yht,0.9,box_ht],'Visible','off')
    rectangle('Position',[0,-0.5,domain_refs(ndomain,2),2.5],'FaceColor','w','EdgeColor','w')
    hold on
    % shift the text for E, 6, 10 up (5, 7 9)
    for k=1:ndomain
            pos=domain_refs(k,:);
            if ~isnan(sum(pos))
                pos1=pos(1);
                wid=pos(2)-pos(1);
                col=domcol(k,:);
                rectangle('Position',[pos1,-0.5,wid,2.5],'FaceColor',col,'EdgeColor',col); 
                if j==1 % only include the gene names for the first non-mut plot
                    xj=pos;
                    if ismember(k, [5 7 9])
                        text(sum(xj)/2,2.5+0.8,domain_name{k},'HorizontalAlignment','center',...
                            'FontWeight','bold','FontSize',11)
                    else
                        text(sum(xj)/2,2.5,domain_name{k},'HorizontalAlignment','center',...
                            'FontWeight','bold','FontSize',11)
                    end
                end
            end
    end
    for k=1:ndomain
            pos=domain_refs(k,:);
            if ~isnan(sum(pos))
                pos1=pos(1);
                if j>1 
                    if k==1
                        xposns=path_diff{j-1,k,1}+pos1+3*oneamod; %modify the posn dependent on when the sequencing started. 
                                       % The 1a start is now the seq_ref
                                       % value and here use the
                                       % multialignment of the root
                                       % sequences.
                                                                   
                    else
                        xposns=path_diff{j-1,k,1}+pos1;
                    end

                    if ~isempty(path_diff{j-1,k,4}) % the AA were calculated
                        harry=path_diff{j-1,k,5};
                        if k==1
                            aanum=path_diff{j-1,k,4}+oneamod;
                        else  
                            aanum=path_diff{j-1,k,4};
                        end
                    else
                        harry1=path_diff{j-1,k,2};
                        harry2=path_diff{j-1,k,3};
                        harry=strings(2,length(harry1));
                        for jj=1:length(harry1)
                            harry(1,jj)=harry1(jj);
                            harry(2,jj)=harry2(jj);
                        end
                        if k==1
                            aanum=path_diff{j-1,k,1}+3*oneamod; % the nt numbers
                        else
                            aanum=path_diff{j-1,k,1}; % the nt numbers
                        end
                    end


                    for jk=1:length(xposns)
                        line([xposns(jk) xposns(jk)], [-0.5 2],'Color','k','LineWidth',1)
                        if harry(1,jk)==harry(2,jk)
                            ytxt=strcat(harry(1,jk),num2str(aanum(jk)));
                        else
                            ytxt=strcat(harry(1,jk),num2str(aanum(jk)),harry(2,jk));
                        end
                        text(xposns(jk),0.75,ytxt,'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom','FontWeight','bold','FontSize',6,'Rotation',270)
                          
                    end
                end

            end
    end
    text(domain_refs(ndomain,2)+200,0.75,num2str(n_clone(j)),'FontWeight','bold')
end

end

   

