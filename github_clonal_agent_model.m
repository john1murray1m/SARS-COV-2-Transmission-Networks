%% generate clonal networks

% agents represented by an index for each individual, 
% a clone identifier, a time of infection, the index of the infecting agent,
% the clone of the infecting agent, the nt distance from the infecting
% agent

Ag=[1,1,0,0,0,0]; % start with one infected agent

clone_num=1; % the current largest clone number
% set the time period for each generation
dt=1;

params=[0.5 1.27 0.77 0.99]; % parameters for the two distributions

% specify the distribution of the number of infections for an individual
pd_inf=makedist('NegativeBinomial','p',params(1),'R',params(2));

% specify the distribution for the number of mutations for each
% transmission
pd_mut=makedist('NegativeBinomial','p',params(3),'R',params(4));

ngen=24; % number of generations
igen=1;
while igen <= ngen && size(Ag,1)<8000 % the network can grow quickly so bound the number of agents
    % for each agent established at the previous generation:
    iiold=find(Ag(:,3)==igen-dt);
    if ~isempty(iiold)
        % 1: determine the number of new agents established
        Ag_new=random(pd_inf,length(iiold),1); % a vector of no. new agents for each of the old
        num_new=sum(Ag_new);
        if num_new>0
            % 2: determine the number of mutations from its infecting agent
            Mut_new=random(pd_mut,num_new,1);
            Ag_new1=zeros(num_new,6);
            Ag_new1(:,1)=size(Ag,1)+(1:num_new)'; % the new indices
            Ag_new1(:,3)=igen; % the time of establishment
            Ag_new1(:,6)=Mut_new; % the mutations

            % get the agents that established non-zero new ones
            ii1=Ag_new>0;
            iiold1=iiold(ii1);
            Agg=Ag_new(ii1);
            icount=0;
            cc=Ag(iiold1,2); % The clone numbers of the infecting agents
            for inew=1:length(iiold1)
                Ag_new1(icount+(1:Agg(inew)),4)=iiold1(inew); % the infecting agent
                Ag_new1(icount+(1:Agg(inew)),2)=cc(inew);
                Ag_new1(icount+(1:Agg(inew)),5)=cc(inew);
                icount=icount+Agg(inew);
            end
            % now fix up the new clone numbers
            icc=Mut_new>0;
            if sum(icc)>0 % some new clone numbers need to be generated
                newcc=clone_num+(1:sum(icc))';
                Ag_new1(icc,2)=newcc;
                clone_num=clone_num+sum(icc); % the new max clone number
            end
            Ag=[Ag; Ag_new1];
        end
    end
    igen=igen+1;
end

[Aedges,ic,ia]=unique(Ag(:,[5 2]),'rows');
weight=Ag(ic,6);
% remove self connections
iww=weight>0;
Aedges=Aedges(iww,:);
weight=weight(iww);

% size of clones
tabs=tabulate(Ag(:,2));
clone_size=tabs(:,2);

Gmodel=digraph(Aedges(:,1),Aedges(:,2),weight);
Gmodel.Nodes.Size=clone_size;

figure(1)
clf
p=plot(Gmodel,'NodeLabel',[],'MarkerSize',sqrt(Gmodel.Nodes.Size));
layout(p,'force','UseGravity',true)
tt(1)={['Gens = ',num2str(ngen),', clones = ',num2str(clone_num),', agents = ',num2str(size(Ag,1))]};
tt(2)={['params = ',num2str(params)]};
title(tt)

% estimate clone size dist
sally=tabulate(Gmodel.Nodes.Size);
ii1=sally(:,2)>0;
sally=sally(ii1,:);

s1=sally(:,1);
s2=cumsum(sally(:,2),'reverse')/sum(sally(:,2));
% generate the linear best fit on a log-log basis
[bp,~,~,~,stats]=regress(log(s2(2:end)),[ones(size(s1(2:end))) log(s1(2:end))]);
alpha=-bp(2); %the power law exponent x^(-alpha)
R2=stats(1);
ss=(1:max(s1));
figure(2)
clf
loglog(s1,s2,'bo','MarkerFaceColor','b');
ylabel('P(Clone size>=x)','Fontweight','bold')
xlabel('x','Fontweight','bold')

hold on
loglog(s1,exp(bp(1)+bp(2)*log(s1)),'b--')
text(max([2,0.1*max(s1)]), 0.5,['\alpha = ',num2str(-bp(2),2)])
