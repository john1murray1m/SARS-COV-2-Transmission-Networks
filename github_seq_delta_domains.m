%% estimate the 12 ORF domains for each delta sequence through pairwise alignment with the ref sequence
% this will be improved (hopefully) in the next few steps. 
% Also correct some of the sequences that have missed the first A in the start codon
% through the sequencing procedure

% % save('Delta_seq_data_all','cov_delta_meta','no_seq_delta_all',...
% %     'seqs_delta','seq_ref','save_codes')
% load('Delta_seq_data_all')

% the reference sequence ORFs start and stop 
domain_refs=[266 13483; 13768 21555; 21563 25384; 25393 26220;...
            26245 26472; 26523 27191; 27202 27387; 27394 27759;...
            27756 27887; 27894 28259; 28274 29533; 29558 29674];
ndomain=size(domain_refs,1);

dr=domain_refs';
dom_align_calcs=cell(no_seq_delta_all,1); % the domain alignment endpoints for each sequence
parfor i=1:no_seq_delta_all
    [~,alignment]=nwalign(seqs_delta(i),seq_ref);
    ii1=alignment(1,:)~='-'; % where the seq alignment is not a gap
    ii2=alignment(3,:)~='-'; % where the ref_seq alignment is not a gap
    ii1sum=cumsum(ii1);
    ii2sum=cumsum(ii2);
    fred=zeros(2*ndomain,1);
    for j=1:(2*ndomain)
        ii3=find(ii2sum==dr(j));
        fred(j)=ii1sum(ii3(1));
    end
    if fred(1)==0 && strcmp(seqs_delta(i).Sequence(1:6),'TGGAGA')% some of the sequences are missing the initial 'A'
                % so assume the sequencing process missed this and insert
                % it at the start of the sequence and therefore move
                % everythin along by 1
                fred=fred+1;
                seqs_delta(i).Sequence=strcat('A',seqs_delta(i).Sequence);
    end
        
        fred=reshape(fred,2,ndomain);
    dom_align_calcs(i)={fred'};
end
 
save('Delta_seq_domains','dom_align_calcs','domain_refs',"ndomain",'seqs_delta')