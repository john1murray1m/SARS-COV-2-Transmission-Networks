%% determine better estimates of the domains using the initial domain estimates from the pair alignments
% 1/2/22 This is for delta

% % save('Delta_seq_data_all','cov_delta_meta','no_seq_delta_all',...
% %     'seqs_delta','seq_ref','save_codes')
% load('Delta_seq_data_all')
% 
% % save('Delta_seq_domains','dom_align_calcs','domain_refs',"ndomain",'seqs_delta')
% load('Delta_seq_domains')

    [dom_nums,seq_group_probs]=github_seq_delta_domain_orfs(seqs_delta,dom_align_calcs); 
    for i=1:no_seq_delta_all
        seqs_delta(i).Domains=dom_nums(:,:,i);
    end


    save('Delta_seqs_domain_orfs',"seqs_delta",'seq_group_probs')