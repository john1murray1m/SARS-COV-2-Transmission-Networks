%% load in the Australian Delta sequence and metadata up to
% 31/10/21 date of submission

save_folder='GISAID_SEQ/DELTA'; % folder for the downloaded data

save_codes={'1643523228981', '1643523710527', '1643526481436', '1643527562911'}; % downloaded delta sequence filenames

icount=0;
for ifiles=1:length(save_codes)
    % choose the sequence fasta file and meta file
    fastafile=['./',save_folder,'/',save_codes{ifiles},'.sequences.fasta'];
    metafile=['./',save_folder,'/',save_codes{ifiles},'.metadata.txt'];


    % sequences
    [cov_delta_head,cov_delta_seq]=fastaread(fastafile);
    
    % metadata
    opts=detectImportOptions(metafile,'FileType','text','Delimiter','\t');
    cov_delta_meta0=readtable(metafile,opts);
    
    % change age to cell of strings for the last file since some are
    % numbers and others are text
    if ifiles==4
        cov_delta_meta0.age=arrayfun(@num2str,cov_delta_meta0.age, 'UniformOutput',0);
    end

    no_seqs_delta=length(cov_delta_head);
    for i=1:no_seqs_delta
        seqs_delta(icount+i).Header=cov_delta_head{i};
        seqs_delta(icount+i).Sequence=cov_delta_seq{i};
        % get the correct line in the meta file - this can be out of alignment
        imeta=find(strcmp(cov_delta_meta0{:,'strain'},cov_delta_head{i}));
        seqs_delta(icount+i).Strain=cov_delta_meta0{imeta,'strain'};
        seqs_delta(icount+i).Date=cov_delta_meta0{imeta,'date'};
        seqs_delta(icount+i).Sex=cov_delta_meta0{imeta,'sex'};
        seqs_delta(icount+i).Age=cov_delta_meta0{imeta,'age'};
        seqs_delta(icount+i).Lineage=cov_delta_meta0{imeta,'pangolin_lineage'};
        seqs_delta(icount+i).Clade=cov_delta_meta0{imeta,'GISAID_clade'};
        seqs_delta(icount+i).Division=cov_delta_meta0{imeta,'division'};
        seqs_delta(icount+i).Location=cov_delta_meta0{imeta,'location'};
        seqs_delta(icount+i).Laboratory=cov_delta_meta0{imeta,'originating_lab'};
        
    end
    icount=icount+no_seqs_delta;

    % gather the meta data
    if ifiles==1
        cov_delta_meta=cov_delta_meta0;
    else
        cov_delta_meta=[cov_delta_meta; cov_delta_meta0];
    end
end
no_seq_delta_all=icount; % the total number of sequences

clear cov_delta_meta0 cov_delta_head cov_delta_seq

%% load in the reference sequence

% read in ref sequence from Wang_JMV_2020
[ref_head,ref_seq]=fastaread('EPI_ISL_412026.fasta');
seq_ref.Header=ref_head;
seq_ref.Sequence=ref_seq; 

save('Delta_seq_data_all','cov_delta_meta','no_seq_delta_all',...
    'seqs_delta','seq_ref','save_codes')