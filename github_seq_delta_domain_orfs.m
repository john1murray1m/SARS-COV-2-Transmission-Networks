%% input a set of sequences, the estimates of the ORF start and stop from the previous pairwise alignment calculations.
% Output the best estimates for the domains and where there are problems


function [dom_nums,seq_group_probs]=github_seq_delta_domain_orfs(seqs,dom_align_calcs)


startCodons = {'ATG'}; % Defaults to standard genetic code
stopCodons = {'TAA','TAG','TGA'}; % Defaults to standard genetic code

ndomain=12; % the number of domains 

% iframes=zeros(ndomain,7); % the orfs (possibly -not used), start range, length range, likely length and start of each gene
iframes=[2 1 700 10000 16000 13218 266;... % 1a 
    1 12000 15000 4000 9000 8088 13768; ...% 1b
     2 17000 24000 3000 5000 3822 21563; ... % S 
     1 24000 27000 600 1200 828 25393; ...% 3a
     1 25000 27000 130 350 228 26245; ...% E
     3 25500 27500 400 900 669 26523; ...% M
     1 26500 27700 100 280 186 27202;... % 6
     1 26500 28000 200 500 366 27394; ... % 7a
     3 26400 28500 90 200 132 27756; ... % 7b
     3 26500 28500 300 450 366 27894; ... % 8 
     2 27000 28500 750 1400 1260 28274; ...% N
     2 28000 30000 90 300 117 29558]; % 10
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

% usual starting AA for these ORF
protAA={'MESLVP', 'MVPHIS', 'MFVFLV', 'MDLFMR','MYSFVS', 'MADSNG', 'MFHLVD', ...
    'MKIILF', 'MIELSL', 'MKFLVF', 'MSDNGP', 'MGYINV'};
protnt=aa2nt(protAA); % if the AA version causes problems when there are Ns then use nt values

protAAint=NaN(length(protAA),6);
protntint=NaN(length(protnt),6*3);
for i=1: length(protAA)
    protAAint(i,:)=aa2int(protAA{i});
    protntint(i,:)=nt2int(protnt{i});
end

    dom_nums=NaN(ndomain,2,length(seqs)); % the possible start and stop nt posns for each sequence
    dom_nums_check=zeros(ndomain,length(seqs)); % this will be set to 1 if same as aligned ORF or matches the protein start
    dom_nums_check_sparse=sparse(ndomain,length(seqs)); % now this will be 1-dom_nums_check
                                                    % a fraction relative
                                                    % to the amount of
                                                    % match otherwise


    for k=1:length(seqs)
        seq0=seqs(k).Sequence;
        dom_nums0=dom_align_calcs{k}; % possible domain start and end 
        % will check that it matches with start and stop. Start is the
        % start of the orf but Stop is the start of the Stop codon. So for
        % the end of the domain, need to add 2 to Stop
        domain_vals=NaN(ndomain,2); % possible domain start and stop if don't match the ref seq
        domain_vals_match=zeros(ndomain,1); % measures match if the ORFS don't align with the ref seq

        orfs=seqshoworfs(seq0,'nodisplay',1,'MinimumLength',30);
        % fix if there is no stop for a start by using the end of the seq
        % -2
        for ii1=1:3
                if length(orfs(ii1).Start)~=length(orfs(ii1).Stop)
                    orfs(ii1).Stop=[orfs(ii1).Stop length(seqs)-2];
                end
        end
        starts=[orfs.Start];
        stops=[orfs.Stop];

        for kk1=1:ndomain
            if kk1~=3 % the start of S needs to be outside ORF1ab
                ii3=find(starts==dom_nums0(kk1,1));
                if ~isempty(ii3)
                    st1=starts(ii3);
                    st2=stops(ii3);
                    if st2+2==dom_nums0(kk1,2)
                        dom_nums(kk1,:,k)=[st1 st2+2];
                        dom_nums_check(kk1,k)=1;
                    elseif st2-st1 >= iframes(kk1,4) && st2-st1<= iframes(kk1,5) % stop not the same but length OK
                        dom_nums(kk1,:,k)=[st1 st2+2];
                        dom_nums_check(kk1,k)=1;
                    end
                end

            % the above works if the start numbering match with ref sequence and the
            % stops numbering  match or  give a length within the right bounds.
            % Otherwise need to search for protein matches at the start
            if dom_nums_check(kk1,k)==0

                if kk1>1 && any(dom_nums_check(1:(kk1-1),k)==1)
                    [dd1,imm]=nanmax(dom_nums(1:(kk1-1),2,k));
                    if imm==kk1-1
                        ii2=find(starts>= dd1-100 & starts<= iframes(kk1,3) & ...
                            stops-starts >= iframes(kk1,4) & stops-starts <= iframes(kk1,5) ...
                            & starts<= dd1+100);
                    else
                         ii2=find(starts>= dd1-100 & starts<= iframes(kk1,3) & ...
                            stops-starts >= iframes(kk1,4) & stops-starts <= iframes(kk1,5));
                    end

                else

                        ii2=find(starts>= iframes(kk1,2) & starts<= iframes(kk1,3) & ...
                            stops-starts >= iframes(kk1,4) & stops-starts <= iframes(kk1,5));
                end
                if ~isempty(ii2)
                    dvall=starts(ii2)';
                                             
                    dvmatch=zeros(size(dvall,1),1); % the number of AA matches for each
                    % if this is the only one then accept
                    if size(dvall,1)==1
                          dom_nums(kk1,:,k)=[starts(ii2), stops(ii2)+2];
                          dom_nums_check(kk1,k)=0.9999;
                    else
                      
                        for jj=1:size(dvall,1)
                            if dvall(jj,1)+17<length(seq0)
                                ss=seq0(dvall(jj,1):(dvall(jj,1)+17));
                               if ~any(nt2int(ss)>4)
                                    dvmatch(jj)=sum(protAAint(kk1,:)==aa2int(nt2aa(ss)))/6;
                               else
                                   dvmatch(jj)=sum(protntint(kk1,:)==nt2int(ss))/18;
                               end
                            end
                        end
                        % get the closest in AAmatch
                        [fred,ii7]=max(dvmatch);
                        start_best=dvall(ii7(1),1);
                        % now get the stop associated with this
                        st2=stops(ii2(ii7(1)));
                        if fred>0.5
                            dom_nums(kk1,:,k)=[start_best, st2+2];
                        end
                        dom_nums_check(kk1,k)=fred;
                    end
                end

            end
            % now repeat this for S since its start is not in orf since it
            % is within ORF1ab
            else % kk1==3
                startsr = sort(cell2mat(regexpi(seq0,startCodons)));
                stopsr = sort(cell2mat(regexpi(seq0,stopCodons)));
                ii3=find(startsr==dom_nums0(kk1,1));
                if ~isempty(ii3)
                    st1=startsr(ii3);
                    ii4=find(stopsr>st1 & mod(stopsr-st1,3)==0 ); % need to be in same reading frame
                    if ~isempty(ii4) 
                        st2=stopsr(ii4(1));
                        if st2-st1>=iframes(kk1,4) && st2-st1<=iframes(kk1,5)
                            dom_nums(kk1,:,k)=[st1 st2+2];
                            dom_nums_check(kk1,k)=1;
                        end
                    end
                end
                if dom_nums_check(kk1,k)==0
                    if dom_nums_check(kk1-1,k)>0
                      ii2=find(startsr>= iframes(kk1,2) & startsr<= iframes(kk1,3) & ...
                          startsr>dom_nums(kk1-1,2,k));
                    else
                       ii2=find(startsr>= iframes(kk1,2) & startsr<= iframes(kk1,3));
                       
                    end
                   if ~isempty(ii2)
                        dvall=[];
                        for ii3=1:length(ii2)
                            st=startsr(ii2(ii3)); % current start being tested
                            ii4=find(stopsr>st  & mod(stopsr-st,3)==0); % the next stop in the same frame
                            % want first stop 
                            % length reasonable
                            if  stopsr(ii4(1))-st>= iframes(kk1,4) && stopsr(ii4(1))-st<= iframes(kk1,5)
                                dvall=[dvall;st stopsr(ii4(1))];
                            end
                        end
                        if ~isempty(dvall)
                            % now get the first start with each stop
                            [~,ia,~]=unique(dvall(:,2));
                            dvall=dvall(ia,:);
                            if size(dvall,1)==1 % accept this
                                dom_nums(kk1,:,k)=[dvall(1,1), dvall(1,2)+2];
                                dom_nums_check(kk1,k)=0.9999;
                            else            
                                dvmatch=zeros(size(dvall,1),1); % the number of AA matches for each
                                for jj=1:size(dvall,1)
                                    if dvall(jj,1)+17<length(seq0)
                                        ss=seq0(dvall(jj,1):(dvall(jj,1)+17));
                                       if ~any(nt2int(ss)>4)
                                            dvmatch(jj)=sum(protAAint(kk1,:)==aa2int(nt2aa(ss)))/6;
                                       else
                                           dvmatch(jj)=sum(protntint(kk1,:)==nt2int(ss))/18;
                                       end
                                    end
                                 end
                                 % get the closest in AAmatch
                                 [fred,ii7]=max(dvmatch);
                                 start_best=dvall(ii7(1),1);
                                 st2=dvall(ii7(1),2);
                                 if fred>0.5
                                    dom_nums(kk1,:,k)=[start_best, st2+2];
                                 end
                                 dom_nums_check(kk1,k)=fred;
                            end
                        end
                    end
                end
            end
        end
    end

    ii1=dom_nums_check(:)<1;
    dom_nums_check_sparse(ii1)=1-dom_nums_check(ii1);
    seq_group_probs={dom_nums_check_sparse};
 end

