%% determine the mutation differences in a domain for clone seqs
function path_diffk=github_mut_seqs_dom(seqs1,kdom,flag_aant)

% path_diffk is the cell of path differences for these sequences in domain
% kdom
% for each domain determine the nt that differ from the root excluding
% nonACTG
path_diffk=cell(length(seqs1)-1,5); % only the first 3 are determined if flag_aant==0;

% input:
% seqs1 are the main sequences along the path,
% kdom is the domain number
% flag_aant is 1 if can calculate the AA values otherwise just the first 3
% components of the output

if flag_aant
    fred=arrayfun(@(x) x.Domains(kdom,:),seqs1,'UniformOutput',false);
else
    fred=arrayfun(@(x) [x.Domains(kdom-1,2), x.Domains(kdom+1,1)],seqs1,'UniformOutput',false);
end
dom_bounds=cell2mat(fred');

iirootd=dom_bounds(1,:);
fred1=seqs1(1).Sequence(iirootd(1):iirootd(2));
numnon1=(nt2int(fred1)>4 & nt2int(fred1)<16);
iiroot=numnon1==0; % the positions in the kth domain for the root sequence that are ACTG
for i=1:(length(seqs1)-1)
    iirootsd=dom_bounds(1+i,:);
    if ~isnan(sum(iirootsd))
        fred2=seqs1(1+i).Sequence(iirootsd(1):iirootsd(2));
        numnons=(nt2int(fred2)>4 & nt2int(fred2)<16);
        iiroots=numnons==0; % the positions in the kth domain for the root sequence that are ACTG
        iicommon=iiroot&iiroots; % the common ACTG posns for both sequences
        iidiff0=(fred1~=fred2);
        iidiff=find(iicommon & iidiff0); % the numbers within the domain where nt differ and are not nonACTG
        if ~isempty(iidiff)
            path_diffk(i,1)={iidiff};
            path_diffk(i,2)={fred1(iidiff)}; % the root nt values
            path_diffk(i,3)={fred2(iidiff)}; % the i+1'st seq nt values
            if flag_aant
                ii0=floor((iidiff-1)/3)+1;
                path_diffk(i,4)={ii0}; % the AA posns
                harry=strings(2,length(ii0));
                for j=1:length(ii0)
                    ii1=3*(ii0(j)-1)+(1:3); % the codon nt posns
                    sally1=fred1(ii1);
                    sally2=fred2(ii1);
    
                    % original AA
                    if nt2int(sally1)<=4 
                        aaorig=nt2aa(sally1);
                    else
                        aaorig='X';
                    end
                    % mutated sequence codon
                    if nt2int(sally2)<=4 
                        aanew=nt2aa(sally2);
                    else
                        aanew='X';
                    end
                    harry(1,j)=aaorig;
                    harry(2,j)=aanew;
                end
                path_diffk(i,5)={harry}; % the collection of AA for the original and i+1'st sequence (in rows)
            end

        end
    end
end
