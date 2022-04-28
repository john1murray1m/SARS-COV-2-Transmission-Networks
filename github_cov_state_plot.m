%% plot a group of sequences over time separated into their states
function github_cov_state_plot(seqs,fignum)

smarkers={'s','o','d','^','v','>','<','p','*'};
scols=colormap(hsv(9));

States={'New South Wales', 'Victoria', 'Queensland', 'Western Australia', 'Tasmania',...
        'Northern Territory', 'Australian Capital Territory', 'Unknown'};

figure(fignum)
clf
tab_d=tabulate([seqs.Date]);
tab_nums=cell2mat(tab_d(:,2));
max_tab=max(tab_nums);
y_state_gap=max_tab/30; % the gap between state markers
    plot(datetime(tab_d(:,1)),tab_nums,'ko-')
    hold on
    % now plot markers for each state at the times they occur
    ileg=1;
    ikleg=[];
    for k=1:length(smarkers)
        if k<length(smarkers)
            ii1=find(strcmp([seqs.Division],States{k}));
            if ~isempty(ii1)
                tfred=datetime([seqs(ii1).Date]);
                h(ileg)=plot(tfred,-k*y_state_gap*ones(size(tfred)),'LineStyle','none','Marker',smarkers{k},...
                    'Color',scols(k,:),'MarkerFaceColor',scols(k,:),'DisplayName',States{k});
                ileg=ileg+1;
            end
        else
            ii1=find(arrayfun(@(x) isempty(x), [seqs.Division]));
            if ~isempty(ii1)
                tfred=datetime([seqs(ii1).Date]);
                h(ileg)=plot(tfred,-k*y_state_gap*ones(size(tfred)),'LineStyle','none','Marker',smarkers{k},...
                    'Color',scols(k,:),'MarkerFaceColor',scols(k,:),'DisplayName',States{k});
            end
        end
            
    end
    legend(h,'Location','northwest')
end
