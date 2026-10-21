function plot_current_pl_pe_ksigma(export_figures)
if nargin < 1, export_figures = false; end
% Current nine-run paper figures. Change paths only; values are saved results.
result_dir = 'F:/RTKLib-LAB/tmp/gls_pfa4_pmd3_20261014';
sigma_dir = fullfile(result_dir, 'figures_sigma_contrast_log');
output_dir = 'F:/RTKLib-LAB/matlab/figures/current_pl_pe_ksigma';
if ~isfolder(output_dir), mkdir(output_dir); end
sites = {'JC01','JC02','JC04'};
methods = {'pure_fde1','igg3','ransac'};
names = {'ARAIM','IGGIII + ARAIM','RANSAC + ARAIM'};
data = cell(3,3);
limits = zeros(3,2,2);
for s=1:3
    for m=1:3
        T = readtable(fullfile(result_dir,[sites{s} '_' methods{m} '_epochs.csv']));
        S = readtable(fullfile(sigma_dir,[sites{s} '_' methods{m} '_sigma.csv']));
        t = T.time_gpst; ts = S.time_gpst;
        if ~isdatetime(t), t=datetime(t,'InputFormat','yyyy-MM-dd HH:mm:ss'); end
        if ~isdatetime(ts), ts=datetime(ts,'InputFormat','yyyy-MM-dd HH:mm:ss'); end
        [found, idx] = ismember(t,ts);
        kh=nan(height(T),1); kv=kh;
        kh(found)=S.Kmd_sigma_H_m(idx(found)); kv(found)=S.Kmd_sigma_V_m(idx(found));
        data{s,m} = struct('t',hours(t-dateshift(t,'start','day')), ...
            'H',[T.HPE_m,T.HPL_m,kh], 'V',[T.VPE_m,T.VPL_m,kv], 'q',T.Q);
    end
    for d=1:2
        field = 'H'; if d==2, field='V'; end
        values=[];
        for m=1:3, a=data{s,m}.(field); values=[values;a(:)]; end %#ok<AGROW>
        values=values(isfinite(values) & values>0);
        limits(s,d,:)=[min(values)*.8,max(values)*1.3];
    end
end
for m=1:3
    fig=figure('Visible','on','Color','w','Position',[80 80 1600 800]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for d=1:2
        field='H'; direction='Horizontal (m)';
        if d==2, field='V'; direction='Vertical (m)'; end
        for s=1:3
            nexttile; hold on; a=data{s,m}.(field); t=data{s,m}.t;
            hpl=scatter(t,a(:,2),5,[23 70 162]/255,'filled');
            hks=scatter(t,a(:,3),5,[0 158 56]/255,'o','filled');
            hpe=scatter(t,a(:,1),5,[228 26 28]/255,'filled');
            set(gca,'YScale','log','FontSize',11);
            ylim(squeeze(limits(s,d,:))'); xlim([0 24]); xticks(0:6:24); grid on; box on;
            bad=isfinite(data{s,m}.q) & ~isfinite(a(:,2)); yl=ylim;
            hu=scatter(t(bad),repmat(yl(2)/1.03,nnz(bad),1),14,[.47 .47 .47],'|');
            ex=isfinite(a(:,2)) & a(:,1)>a(:,2);
            scatter(t(ex),a(ex,1),15,[.42 .09 .14],'x','HandleVisibility','off');
            if s==1, ylabel(direction); end
            if d==1, title(sites{s}); else, xlabel('GPST (h)'); end
            if s==1 && d==1
                legend([hpe,hpl,hks,hu],{'PE','PL','K_{MD}\sigma','Unavailable'},'Location','best','FontSize',10);
            end
        end
    end
    sgtitle(names{m},'FontSize',14);
    stem=fullfile(output_dir,[methods{m} '_PL_PE_Ksigma_log']);
    if export_figures
        exportgraphics(fig,[stem '.png'],'Resolution',300);
        exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
        savefig(fig,[stem '.fig']);
    end
    drawnow;
end
fprintf('PL/PE/Ksigma figure windows are open; export enabled: %d.\n',export_figures);
end
