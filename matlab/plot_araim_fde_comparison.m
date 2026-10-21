function plot_araim_fde_comparison(export_figures)
% Paired FDE-toggle experiment; default keeps visible windows for copying.
if nargin<1, export_figures=false; end
result_dir='F:/RTKLib-LAB/tmp/robust_araim_fde_on_20261015';
out_dir=fullfile(result_dir,'figures');
T=readtable(fullfile(result_dir,'comparison.csv'),'VariableNamingRule','preserve');
sites={'JC01','JC02','JC04'}; methods={'igg3','ransac'};
fields={'availability_pct','H_H_V_coverage_pct','H_H_V_RMS_m'};
labels={'PL availability (%)','Horizontal coverage (%)','Horizontal RMS (m)'};
fig=figure('Visible','on','Color','w','Position',[60 60 1550 950]);
layout=tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for s=1:3
    for metric=1:3
        ax=nexttile; hold(ax,'on'); vals=zeros(2,2);
        for m=1:2
            for enabled=0:1
                ix=strcmp(T.site,sites{s}) & strcmp(T.method,methods{m}) & T.FDE==enabled;
                assert(nnz(ix)==1);
                vals(m,enabled+1)=T.(fields{metric})(ix);
            end
        end
        b=bar(ax,vals,'grouped');
        b(1).FaceColor=[.65 .70 .76]; b(2).FaceColor=[.10 .32 .69];
        pure=strcmp(T.site,sites{s}) & strcmp(T.method,'pure_fde1');
        hp=yline(ax,T.(fields{metric})(pure),'--','Color',[.65 .17 .13],'LineWidth',1.3);
        xticks(ax,1:2); xticklabels(ax,{'IGGIII','RANSAC'});
        if metric<3, ylim(ax,[0 105]); else, ylim(ax,[0 max([vals(:);T.(fields{metric})(pure)])*1.25]); end
        set(ax,'FontSize',12); grid(ax,'on'); box(ax,'on');
        if s==1, title(ax,labels{metric},'FontSize',13); end
        if metric==1, ylabel(ax,sites{s},'FontSize',14); end
        for m=1:2
            for enabled=1:2
                if metric<3, textLabel=sprintf('%.2f',vals(m,enabled)); else, textLabel=sprintf('%.4f',vals(m,enabled)); end
                text(ax,b(enabled).XEndPoints(m),vals(m,enabled),textLabel,...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
            end
        end
        if s==1 && metric==1
            lg=legend(ax,[b(1),b(2),hp],{'ARAIM FDE off','ARAIM FDE on','Pure ARAIM (FDE on)'},'Orientation','horizontal');
            lg.Layout.Tile='north';
        end
    end
end
if export_figures
    if ~isfolder(out_dir), mkdir(out_dir); end
    stem=fullfile(out_dir,'FDE_on_off_comparison');
    exportgraphics(fig,[stem '.png'],'Resolution',300);
    exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
    savefig(fig,[stem '.fig']);
end
drawnow;
end
