function plot_verified_paper(export_figures, show_windows)
% Current verified FDE results only. Default: visible, editable figure windows.
% Run prepare_verified_paper_data.py once after changing the experiment paths.
if nargin<1, export_figures=false; end
if nargin<2, show_windows=true; end
source_dir='F:/RTKLib-LAB/matlab/figures/verified_paper/source_data';
output_dir='F:/RTKLib-LAB/matlab/figures/verified_paper';
if ~isfolder(output_dir), mkdir(output_dir); end
visibility='off'; if show_windows, visibility='on'; end
sites={'JC01','JC02','JC04'};
methods={'pure_fde1','igg3','ransac'};
names={'ARAIM','IGGIII + ARAIM','RANSAC + ARAIM'};
short={'ARAIM','IGGIII','RANSAC'};
blue=[23 70 162]/255; green=[0 150 65]/255; red=[228 26 28]/255;
method_colors=[blue;green;red];
metrics=readtable(fullfile(source_dir,'metrics.csv'),'TextType','string');
D=cell(3,3); limits=zeros(3,2,2);
for s=1:3
    for m=1:3
        T=readtable(fullfile(source_dir,[sites{s} '_' methods{m} '.csv']));
        t=T.time_gpst;
        if ~isdatetime(t), t=datetime(t,'InputFormat','yyyy-MM-dd HH:mm:ss'); end
        D{s,m}=struct('t',hours(t-dateshift(t,'start','day')),'q',T.Q,...
            'H',[T.HPE_m T.HPL_m T.Ksigma_H_m],...
            'V',[T.VPE_m T.VPL_m T.Ksigma_V_m]);
    end
    for d=1:2
        field='H';if d==2, field='V';end
        values=[];
        for m=1:3, a=D{s,m}.(field);values=[values;a(:)];end %#ok<AGROW>
        values=values(isfinite(values)&values>0);
        limits(s,d,:)=[min(values)*.7 max(values)*1.4];
    end
end

% Main temporal evidence: no connecting lines across missing/state-excluded epochs.
for s=1:3
    for q=1:2
        state='Fixed';if q==2, state='Float';end
        fig=figure('Visible',visibility,'Color','w','Position',[60 60 1440 760]);
        tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        for d=1:2
            field='H';if d==2,field='V';end
            for m=1:3
                ax=nexttile;hold(ax,'on');a=D{s,m}.(field);t=D{s,m}.t;selected=D{s,m}.q==q;
                good=@(v) selected&isfinite(v)&v>0;
                hp=scatter(ax,t(good(a(:,2))),a(good(a(:,2)),2),8,blue,'filled');
                hk=scatter(ax,t(good(a(:,3))),a(good(a(:,3)),3),9,green,'o','filled');
                he=scatter(ax,t(good(a(:,1))),a(good(a(:,1)),1),8,red,'filled');
                set(ax,'YScale','log');ylim(ax,squeeze(limits(s,d,:))');
                missing=selected&~isfinite(a(:,2));yl=ylim(ax);
                hu=scatter(ax,t(missing),repmat(yl(2)/1.08,nnz(missing),1),30,[.4 .4 .4],'|');
                format_time(ax);
                if d==1,title(ax,names{m},'FontWeight','normal');end
                if m==1,ylabel(ax,[field ' (m)']);end
                if d==2,xlabel(ax,'GPST (h)');end
                if d==1&&m==1
                    lg=legend(ax,[he hp hk hu],{'PE','PL','K_{MD}\sigma','PL unavailable'},'Orientation','horizontal');
                    lg.Layout.Tile='north';
                end
            end
        end
        sgtitle([sites{s} '  |  ' state],'FontName','Arial','FontSize',15,'FontWeight','normal');
        deliver(fig,[sites{s} '_' lower(state) '_envelope']);
    end
end

% Equal-epoch denominator for PL versus Ksigma; stratification is also exported.
for category={'all','fixed','float'}
    group=category{1};
    fig=figure('Visible',visibility,'Color','w','Position',[80 80 1320 700]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for d=1:2
        direction='H';if d==2,direction='V';end
        for s=1:3
            ax=nexttile;vals=zeros(3,2);
            for m=1:3
                row=metrics.site==sites{s}&metrics.method==methods{m}&metrics.group==group&metrics.direction==direction;
                assert(nnz(row)==1);
                vals(m,:)=[metrics.PL_coverage_pct(row),metrics.Ksigma_coverage_pct(row)];
            end
            b=bar(ax,vals,'grouped');b(1).FaceColor=blue;b(2).FaceColor=green;
            hold(ax,'on');yline(ax,100,':','Color',[.55 .55 .55],'HandleVisibility','off');
            ylim(ax,[0 110]);yticks(ax,0:25:100);xticks(ax,1:3);xticklabels(ax,short);
            set(ax,'FontName','Arial','FontSize',12);box(ax,'on');grid(ax,'on');ax.GridAlpha=.12;
            if d==1,title(ax,sites{s},'FontWeight','normal');end
            if s==1,ylabel(ax,[direction ' coverage (%)']);end
            for k=1:2
                for m=1:3
                    text(ax,b(k).XEndPoints(m),vals(m,k)+1.4,sprintf('%.1f',vals(m,k)),...
                        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'FontName','Arial');
                end
            end
            if d==1&&s==1
                lg=legend(ax,b,{'PL','K_{MD}\sigma'},'Orientation','horizontal');lg.Layout.Tile='north';
            end
        end
    end
    sgtitle(['Coverage  |  ' group],'FontName','Arial','FontSize',15,'FontWeight','normal');
    deliver(fig,['coverage_' group]);
end

fig=figure('Visible',visibility,'Color','w','Position',[100 100 1320 410]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for s=1:3
    ax=nexttile;vals=zeros(3,2);
    for m=1:3
        row=metrics.site==sites{s}&metrics.method==methods{m}&metrics.group=='all'&metrics.direction=='H';
        vals(m,:)=[metrics.availability_pct(row),metrics.fixed_pct(row)];
    end
    b=bar(ax,vals,'grouped');b(1).FaceColor=blue;b(2).FaceColor=[.68 .71 .75];
    ylim(ax,[0 110]);yticks(ax,0:25:100);xticks(ax,1:3);xticklabels(ax,short);
    set(ax,'FontName','Arial','FontSize',12);grid(ax,'on');box(ax,'on');ax.GridAlpha=.12;
    title(ax,sites{s},'FontWeight','normal');if s==1,ylabel(ax,'Rate (%)');end
    for k=1:2
        for m=1:3
            text(ax,b(k).XEndPoints(m),vals(m,k)+1.4,sprintf('%.1f',vals(m,k)),...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontName','Arial','FontSize',10);
        end
    end
    if s==1
        lg=legend(ax,b,{'PL availability','Fixed epochs'},'Orientation','horizontal');lg.Layout.Tile='north';
    end
end
deliver(fig,'availability_fixed_rate');

fig=figure('Visible',visibility,'Color','w','Position',[100 100 1320 700]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for d=1:2
    field='H';if d==2,field='V';end
    for s=1:3
        ax=nexttile;hold(ax,'on');h=gobjects(3,1);
        for m=1:3
            a=D{s,m}.(field);valid=isfinite(a(:,1))&a(:,1)>0;
            h(m)=scatter(ax,D{s,m}.t(valid),a(valid,1),7,method_colors(m,:),'filled');
        end
        set(ax,'YScale','log');format_time(ax);
        if d==1,title(ax,sites{s},'FontWeight','normal');end
        if s==1,ylabel(ax,[field ' PE (m)']);end
        if d==2,xlabel(ax,'GPST (h)');end
        if d==1&&s==1
            lg=legend(ax,h,names,'Orientation','horizontal');lg.Layout.Tile='north';
        end
    end
end
deliver(fig,'position_error');
fprintf('Verified paper figures complete. Export: %d; visible windows: %d.\n',export_figures,show_windows);

    function format_time(ax)
        xlim(ax,[0 24]);xticks(ax,0:6:24);grid(ax,'on');box(ax,'on');
        set(ax,'FontName','Arial','FontSize',12);ax.GridAlpha=.12;
    end
    function deliver(fig,stem)
        drawnow;
        if export_figures
            exportgraphics(fig,fullfile(output_dir,[stem '.png']),'Resolution',180);
            exportgraphics(fig,fullfile(output_dir,[stem '.pdf']),'ContentType','vector');
            savefig(fig,fullfile(output_dir,[stem '.fig']));
        end
    end
end
