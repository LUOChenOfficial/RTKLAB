function figures = plot_robustSS_paper(export_figures, show_windows)
%PLOT_ROBUSTSS_PAPER Draw fixed/float envelopes from the verified FDE runs.
%   plot_robustSS_paper             opens six editable figure windows.
%   plot_robustSS_paper(true)       also exports PNG, vector PDF and FIG.
%   plot_robustSS_paper(true,false) exports without showing windows.
%   Fig. 5: three fixed-epoch PL/PE/Ksigma figures, one per site.
%   Fig. 6: three float-epoch PL/PE/Ksigma figures, one per site.
%   Saved CSVs already contain posterior-scaled Ksigma; do not scale again.
%   All available epochs are used. Scatter marks never connect missing epochs.
close all;
if nargin<1, export_figures=false; end
if nargin<2, show_windows=true; end

% Paths: change these when selecting a different prepared result dataset.
script_dir=fileparts(mfilename('fullpath'));
data_dir=fullfile(script_dir,'figures','verified_paper','source_data');
output_dir=fullfile(script_dir,'figures','robustSS_paper');

sites={'JC01','JC02','JC04'};
methods={'pure_fde1','igg3','ransac'};
names={'ARAIM','IGGIII + ARAIM','RANSAC + ARAIM'};
blue=[23 70 162]/255;green=[0 150 65]/255;red=[228 26 28]/255;
visibility='off';if show_windows,visibility='on';end
if export_figures && ~isfolder(output_dir),mkdir(output_dir);end
figures=gobjects(6,1);
data=cell(3,3);limits=zeros(3,2,2);
for s=1:3
    for m=1:3
        path=fullfile(data_dir,[sites{s} '_' methods{m} '.csv']);
        assert(isfile(path),'Missing prepared data: %s',path);
        T=readtable(path);
        t=T.time_gpst;
        if ~isdatetime(t),t=datetime(t,'InputFormat','yyyy-MM-dd HH:mm:ss');end
        data{s,m}=struct('t',hours(t-dateshift(t,'start','day')),'q',T.Q,...
            'H',[T.HPE_m T.HPL_m T.Ksigma_H_m],...
            'V',[T.VPE_m T.VPL_m T.Ksigma_V_m]);
    end
    for d=1:2
        field='H';if d==2,field='V';end
        values=cell(3,1);
        for m=1:3,a=data{s,m}.(field);values{m}=a(:);end
        values=vertcat(values{:});values=values(isfinite(values)&values>0);
        assert(~isempty(values),'No positive %s data for %s.',field,sites{s});
        limits(s,d,:)=[min(values)*.7 max(values)*1.4];
    end
end

% Figures 5/6: same per-site/direction limits across both states and methods.
index=0;
for q=1:2
    state='Fixed';number=5;
    if q==2,state='Float';number=6;end
    for s=1:3
        index=index+1;
        figures(index)=new_figure(sprintf('Fig. %d | %s | %s',number,sites{s},state));
        for d=1:2
            field='H';if d==2,field='V';end
            for m=1:3
                ax=nexttile;hold(ax,'on');a=data{s,m}.(field);t=data{s,m}.t;
                selected=data{s,m}.q==q;
                good=@(v) selected&isfinite(v)&v>0;
                hp=scatter(ax,t(good(a(:,2))),a(good(a(:,2)),2),8,blue,'filled');
                hk=scatter(ax,t(good(a(:,3))),a(good(a(:,3)),3),9,green,'o','filled');
                he=scatter(ax,t(good(a(:,1))),a(good(a(:,1)),1),8,red,'filled');
                set(ax,'YScale','log');ylim(ax,squeeze(limits(s,d,:))');
                format_time(ax);
                if d==1,title(ax,names{m},'FontWeight','normal');end
                if m==1,ylabel(ax,[field ' (m)']);end
                if d==2,xlabel(ax,'GPST (h)');end
                if d==1&&m==1
                    lg=legend(ax,[he hp hk],{'PE','PL','K_{MD}\sigma'},'Orientation','horizontal');
                    lg.Layout.Tile='north';
                end
            end
        end
        sgtitle([sites{s} '  |  ' state],'FontName','Arial','FontSize',15,'FontWeight','normal');
        deliver(figures(index),sprintf('Fig%d_%s_%s_envelope',number,sites{s},lower(state)));
    end
end
if show_windows,figure(figures(1));end
fprintf('Fixed/float envelopes: %d windows; export=%d; visible=%d.\n',numel(figures),export_figures,show_windows);

    function fig=new_figure(name)
        fig=figure('Visible',visibility,'Color','w','Position',[60 60 1440 760],...
            'Name',name,'NumberTitle','off');
        tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    end
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
