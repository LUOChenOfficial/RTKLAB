function plot_satellite_overview(export_figures)
if nargin < 1, export_figures = false; end
% Display saved original satellite data; break trajectories across time gaps.
input_dir='F:/RTKLib-LAB/matlab/figures/satellites';
output_dir=input_dir;
sites={'JC01','JC02','JC04'};
colors=[0 .32 .78; .88 .38 0];
data=cell(3,1);
for k=1:3, data{k}=load(fullfile(input_dir,sites{k},'satellite_plot_data.mat')); end
f=figure('Visible','on','Color','w','Position',[60 60 1500 900]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for k=1:3
    for sys=1:2
        nexttile; hold on;
        if sys==1, s=data{k}.sat_counts.gps; else, s=data{k}.sat_counts.bds; end
        t=hours(s.time-dateshift(s.time,'start','day'));
        scatter(t,s.count_1,7,colors(1,:),'filled');
        scatter(t,s.count_2,7,colors(2,:),'filled');
        xlim([0 24]); xticks(0:6:24); ylim([0 max([s.count_1;s.count_2])+2]);
        grid on; box on; set(gca,'FontSize',12);
        title([sites{k} ' / ' s.system_name]); ylabel('Number of satellites');
        if k==3, xlabel('GPST (h)'); end
        legend(s.labels,'Location','southwest','Orientation','horizontal');
    end
end
if export_figures, writefig(f,fullfile(output_dir,'fig2_satellite_counts_JC_all')); end
f=figure('Visible','on','Color','w','Position',[50 50 1800 650]);
for k=1:3
    ax=polaraxes(f,'Position',[(k-1)/3+.045 .16 .25 .70]); hold(ax,'on');
    ax.ThetaZeroLocation='top'; ax.ThetaDir='clockwise'; ax.FontSize=11;
    ax.ThetaTick=0:45:315;
    ax.ThetaTickLabel={'N','45°','E','135°','S','225°','W','315°'};
    rlim(ax,[0 90]); rticks(ax,[30 60 90]); rticklabels(ax,{'60°','30°','0°'});
    title(ax,sites{k},'FontSize',14);
    tracks=data{k}.skyplot_data.tracks;
    for j=1:numel(tracks)
        tr=tracks(j); which=find('GC'==tr.system);
        if isempty(which), continue; end
        t=tr.time(:); az=tr.az(:); el=tr.el(:);
        % Keep every observed point. NaNs prevent interpolation across gaps.
        dt=seconds(diff(t)); step=median(dt(dt>0));
        if isempty(step) || ~isfinite(step), step=30; end
        breaks=find(dt>2*step);
        theta=deg2rad(az); r=90-el;
        for b=flip(breaks(:)')
            theta=[theta(1:b);nan;theta(b+1:end)]; %#ok<AGROW>
            r=[r(1:b);nan;r(b+1:end)]; %#ok<AGROW>
        end
        polarplot(ax,theta,r,'Color',colors(which,:),'LineWidth',1.2);
    end
    hg=polarplot(ax,nan,nan,'Color',colors(1,:),'LineWidth',1.5);
    hc=polarplot(ax,nan,nan,'Color',colors(2,:),'LineWidth',1.5);
    legend(ax,[hg hc],{'GPS','BDS'},'Location','southoutside','Orientation','horizontal');
end
if export_figures, writefig(f,fullfile(output_dir,'fig3_skyplots_JC_all')); end
drawnow;
fprintf('Combined satellite-count and gap-aware skyplot figures complete.\n');
end

function writefig(f,stem)
exportgraphics(f,[stem '.png'],'Resolution',300);
exportgraphics(f,[stem '.pdf'],'ContentType','vector');
savefig(f,[stem '.fig']);
end
