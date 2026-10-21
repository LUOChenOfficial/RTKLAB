function plot_jc_errors
% Preview stored JC solutions without rerunning the positioning solver.
root = 'F:/RTKLib-LAB';
out = 'F:/RTKLib-LAB/tmp/jc_error_preview_20261012';
sites = {'JC01','JC02','JC04'};
truth = [-1438180.4067 4627043.3969 4135258.3517; ...
         -1438238.7490 4627034.2716 4135248.5434; ...
         -1437589.8951 4627633.0746 4134807.3448];
components = {'East','North','Up'};
colors = [0.15 0.40 0.70; 0.80 0.38 0.16; 0.23 0.57 0.40];
data = cell(3,1);
for s = 1:3
    fid = fopen(fullfile(root,'result',sites{s},'test.pos'),'r');
    assert(fid >= 0, 'Cannot open input');
    cleaner = onCleanup(@() fclose(fid));
    time = []; values = [];
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line,'%') || startsWith(line,'#'), continue; end
        parts = strsplit(line);
        assert(numel(parts)>=20, 'Missing float coordinate columns');
        time(end+1,1) = datenum([parts{1} ' ' parts{2}], 'yyyy/mm/dd HH:MM:SS.FFF'); %#ok<AGROW>
        values(end+1,:) = str2double(parts(3:20)); %#ok<AGROW>
    end
    clear cleaner
    xyz = truth(s,:);
    lon = atan2(xyz(2),xyz(1)); p = hypot(xyz(1),xyz(2));
    e2 = 6.6943799901413165e-3; a = 6378137;
    lat = atan2(xyz(3),p*(1-e2));
    for k = 1:12
        n = a/sqrt(1-e2*sin(lat)^2);
        lat = atan2(xyz(3)+e2*n*sin(lat),p);
    end
    rot = [-sin(lon) cos(lon) 0; ...
           -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); ...
           cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
    d.hours = (time-floor(time(1)))*24;
    d.fixed = (values(:,1:3)-xyz)*rot';
    d.fixed(values(:,4)~=1,:) = NaN;
    d.float = (values(:,16:18)-xyz)*rot';
    validFloat = all(isfinite(values(:,16:18)),2) & vecnorm(values(:,16:18),2,2)>1e6;
    d.float(~validFloat,:) = NaN;
    data{s} = d;
    fprintf('%s: total=%d fixed=%d float_coordinates=%d\n',sites{s},numel(time),sum(values(:,4)==1),sum(validFloat));
    fprintf('  fixed ENU RMS (cm): %.4f %.4f %.4f\n',sqrt(mean(d.fixed.^2,1,'omitnan'))*100);
    fprintf('  float ENU RMS (m): %.4f %.4f %.4f\n',sqrt(mean(d.float.^2,1,'omitnan')));
end
for mode = {'fixed','float'}
    name = mode{1};
    if strcmp(name,'fixed'), scale=100; unit='cm'; titleText='Fixed solution (Q = 1)';
    else, scale=1; unit='m'; titleText='Float estimates (all available epochs)'; end
    for s = 1:3
        allErrors = data{s}.(name)*scale;
        maxAbs = max(abs(allErrors(:)),[],'omitnan');
        % Share a symmetric range across ENU only within this station and mode.
        magnitude = 10^floor(log10(max(maxAbs,eps)));
        axisLimit = ceil(1.05*maxAbs/magnitude)*magnitude;
        f = figure('Visible','off','Color','w','Position',[80 80 1050 720]);
        layout = tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
        title(layout,[sites{s} ' | ' titleText],'FontName','Arial','FontSize',17,'FontWeight','bold');
        errors = data{s}.(name)*scale;
        for c = 1:3
            ax=nexttile(layout);
            plot(ax,data{s}.hours,errors(:,c),'.-','Color',colors(c,:), ...
                'LineWidth',0.6,'MarkerSize',3);
            yline(ax,0,':','Color',[0.5 0.5 0.5]);
            xlim(ax,[0 24]); xticks(ax,0:4:24); grid(ax,'on');
            ylim(ax,[-axisLimit axisLimit]);
            ax.FontName='Arial'; ax.FontSize=11;
            ylabel(ax,sprintf('%s error (%s)',components{c},unit));
            if c==3, xlabel(ax,'GPST hour, 2024-10-29'); end
        end
        exportgraphics(f,fullfile(out,[sites{s} '_ENU_' name '.png']),'Resolution',160);
        savefig(f,fullfile(out,[sites{s} '_ENU_' name '.fig']));
        close(f);
    end
end
end
