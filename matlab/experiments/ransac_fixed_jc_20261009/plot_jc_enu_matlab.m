% JC01/JC02/JC04：两边开启RANSAC和连续固定的ENU误差，按解状态显示。
% 直接读取现有POS；不重新解算，不平滑、不删除异常点。
% 仅比较共同输出历元；Q=2浮点，Q=1固定，Q=4差分。
% 参考坐标从result/JCxx/plot_enu_error.m读取，并非独立测量真值。
clear; close all; clc;
baseDir = 'F:/RTKLib-LAB/tmp/ransac_fixed_jc_20261009';
repoDir = fileparts(fileparts(baseDir));
datasets = {'JC01','JC02','JC04'};
labels = {'原版 RANSAC-RAIM','当前 RANSAC+ARAIM'};
colors = [0.82 0.48 0.22; 0.16 0.42 0.61];
qualityCodes = [2 1 4 1];
qualityNames = {'浮点解','固定解','差分解 Q=4','两边同时固定'};
suffixes = {'float','fixed','dgps','both_fixed'};
exportFigures = true; % false：只显示窗口；true：另存PNG/PDF/可编辑FIG
data = cell(1,numel(datasets));
for k = 1:numel(datasets)
    ds = datasets{k};
    refSource = fileread(fullfile(repoDir,'matlab','result',ds,'plot_enu_error.m'));
    token = regexp(refSource,'truthEcef\s*=\s*\[([^\]]+)\]','tokens','once');
    assert(~isempty(token),'未找到%s参考坐标',ds);
    ref = sscanf(strrep(token{1},',',' '),'%f').';
    assert(numel(ref)==3,'参考坐标必须是XYZ三维');
    a = readPos(fullfile(baseDir,ds,'original.pos'));
    b = readPos(fullfile(baseDir,ds,'current.pos'));
    [time,ia,ib] = intersect(a.time,b.time,'sorted');
    E = enuRotation(ref);
    item.time = time;
    item.hour = hours(time-datetime(2024,10,29)); % POS为GPST，不转换为本地时间
    item.enu = {(a.xyz(ia,:)-ref)*E.',(b.xyz(ib,:)-ref)*E.'};
    item.q = {a.q(ia),b.q(ib)};
    data{k} = item;
    fprintf('%s common=%d; original float/fix/dgps=%d/%d/%d; current=%d/%d/%d\n', ...
        ds,numel(time),nnz(item.q{1}==2),nnz(item.q{1}==1),nnz(item.q{1}==4), ...
        nnz(item.q{2}==2),nnz(item.q{2}==1),nnz(item.q{2}==4));
end
for s = 1:numel(qualityCodes)
    qCode = qualityCodes(s);
    fig = figure('Color','w','Name',qualityNames{s},'Position',[70 70 1450 850]);
    layout = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    legendHandles = gobjects(1,2);
    directions = {'东 E','北 N','天 U'};
    for row = 1:3
        for col = 1:numel(datasets)
            ax = nexttile(layout,(row-1)*3+col); hold(ax,'on');
            ax.Toolbar.Visible = 'off'; % 导出时避免悬浮工具栏遮住标题
            item = data{col}; counts = cellfun(@(q)nnz(q==qCode),item.q);
            if s==4
                bothFixed = item.q{1}==1 & item.q{2}==1;
                counts = [nnz(bothFixed),nnz(bothFixed)];
            end
            rmsValues = nan(1,2);
            for method = 1:2
                mask = item.q{method}==qCode;
                if s==4, mask = bothFixed; end
                y = item.enu{method}(:,row); y(~mask) = NaN;
                [x,yplot] = breakGaps(item.hour,y);
                h = plot(ax,x,yplot,'Color',colors(method,:),'LineWidth',0.7, ...
                    'Marker','.','MarkerSize',3,'DisplayName',labels{method});
                if row==1 && col==1, legendHandles(method)=h; end
                if any(mask), rmsValues(method)=sqrt(mean(item.enu{method}(mask,row).^2)); end
            end
            yline(ax,0,'--','Color',[.5 .5 .5],'HandleVisibility','off');
            grid(ax,'on'); ax.GridAlpha=.15; ax.FontName='Microsoft YaHei'; ax.FontSize=10;
            xlim(ax,[0 24]); xticks(ax,0:6:24);
            ylabel(ax,[directions{row} '误差 (m)']);
            if row==1
                title(ax,sprintf('%s · %s历元：%d / %d',datasets{col},qualityNames{s},counts), ...
                    'FontSize',11,'FontWeight','normal');
            end
            if row==3, xlabel(ax,'GPST 时刻 (h)'); end
            if sum(counts)==0
                text(ax,.5,.5,['共同输出历元中两边均无' qualityNames{s}], ...
                    'Units','normalized','HorizontalAlignment','center','FontSize',10);
            else
                text(ax,.02,.97,sprintf('RMS: %.3f / %.3f m',rmsValues), ...
                    'Units','normalized','VerticalAlignment','top','BackgroundColor','w','FontSize',9);
            end
            fprintf('%s %s %s RMS original/current=%.6f/%.6f m\n', ...
                datasets{col},qualityNames{s},directions{row},rmsValues);
        end
    end
    title(layout,['RANSAC + 连续固定 · ' qualityNames{s} ' · 共同输出历元 · 2024-10-29'], ...
        'FontName','Microsoft YaHei','FontSize',15);
    lgd = legend(legendHandles,labels,'Orientation','horizontal','Box','off','FontSize',11);
    lgd.Layout.Tile = 'south';
    drawnow;
    if exportFigures
        output = fullfile(baseDir,['JC_enu_matlab_' suffixes{s}]);
        exportgraphics(fig,[output '.png'],'Resolution',300);
        exportgraphics(fig,[output '.pdf'],'ContentType','vector');
        savefig(fig,[output '.fig']);
    end
end

function sol = readPos(path)
    fid = fopen(path,'r'); assert(fid>=0,'无法打开%s',path);
    cleanup = onCleanup(@()fclose(fid));
    stamps = strings(0,1); xyz=zeros(0,3); q=zeros(0,1);
    while true
        line=fgetl(fid); if ~ischar(line), break; end
        line=strtrim(line);
        if isempty(line)||startsWith(line,'%')||startsWith(line,'#'), continue; end
        p=split(string(line)); assert(numel(p)>=7,'POS列数不足：%s',path);
        stamps(end+1,1)=p(1)+" "+p(2); %#ok<AGROW>
        xyz(end+1,:)=str2double(p(3:5)).'; %#ok<AGROW>
        q(end+1,1)=str2double(p(6)); %#ok<AGROW>
    end
    sol.time=datetime(stamps,'InputFormat','yyyy/MM/dd HH:mm:ss.SSS');
    sol.xyz=xyz; sol.q=q;
    assert(all(isfinite(xyz),'all'),'POS包含非有限XYZ');
end

function E = enuRotation(ref)
    e2=6.6943799901413165e-3; a=6378137;
    lon=atan2(ref(2),ref(1)); p=hypot(ref(1),ref(2));
    lat=atan2(ref(3),p*(1-e2));
    for i=1:12
        n=a/sqrt(1-e2*sin(lat)^2);
        lat=atan2(ref(3)+e2*n*sin(lat),p);
    end
    E=[-sin(lon),cos(lon),0; ...
       -sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat); ...
       cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)];
end

function [xout,yout] = breakGaps(x,y)
    % 缺失历元或解状态变化时不跨段连接；孤立历元由点标记保留。
    xout=[]; yout=[];
    for i=1:numel(x)
        if i>1 && x(i)-x(i-1)>30.1/3600
            xout(end+1,1)=NaN; yout(end+1,1)=NaN; %#ok<AGROW>
        end
        xout(end+1,1)=x(i); yout(end+1,1)=y(i); %#ok<AGROW>
    end
end
