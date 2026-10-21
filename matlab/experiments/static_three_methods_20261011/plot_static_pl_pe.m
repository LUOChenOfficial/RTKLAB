clear; close all; clc;
root = 'F:/RTKLib-LAB/tmp/static_three_methods_20261011';
names = {'pure_gls','ransac'};
labels = {'Original ARAIM','RANSAC-ARAIM'};
colors.pl = [0.10 0.45 0.85];
colors.pe = [0.85 0.20 0.20];
[refTime, refEcef] = readKsxt(fullfile(root,'..','..','result','Static','COM25-comB-1.dat'));
data = cell(1,2);
for k = 1:2
    data{k} = readMethod(fullfile(root,[names{k} '.pos']), fullfile(root,[names{k} '.pint']), refTime, refEcef);
end
fig = figure('Color','w','Position',[80 80 1450 780]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for k = 1:2
    pln=sum(isfinite(data{k}.hpl)&isfinite(data{k}.vpl));
    nexttile; plotPanel(data{k}.elapsed,data{k}.hpl,data{k}.hpe,colors,sprintf('%s | H | PE %d, PL %d',labels{k},sum(isfinite(data{k}.hpe)),pln),'PE / PL (m)');
    nexttile; plotPanel(data{k}.elapsed,data{k}.vpl,data{k}.vpe,colors,sprintf('%s | V | PE %d, PL %d',labels{k},sum(isfinite(data{k}.vpe)),pln),'PE / PL (m)');
end
axleg=axes(fig,'Visible','off','Position',[0.35 0.01 0.30 0.04]); hold(axleg,'on'); h1=scatter(axleg,nan,nan,38,colors.pe,'filled'); h2=scatter(axleg,nan,nan,38,colors.pl,'filled'); legend(axleg,[h1 h2],{'PE','PL'},'Orientation','horizontal','Location','north');
exportgraphics(fig,fullfile(root,'Static_ARAIM_RANSAC_PL_PE.png'),'Resolution',220);
exportgraphics(fig,fullfile(root,'Static_ARAIM_RANSAC_PL_PE.pdf'),'ContentType','vector');
fprintf('Wrote %s\n',fullfile(root,'Static_ARAIM_RANSAC_PL_PE.png'));

function plotPanel(t,pl,pe,c,titleText,ylabelText)
    hold on; scatter(t,pe,8,c.pe,'filled','MarkerFaceAlpha',0.85); scatter(t,pl,8,c.pl,'filled','MarkerFaceAlpha',0.85);
    grid on; box on; xlabel('Elapsed time (h)'); ylabel(ylabelText); title(titleText,'FontWeight','normal');
    finite=[pl(:);pe(:)]; finite=finite(isfinite(finite) & finite>=0); if ~isempty(finite), ylim([0 1.05*max(finite)]); end
end

function d=readMethod(posPath,pintPath,refTime,refEcef)
    P=readNumeric(posPath,20); I=readNumeric(pintPath,20); n=min(size(P,1),size(I,1)); P=P(1:n,:); I=I(1:n,:);
    d.hpl=I(:,10); d.vpl=I(:,11); d.hpe=nan(n,1); d.vpe=nan(n,1); T=readTimes(posPath,n); d.elapsed=hours(T-T(1));
    for i=1:n
        [~,j]=min(abs(seconds(T(i)-refTime))); E=ecef2enu(P(i,1:3).'-refEcef(j,:).',refEcef(j,:)); d.hpe(i)=norm(E(1:2)); d.vpe(i)=abs(E(3));
    end
end

function A=readNumeric(path,ncol)
    fid=fopen(path,'r'); fmt=['%s %s ' repmat('%f ',1,ncol-2)]; C=textscan(fid,fmt,'CommentStyle','%','CollectOutput',true); fclose(fid); A=C{2};
end

function T=readTimes(path,n)
    fid=fopen(path,'r'); a=textscan(fid,'%s %s %*[^\n]','CommentStyle','%'); fclose(fid);
    T=datetime(strcat(a{1},{' '},a{2}),'InputFormat','yyyy/MM/dd HH:mm:ss.SSS'); T=T(1:min(n,numel(T)));
end

function [T,X]=readKsxt(path)
    fid=fopen(path,'r'); T=datetime.empty; X=zeros(0,3);
    while true
        s=fgetl(fid); if ~ischar(s), break; end; if ~startsWith(s,'$KSXT,'), continue; end
        a=split(s,','); if numel(a)<11, continue; end; if str2double(a{11})~=3, continue; end
        ts=a{2}; T(end+1,1)=datetime(ts(1:14),'InputFormat','yyyyMMddHHmmss')+seconds(str2double(ts(15:end)));
        X(end+1,:)=llh2ecef(str2double(a{4}),str2double(a{3}),str2double(a{5}));
    end; fclose(fid);
end

function x=llh2ecef(lat,lon,h)
    a=6378137; e2=6.6943799901413165e-3; lat=lat*pi/180; lon=lon*pi/180; N=a/sqrt(1-e2*sin(lat)^2);
    x=[(N+h)*cos(lat)*cos(lon),(N+h)*cos(lat)*sin(lon),(N*(1-e2)+h)*sin(lat)];
end

function enu=ecef2enu(dx,ref)
    x=ref(1); y=ref(2); z=ref(3); lon=atan2(y,x); p=hypot(x,y); e2=6.6943799901413165e-3; lat=atan2(z,p*(1-e2));
    for k=1:12, N=6378137/sqrt(1-e2*sin(lat)^2); lat=atan2(z+e2*N*sin(lat),p); end
    sl=sin(lon); cl=cos(lon); sp=sin(lat); cp=cos(lat); R=[-sl cl 0;-sp*cl -sp*sl cp;cp*cl cp*sl sp]; enu=R*dx;
end
