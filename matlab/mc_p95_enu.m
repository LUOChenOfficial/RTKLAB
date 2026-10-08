function out = mc_p95_enu(debugFile, outCsv, outPng, posFile, quantile)
if nargin<1, debugFile='F:/RTKLib-LAB/result/Open/test_ambfix.rtk_debug'; end
if nargin<2, outCsv='F:/RTKLib-LAB/result/Open/mc_p95_enu.csv'; end
if nargin<3, outPng='F:/RTKLib-LAB/result/Open/mc_p95_enu.png'; end
if nargin<4, posFile=strrep(debugFile,'test_ambfix.rtk_debug','test_ambfix.pos'); end
if nargin<5, quantile=0.95; end
assert(quantile>0 && quantile<1, 'quantile must be in (0,1)');
qLabel=round(100*quantile);
fid=fopen(debugFile,'r'); assert(fid>0);
ids={}; t={}; E={}; N={}; U={}; W={}; nbv=[];
while true
    line=fgetl(fid); if ~ischar(line), break; end
    a=strsplit(strtrim(line));
    if numel(a)<10 || ~strcmp(a{1},'MC_MODE'), continue; end
    id=a{2}; nb=str2double(a{5}); base=6+nb;
    if numel(a)<base+5, continue; end
    idx=find(strcmp(ids,id),1);
    if isempty(idx)
        idx=numel(ids)+1; ids{idx}=id; t{idx}=[a{3} ' ' a{4}]; nbv(idx)=nb;
        E{idx}=[]; N{idx}=[]; U{idx}=[]; W{idx}=[];
    end
    E{idx}(end+1)=str2double(a{base+3});
    N{idx}(end+1)=str2double(a{base+4});
    U{idx}(end+1)=str2double(a{base+5});
    W{idx}(end+1)=str2double(a{base+1});
end
fclose(fid);
n=numel(ids); out=table('Size',[n 9],'VariableTypes',{'string','string','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'id','time','nb','E_bound','N_bound','U_bound','norm_bound','acceptedModes','acceptedWeight'});
for i=1:n
    w=W{i}; w=w/sum(w); e=abs(E{i}); nn=abs(N{i}); u=abs(U{i}); d=sqrt(E{i}.^2+N{i}.^2+U{i}.^2);
    out.id(i)=string(ids{i}); out.time(i)=string(t{i}); out.nb(i)=nbv(i);
    out.E_bound(i)=wquant(e,w,quantile); out.N_bound(i)=wquant(nn,w,quantile);
    out.U_bound(i)=wquant(u,w,quantile); out.norm_bound(i)=wquant(d,w,quantile);
    out.acceptedModes(i)=numel(w); out.acceptedWeight(i)=sum(W{i});
end
writetable(out,outCsv);
x=(1:n)';
allTimes={}; pf=fopen(posFile,'r');
if pf>0
    while true
        s=fgetl(pf); if ~ischar(s), break; end
        s=strtrim(s); if isempty(s)||startsWith(s,'%'), continue; end
        a=strsplit(s); if numel(a)>=2, allTimes{end+1}=[a{1} ' ' a{2}]; end
    end
    fclose(pf);
end
if ~isempty(allTimes)
    xt=nan(numel(allTimes),1); fixed=containers.Map(out.time,1:height(out));
    for j=1:numel(allTimes), if isKey(fixed,allTimes{j}), xt(j)=fixed(allTimes{j}); end, end
    xfull=1:numel(allTimes); ye=nan(size(xfull)); yn=ye; yu=ye;
    for j=1:numel(allTimes), if isfinite(xt(j)), ye(j)=out.E_bound(xt(j)); yn(j)=out.N_bound(xt(j)); yu(j)=out.U_bound(xt(j)); end, end
else
    xfull=x; ye=out.E_bound; yn=out.N_bound; yu=out.U_bound;
end
% Exclude the first ten epochs used for filter initialization.
skip=10;
if numel(xfull)>skip
    xfull=xfull(skip+1:end); ye=ye(skip+1:end); yn=yn(skip+1:end); yu=yu(skip+1:end);
end
ye=100*ye; yn=100*yn; yu=100*yu;
figure('Color','w','Position',[120 80 1150 720]);
tl=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; scatter(xfull,ye,16,[0.05 0.35 0.75],'filled'); grid on; box on; ylabel(sprintf('E_{%d} (cm)',qLabel)); title(sprintf('%d%% ENU offset bound per epoch',qLabel));
nexttile; scatter(xfull,yn,16,[0.85 0.25 0.10],'filled'); grid on; box on; ylabel(sprintf('N_{%d} (cm)',qLabel));
nexttile; scatter(xfull,yu,16,[0.15 0.55 0.25],'filled'); grid on; box on; ylabel(sprintf('U_{%d} (cm)',qLabel)); xlabel('Epoch index');
set(findall(gcf,'-property','FontName'),'FontName','Arial');
exportgraphics(gcf,outPng,'Resolution',180);
fprintf('epochs=%d E%d max=%.4g N%d max=%.4g U%d max=%.4g norm%d max=%.4g\n',n,qLabel,max(out.E_bound),qLabel,max(out.N_bound),qLabel,max(out.U_bound),qLabel,max(out.norm_bound));
end
function q=wquant(x,w,p)
[x,ix]=sort(x(:)); w=w(ix); c=cumsum(w)/sum(w); j=find(c>=p,1,'first'); q=x(j);
end
