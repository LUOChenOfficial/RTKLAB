function plot_mc_enu(debugFile, posFile, outPng)
if nargin<1, debugFile='F:/RTKLib-LAB/result/Static/test_ambfix.rtk_debug'; end
if nargin<2, posFile='F:/RTKLib-LAB/result/Static/test_ambfix.pos'; end
if nargin<3, outPng='F:/RTKLib-LAB/result/Static/mc_enu_offsets.png'; end
p=fopen(posFile,'r'); assert(p>0); xyz=[];
while true
    s=fgetl(p); if ~ischar(s), break; end
    if isempty(strtrim(s)) || startsWith(strtrim(s),'%'), continue; end
    tok=strsplit(strtrim(s));
    if numel(tok)>=4
        xyz=[str2double(tok{2});str2double(tok{3});str2double(tok{4})];
        if all(isfinite(xyz)), break; end
    end
end
fclose(p); assert(numel(xyz)==3);
% Reference latitude/longitude from the first valid Static ECEF solution.
% This Static run is around 23.1184 deg N, 114.1272 deg E.
[lat,lon]=ecef2ll(xyz);
if ~isfinite(lat) || ~isfinite(lon)
    lat=23.1184112353444*pi/180; lon=114.127162365773*pi/180;
end
T=[-sin(lon), cos(lon), 0; -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat); cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];
f=fopen(debugFile,'r'); assert(f>0);
E=[]; N=[]; U=[]; P=[]; idx=0;
while true
    s=fgetl(f); if ~ischar(s), break; end
    a=strsplit(strtrim(s));
    if numel(a)<10 || ~strcmp(a{1},'MC_MODE'), continue; end
    nb=str2double(a{5});
    base=5+nb;
    if numel(a)<base+5, continue; end
    prob=str2double(a{base+2});
    dxyz=[str2double(a{base+3});str2double(a{base+4});str2double(a{base+5})];
    if any(~isfinite(dxyz)), continue; end
    enu=T*dxyz; idx=idx+1;
    E(idx,1)=enu(1); N(idx,1)=enu(2); U(idx,1)=enu(3); P(idx,1)=prob;
end
fclose(f);
x=(1:numel(E))';
figure('Color','w','Position',[100 100 1200 700]);
subplot(3,1,1); scatter(x,E,10,log10(max(P,1e-12)),'filled'); grid on; ylabel('E offset (m)'); title('MC integer-mode position offsets in ENU');
subplot(3,1,2); scatter(x,N,10,log10(max(P,1e-12)),'filled'); grid on; ylabel('N offset (m)');
subplot(3,1,3); scatter(x,U,10,log10(max(P,1e-12)),'filled'); grid on; ylabel('U offset (m)'); xlabel('MC_MODE record index');
colormap turbo; cb=colorbar; cb.Label.String='log_{10}(joint probability)';
exportgraphics(gcf,outPng,'Resolution',180);
fprintf('records=%d, E=[%.4g,%.4g] N=[%.4g,%.4g] U=[%.4g,%.4g]\n',numel(E),min(E),max(E),min(N),max(N),min(U),max(U));
end
function [lat,lon]=ecef2ll(x)
a=6378137; e2=6.6943799901413165e-3;
lon=atan2(x(2),x(1)); p=hypot(x(1),x(2)); lat=atan2(x(3),p*(1-e2));
for k=1:8
    v=a/sqrt(1-e2*sin(lat)^2);
    lat=atan2(x(3)+e2*v*sin(lat),p);
end
end
