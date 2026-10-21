clear; close all; clc;
root='F:/RTKLib-LAB/tmp/static_three_methods_20261011'; outdir=fullfile(root,'figures_p1e4'); if ~exist(outdir,'dir'),mkdir(outdir);end
datasets={'JC01','JC02','JC04','Static'}; methods={'slope','orig','ransac'}; labels={'Pure slope','Original ARAIM','RANSAC-ARAIM'};
cpe=[0.92 0.08 0.08]; cpl=[0.08 0.20 0.90];
for id=1:numel(datasets)
  ds=datasets{id}; D=cell(1,3);
  for im=1:3, D{im}=jsondecode(fileread(fullfile(root,sprintf('%s_%s_p1e4.json',ds,methods{im})))); end
  for st=[1 2]
    if st==1, qname='fixed'; else, qname='float'; end
    fig=figure('Color','w','Position',[80 60 1700 1200]); tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    handles=[];
    for im=1:3
      q=[D{im}.q]; if st==1, mask=(q==1); else, mask=(q~=1); end
      if ~any(mask), mask=false(size(q)); end
      if any(mask)
        t=[D{im}.sec]; t=(t-t(find(mask,1)))/3600; t=t(mask);
        hp=[D{im}.hpl]; vp=[D{im}.vpl]; he=[D{im}.hpe]; ve=[D{im}.vpe]; hp=hp(mask);vp=vp(mask);he=he(mask);ve=ve(mask);
      else
        t=[]; hp=[];vp=[];he=[];ve=[];
      end
      nexttile; plot_one(t,hp,he,cpl,cpe,sprintf('%s | H | PE %d, PL %d',labels{im},sum(isfinite(he)),sum(isfinite(hp))));
      nexttile; plot_one(t,vp,ve,cpl,cpe,sprintf('%s | V | PE %d, PL %d',labels{im},sum(isfinite(ve)),sum(isfinite(vp))));
    end
    axleg=axes(fig,'Visible','off','Position',[0.38 0.005 0.24 0.035]); hold(axleg,'on'); h1=scatter(axleg,nan,nan,38,cpe,'filled');h2=scatter(axleg,nan,nan,38,cpl,'filled');legend(axleg,[h1 h2],{'PE','PL'},'Orientation','horizontal','Location','north');
    sgtitle(sprintf('%s | %s | p_{obs}=p_{wrong fix}=10^{-4}',ds,qname),'FontWeight','bold');
    exportgraphics(fig,fullfile(outdir,sprintf('%s_%s_p1e4.png',ds,qname)),'Resolution',220); exportgraphics(fig,fullfile(outdir,sprintf('%s_%s_p1e4.pdf',ds,qname)),'ContentType','vector'); close(fig);
  end
end
fprintf('Wrote figures to %s\n',outdir);
function plot_one(t,pl,pe,cpl,cpe,ttl)
  hold on; scatter(t,pe,8,cpe,'filled','MarkerFaceAlpha',0.85); scatter(t,pl,8,cpl,'filled','MarkerFaceAlpha',0.85); grid on;box on; xlabel('Elapsed time (h)');ylabel('PE / PL (m)');title(ttl,'FontWeight','normal');
  z=[pl(:);pe(:)];z=z(isfinite(z)&z>=0);if ~isempty(z),ylim([0 1.05*max(z)]);end
end
