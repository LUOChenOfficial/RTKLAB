function results = mc_amb_ratio(debugFile, outFile, numSamples, ratioThreshold, maxGroups)
% Monte Carlo ambiguity validation from RTKLib-LAB AMB_FLOAT/AMB_COV records.
% The true integer ambiguity is zero; accepted offset is k = z1.
if nargin < 1 || isempty(debugFile), debugFile = 'F:/RTKLib-LAB/result/Static/test_ambfix.rtk_debug'; end
if nargin < 2 || isempty(outFile), outFile = 'F:/RTKLib-LAB/result/Static/mc_ratio3.mat'; end
if nargin < 3 || isempty(numSamples), numSamples = 1000; end
if nargin < 4 || isempty(ratioThreshold), ratioThreshold = 3; end
if nargin < 5 || isempty(maxGroups), maxGroups = inf; end
rng(20260921, 'twister');
groups = read_amb_groups(debugFile);
ng = min(numel(groups), maxGroups);
total = 0; accepted = 0; correct = 0;
offsetCounts = containers.Map('KeyType', 'char', 'ValueType', 'double');
groupSummary = repmat(struct('id',0,'time','','nb',0,'samples',0,'accepted',0,'correct',0,'valid',true), ng, 1);
for g = 1:ng
    Q = (groups(g).Q + groups(g).Q') / 2;
    [L,p] = chol(Q, 'lower');
    if p ~= 0
        [L,p] = chol(Q + 1e-10*max(1,trace(Q)/groups(g).nb)*eye(groups(g).nb), 'lower');
    end
    groupSummary(g).id = groups(g).id; groupSummary(g).time = groups(g).time;
    groupSummary(g).nb = groups(g).nb; groupSummary(g).samples = numSamples;
    if p ~= 0, groupSummary(g).valid = false; continue; end
    eta = L * randn(groups(g).nb, numSamples);
    ga = 0; gc = 0;
    for m = 1:numSamples
        [z1,~,J1,J2] = ils_two_best(eta(:,m), Q, ratioThreshold);
        ratio = J2 / max(J1, realmin);
        if ratio + 1e-12 >= ratioThreshold
            accepted = accepted + 1; ga = ga + 1; k = z1;
            if all(k == 0), correct = correct + 1; gc = gc + 1; end
            key = offset_key(k);
            if isKey(offsetCounts, key), offsetCounts(key) = offsetCounts(key) + 1;
            else, offsetCounts(key) = 1; end
        end
        total = total + 1;
    end
    groupSummary(g).accepted = ga; groupSummary(g).correct = gc;
    if mod(g, max(1,floor(ng/20))) == 0 || g == ng
        fprintf('group %d/%d, accepted %.4f, correct %.4f\n', g,ng,accepted/max(total,1),correct/max(total,1));
    end
end
keys = offsetCounts.keys; counts = zeros(numel(keys),1);
for i=1:numel(keys), counts(i)=offsetCounts(keys{i}); end
[counts,ix] = sort(counts,'descend'); keys=keys(ix);
offsetTable = table(keys(:), counts, counts/max(accepted,1), 'VariableNames', {'offset','count','conditionalProbability'});
results = struct('debugFile',debugFile,'numGroups',ng,'numSamplesPerGroup',numSamples, ...
    'ratioThreshold',ratioThreshold,'totalSamples',total,'acceptedSamples',accepted, ...
    'correctAcceptedSamples',correct,'acceptanceRate',accepted/max(total,1), ...
    'unconditionalCorrectRate',correct/max(total,1),'conditionalCorrectRate',correct/max(accepted,1), ...
    'conditionalWrongRate',(accepted-correct)/max(accepted,1),'groupSummary',groupSummary,'offsetTable',offsetTable);
save(outFile,'results','-v7.3');
fprintf('saved %s\n',outFile);
fprintf('total=%d accepted=%d correct=%d acceptance=%.8g conditional_correct=%.8g\n',total,accepted,correct,results.acceptanceRate,results.conditionalCorrectRate);
end

function groups = read_amb_groups(fileName)
fid=fopen(fileName,'r'); assert(fid>=0,'Cannot open debug file: %s',fileName);
tmp=containers.Map('KeyType','char','ValueType','any');
while true
    line=fgetl(fid); if ~ischar(line), break; end
    a=strsplit(strtrim(line)); if isempty(a), continue; end
    if strcmp(a{1},'AMB_FLOAT') && numel(a)>=7
        id=a{2}; nb=str2double(a{5}); idx=str2double(a{6}); val=str2double(a{7});
        if ~isKey(tmp,id), tmp(id)=struct('id',str2double(id),'time',[a{3} ' ' a{4}],'nb',nb,'a',nan(nb,1),'Q',nan(nb,nb)); end
        s=tmp(id); s.a(idx)=val; tmp(id)=s;
    elseif strcmp(a{1},'AMB_COV') && numel(a)>=8
        id=a{2}; nb=str2double(a{5}); row=str2double(a{6}); col=str2double(a{7}); val=str2double(a{8});
        if ~isKey(tmp,id), tmp(id)=struct('id',str2double(id),'time',[a{3} ' ' a{4}],'nb',nb,'a',nan(nb,1),'Q',nan(nb,nb)); end
        s=tmp(id); s.Q(row,col)=val; tmp(id)=s;
    end
end
fclose(fid); ids=tmp.keys; groups=cell2mat(values(tmp,ids));
[~,ix]=sort([groups.id]); groups=groups(ix);
groups=groups(arrayfun(@(s) all(isfinite(s.a)) && all(isfinite(s.Q(:))),groups));
end

function key=offset_key(k)
key=sprintf('%d,',k(:));
end

function [z1,z2,J1,J2]=ils_two_best(z,Q,ratioThreshold)
n=numel(z); R=chol(inv(Q)); % Q^{-1}=R'*R
a=zeros(n,1); d=0;
for j=n:-1:1
    c=z(j)+R(j,j+1:n)*(z(j+1:n)-a(j+1:n))/R(j,j);
    a(j)=round(c);
    e=R(j,j)*(z(j)-a(j))+R(j,j+1:n)*(z(j+1:n)-a(j+1:n));
    d=d+e*e;
end
radius=max(ratioThreshold*d,1e-12);
z1=zeros(n,1); z2=zeros(n,1); J1=inf; J2=inf; cur=zeros(n,1);
search(n,0);
    function search(j,partial)
        if j==0
            if partial<J1, J2=J1; z2=z1; J1=partial; z1=cur;
            elseif partial<J2 && any(cur~=z1), J2=partial; z2=cur; end
            return
        end
        s=R(j,j+1:n)*(z(j+1:n)-cur(j+1:n)); c=z(j)+s/R(j,j);
        rem=radius-partial; if rem<0, return; end
        h=sqrt(rem)/abs(R(j,j)); lo=ceil(c-h); hi=floor(c+h);
        for v=lo:hi
            cur(j)=v; e=R(j,j)*(z(j)-v)+s; search(j-1,partial+e*e);
        end
    end
end
