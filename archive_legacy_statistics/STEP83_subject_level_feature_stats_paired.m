%% STEP83_subject_level_feature_stats_paired.m

cfg = get_project_config();

clear; clc;

cfg = struct();
cfg.ROOT = cfg.outputRoot;
cfg.datasetPath = '';
cfg.candidatePaths = { ...
    fullfile(cfg.ROOT,'_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'Subjects','_wm_ml','dataset.mat')
    fullfile(cfg.ROOT,'_wm_dataset','dataset.mat')
    fullfile(cfg.ROOT,'dataset.mat')
};
cfg.montageFile = fullfile(cfg.ROOT,'ipm2.ced.tsv');
cfg.outDir = fullfile(cfg.ROOT,'_wm_STEP83_subject_level_feature_stats');

cfg.runList = [1 2 3];
cfg.phaseList = {'stim','maint','retr'};
cfg.aggregateMethod = 'median';
cfg.alphaFDR = 0.05;
cfg.minSubjectsBinary = 8;
cfg.minSubjectsThreeClass = 8;

cfg.labelMap.color       = {'1','color'};
cfg.labelMap.orientation = {'2','orientation'};
cfg.labelMap.conjunction = {'3','conjunction'};

contrasts = struct([]);
contrasts(1).name='orientation_vs_conjunction'; contrasts(1).A=cfg.labelMap.orientation; contrasts(1).B=cfg.labelMap.conjunction;
contrasts(2).name='color_vs_orientation';       contrasts(2).A=cfg.labelMap.color;       contrasts(2).B=cfg.labelMap.orientation;
contrasts(3).name='color_vs_conjunction';       contrasts(3).A=cfg.labelMap.color;       contrasts(3).B=cfg.labelMap.conjunction;

if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
diary(fullfile(cfg.outDir,'STEP83_subject_level_feature_stats_log.txt'));

fprintf('\n=== STEP83 subject-level paired feature statistics ===\n');
fprintf('Started: %s\n', datestr(now));

[T,col,featureNames] = local_load_data(cfg);
labels = local_load_labels(cfg.montageFile);

subjVals  = cleanstr(T.(col.subj));
condVals  = cleanstr(T.(col.cond));
runVals   = double(T.(col.run));
phaseVals = cleanstr(T.(col.phase));

allBin = {};
allFri = {};
jobs = {};

for ir = 1:numel(cfg.runList)
    runNum = cfg.runList(ir);
    for ip = 1:numel(cfg.phaseList)
        phaseName = cfg.phaseList{ip};
        idxBase = runVals==runNum & phaseVals==cleanstr({phaseName});

        fprintf('\n--- run %d | phase %s | rows=%d ---\n', runNum, phaseName, sum(idxBase));

        if sum(idxBase) < 20, continue; end

        %% Friedman 3-class: color/orientation/conjunction
        [X3, subj3] = subj_cond_matrix(T,idxBase,subjVals,condVals,featureNames, ...
            {cfg.labelMap.color,cfg.labelMap.orientation,cfg.labelMap.conjunction}, cfg.aggregateMethod);

        if numel(subj3) >= cfg.minSubjectsThreeClass
            F = run_friedman(X3, featureNames);
            F.Test = repmat("Friedman_3class_subjectLevel",height(F),1);
            F.Contrast = repmat("color_vs_orientation_vs_conjunction",height(F),1);
            F.Run = repmat(runNum,height(F),1);
            F.Phase = repmat(string(phaseName),height(F),1);
            F.N_subjects = repmat(numel(subj3),height(F),1);
            F = movevars(F,{'Test','Contrast','Run','Phase','N_subjects'},'Before',1);
            F = add_localization(F,labels);
            allFri{end+1,1}=F; %#ok<AGROW>

            nSig = sum(F.q_FDR < cfg.alphaFDR,'omitnan');
            minq = min(F.q_FDR,[],'omitnan');
            maxe = max(F.KendallW,[],'omitnan');
            jobs{end+1,1}=table(string("Friedman_3class_subjectLevel"),string("color_vs_orientation_vs_conjunction"),runNum,string(phaseName),numel(subj3),height(F),nSig,minq,maxe, ...
                'VariableNames',{'Test','Contrast','Run','Phase','N_subjects','N_features','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>
            save_outputs(cfg,'Friedman_3class','color_vs_orientation_vs_conjunction',runNum,phaseName,F,'friedman');
            fprintf('Friedman: subjects=%d | Nsig=%d | min q=%.3g\n',numel(subj3),nSig,minq);
        end

        %% Binary paired signrank
        for ic = 1:numel(contrasts)
            C = contrasts(ic);
            [X2, subj2] = subj_cond_matrix(T,idxBase,subjVals,condVals,featureNames,{C.A,C.B},cfg.aggregateMethod);

            if numel(subj2) < cfg.minSubjectsBinary
                fprintf('Signrank %s skipped: subjects=%d\n', C.name, numel(subj2));
                continue;
            end

            B = run_signrank(X2, featureNames);
            B.Test = repmat("Signrank_binary_subjectLevel",height(B),1);
            B.Contrast = repmat(string(C.name),height(B),1);
            B.Run = repmat(runNum,height(B),1);
            B.Phase = repmat(string(phaseName),height(B),1);
            B.N_subjects = repmat(numel(subj2),height(B),1);
            B = movevars(B,{'Test','Contrast','Run','Phase','N_subjects'},'Before',1);
            B = add_localization(B,labels);
            allBin{end+1,1}=B; %#ok<AGROW>

            nSig = sum(B.q_FDR < cfg.alphaFDR,'omitnan');
            minq = min(B.q_FDR,[],'omitnan');
            maxe = max(abs(B.RankBiserial_paired),[],'omitnan');
            jobs{end+1,1}=table(string("Signrank_binary_subjectLevel"),string(C.name),runNum,string(phaseName),numel(subj2),height(B),nSig,minq,maxe, ...
                'VariableNames',{'Test','Contrast','Run','Phase','N_subjects','N_features','N_sig_FDR','Min_q','MaxEffect'}); %#ok<AGROW>
            save_outputs(cfg,'Signrank_binary',C.name,runNum,phaseName,B,'signrank');
            fprintf('Signrank %s: subjects=%d | Nsig=%d | min q=%.3g\n',C.name,numel(subj2),nSig,minq);
        end
    end
end

if isempty(allBin), Binary_All=table(); else, Binary_All=vertcat(allBin{:}); end
if isempty(allFri), Friedman_All=table(); else, Friedman_All=vertcat(allFri{:}); end
if isempty(jobs), JobSummary=table(); else, JobSummary=sortrows(vertcat(jobs{:}),{'N_sig_FDR','Min_q'},{'descend','ascend'}); end

writetable(Binary_All,fullfile(cfg.outDir,'STEP83_binary_subject_level_signrank_all_features.csv'));
writetable(Friedman_All,fullfile(cfg.outDir,'STEP83_threeclass_subject_level_friedman_all_features.csv'));
writetable(JobSummary,fullfile(cfg.outDir,'STEP83_subject_level_job_summary.csv'));
save(fullfile(cfg.outDir,'STEP83_subject_level_feature_stats.mat'),'Binary_All','Friedman_All','JobSummary','cfg','-v7.3');

fprintf('\nSaved job summary:\n%s\n',fullfile(cfg.outDir,'STEP83_subject_level_job_summary.csv'));
fprintf('Finished: %s\n',datestr(now));
diary off;

%% ===================== functions =====================

function [T,col,featureNames] = local_load_data(cfg)
datasetPath='';
for i=1:numel(cfg.candidatePaths)
    if exist(cfg.candidatePaths{i},'file'), datasetPath=cfg.candidatePaths{i}; break; end
end
if isempty(datasetPath), error('dataset.mat not found.'); end
fprintf('Loading: %s\n',datasetPath);
S=load(datasetPath);
if isfield(S,'DS') && istable(S.DS)
    T=S.DS;
elseif isfield(S,'T') && istable(S.T)
    T=S.T;
else
    error('Could not find table DS or T.');
end
names=T.Properties.VariableNames;
col.subj=findvar(names,{'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID'});
col.run=findvar(names,{'runNum','RunNum','run','Run'});
col.phase=findvar(names,{'phase','Phase'});
col.cond=findvar(names,{'Condition','condition','yCondition','Label','label'});
if isempty(col.subj)||isempty(col.run)||isempty(col.phase)||isempty(col.cond), error('Missing metadata columns.'); end
featureNames=detect_features(T);
end

function v=findvar(names,cands)
v='';
for i=1:numel(cands)
    if any(strcmp(names,cands{i})), v=names{find(strcmp(names,cands{i}),1)}; return; end
end
end

function featureNames=detect_features(T)
names=T.Properties.VariableNames;
isNum=false(1,numel(names));
for j=1:numel(names), x=T.(names{j}); isNum(j)=isnumeric(x)||islogical(x); end
numNames=names(isNum);
meta={'Subject','subject','SubjectID','subjectID','subj','subjID','Subj','SubjID','Run','run','runNum','RunNum','session','Session','Phase','phase','Condition','condition','yCondition','Label','label','Correct','correct','yCorrect','TrialNum','trialNum','TrialIndex','trialIndex','Trial','trial','PatternID','patternID','StartRow','EndRow','Second10Row','RetrRow','Fold','fold','CVFold'};
isMeta=false(size(numNames));
for j=1:numel(numNames), isMeta(j)=any(strcmpi(numNames{j},meta)); end
featureNames=numNames(~isMeta);
end

function s=cleanstr(x)
if isnumeric(x)||islogical(x), s=string(x);
elseif iscell(x), s=string(x);
elseif iscategorical(x), s=string(x);
elseif isstring(x), s=x;
elseif ischar(x), s=string(cellstr(x));
else, s=string(x); end
s=lower(strtrim(s)); s=regexprep(s,'\s+',''); s=regexprep(s,'[^\w]','');
end

function [Xsubj,subjKeep] = subj_cond_matrix(T,idxBase,subjVals,condVals,featureNames,classSets,agg)
subs=unique(subjVals(idxBase)); subs=subs(~ismissing(subs));
nC=numel(classSets); nF=numel(featureNames); Xlist={}; subjKeep=strings(0,1);
for s=1:numel(subs)
    idxS=idxBase & subjVals==subs(s);
    Xi=NaN(nC,nF); ok=true;
    for c=1:nC
        idxC=idxS & ismember(condVals,cleanstr(classSets{c}));
        if sum(idxC)<1, ok=false; break; end
        X=double(table2array(T(idxC,featureNames))); X(~isfinite(X))=NaN;
        if strcmpi(agg,'mean'), Xi(c,:)=mean(X,1,'omitnan'); else, Xi(c,:)=median(X,1,'omitnan'); end
    end
    if ok, Xlist{end+1,1}=Xi; subjKeep(end+1,1)=subs(s); end %#ok<AGROW>
end
Xsubj=NaN(numel(Xlist),nC,nF);
for s=1:numel(Xlist), Xsubj(s,:,:)=Xlist{s}; end
end

function B=run_signrank(Xsubj,featureNames)
nF=size(Xsubj,3);
p=NaN(nF,1); z=NaN(nF,1); rb=NaN(nF,1); medA=NaN(nF,1); medB=NaN(nF,1); medD=NaN(nF,1); nsub=NaN(nF,1);
for j=1:nF
    A=squeeze(Xsubj(:,1,j)); Bv=squeeze(Xsubj(:,2,j)); ok=isfinite(A)&isfinite(Bv);
    A=A(ok); Bv=Bv(ok); d=A-Bv; d=d(isfinite(d)&d~=0); nsub(j)=numel(d);
    if numel(d)<5 || std(d,0,'omitnan')<=eps, continue; end
    try
        [p(j),~,stats]=signrank(A,Bv);
        if isfield(stats,'zval'), z(j)=stats.zval; end
        ranks=tiedrank(abs(d)); Wp=sum(ranks(d>0)); Wn=sum(ranks(d<0)); rb(j)=(Wp-Wn)/max(Wp+Wn,eps);
        medA(j)=median(A,'omitnan'); medB(j)=median(Bv,'omitnan'); medD(j)=median(A-Bv,'omitnan');
    catch
    end
end
q=bh(p);
B=table(string(featureNames(:)),family(string(featureNames(:))),p,q,z,rb,medA,medB,medD,nsub, ...
    'VariableNames',{'Feature','Family','p_signrank','q_FDR','z_signrank','RankBiserial_paired','Median_A_subjectAgg','Median_B_subjectAgg','Median_Diff_AminusB','N_subjects_feature'});
end

function F=run_friedman(Xsubj,featureNames)
nF=size(Xsubj,3);
p=NaN(nF,1); chi=NaN(nF,1); W=NaN(nF,1); m1=NaN(nF,1); m2=NaN(nF,1); m3=NaN(nF,1); nsub=NaN(nF,1);
for j=1:nF
    M=squeeze(Xsubj(:,:,j)); ok=all(isfinite(M),2); M=M(ok,:); n=size(M,1); nsub(j)=n;
    if n<5 || std(M(:),0,'omitnan')<=eps, continue; end
    try
        [p(j),tbl]=friedman(M,1,'off');
        chi(j)=tbl{2,5}; W(j)=chi(j)/max(n*(3-1),1);
        m1(j)=median(M(:,1),'omitnan'); m2(j)=median(M(:,2),'omitnan'); m3(j)=median(M(:,3),'omitnan');
    catch
    end
end
q=bh(p);
F=table(string(featureNames(:)),family(string(featureNames(:))),p,q,chi,W,m1,m2,m3,nsub, ...
    'VariableNames',{'Feature','Family','p_friedman','q_FDR','Friedman_chi2','KendallW','Median_color_subjectAgg','Median_orientation_subjectAgg','Median_conjunction_subjectAgg','N_subjects_feature'});
end

function q=bh(p)
p=double(p(:)); q=NaN(size(p)); ok=isfinite(p)&p>=0&p<=1; p0=p(ok); m=numel(p0);
if m==0, return; end
[ps,ord]=sort(p0); qs=ps.*m./(1:m)';
for i=m-1:-1:1, qs(i)=min(qs(i),qs(i+1)); end
qs(qs>1)=1; tmp=NaN(m,1); tmp(ord)=qs; q(ok)=tmp;
end

function fam=family(fn)
s=lower(string(fn)); fam=repmat("unknown",size(s));
isConn=contains(s,"rie")|contains(s,"pli")|contains(s,"plv")|contains(s,"coh")|contains(s,"conn")|contains(s,"c01")|contains(s,"c02")|contains(s,"c03")|contains(s,"_c");
isRBP=contains(s,"rbp")|contains(s,"relative");
isBP=contains(s,"bp_")|contains(s,"bandpower")|contains(s,"alpha")|contains(s,"beta")|contains(s,"theta")|contains(s,"delta")|contains(s,"gamma");
isTemp=contains(s,"skew")|contains(s,"kurt")|contains(s,"hjorth")|contains(s,"line")|contains(s,"auc")|contains(s,"rms")|contains(s,"mean")|contains(s,"std")|contains(s,"var")|contains(s,"max")|contains(s,"min");
isComp=contains(s,"lz")|contains(s,"lzc")|contains(s,"entropy")|contains(s,"samp")|contains(s,"perm")|contains(s,"higuchi")|contains(s,"fractal");
isTF=contains(s,"morlet")|contains(s,"tf_")|contains(s,"wavelet")|contains(s,"tfr");
fam(isBP)="bandpower"; fam(isRBP)="relative_bandpower"; fam(isTemp)="temporal_statistical"; fam(isComp)="complexity"; fam(isTF)="time_frequency"; fam(isConn)="connectivity";
end

function labels=local_load_labels(file)
labels=containers.Map('KeyType','double','ValueType','char');
if ~exist(file,'file'), return; end
try
    M=readtable(file,'FileType','text','Delimiter','\t');
    if any(strcmp(M.Properties.VariableNames,'Number')), ncol='Number'; else, ncol=M.Properties.VariableNames{1}; end
    if any(strcmp(M.Properties.VariableNames,'labels')), lcol='labels'; else, lcol=M.Properties.VariableNames{2}; end
    nums=double(M.(ncol)); labs=string(M.(lcol));
    for i=1:numel(nums), if isfinite(nums(i)), labels(nums(i))=char(labs(i)); end, end
catch
end
end

function T=add_localization(T,labels)
[A,B,E,CD,PS]=parse_edges(string(T.Feature));
T.ChannelA=A; T.ChannelB=B; T.EdgeLabel_FIXED=E; T.ChannelsDetected_FIXED=CD; T.ParseStatus_FIXED=PS;
T.ChannelA_Label=strings(height(T),1); T.ChannelB_Label=strings(height(T),1); T.Edge_ElectrodeLabels=strings(height(T),1);
for i=1:height(T)
    if isfinite(A(i)), T.ChannelA_Label(i)=lab(labels,A(i)); end
    if isfinite(B(i)), T.ChannelB_Label(i)=lab(labels,B(i)); end
    if isfinite(A(i))&&isfinite(B(i)), T.Edge_ElectrodeLabels(i)=T.ChannelA_Label(i)+"-"+T.ChannelB_Label(i); end
end
end

function [A,B,E,CD,PS]=parse_edges(fn)
n=numel(fn); A=NaN(n,1); B=NaN(n,1); E=strings(n,1); CD=strings(n,1); PS=strings(n,1);
for i=1:n
    f=char(lower(fn(i))); chans=[];
    tok=regexp(f,'c0?(\d{1,2})[_-]c0?(\d{1,2})','tokens');
    if ~isempty(tok), chans=[str2double(tok{1}{1}) str2double(tok{1}{2})]; PS(i)="cXX_cYY"; end
    if isempty(chans)
        tok=regexp(f,'ch0?(\d{1,2})[_-]ch0?(\d{1,2})','tokens');
        if ~isempty(tok), chans=[str2double(tok{1}{1}) str2double(tok{1}{2})]; PS(i)="chXX_chYY"; end
    end
    chans=unique(chans(isfinite(chans)&chans>=1&chans<=128),'stable');
    if isempty(chans), PS(i)="no_channel_found"; continue; end
    CD(i)=strjoin("C"+compose("%02d",chans),";");
    if numel(chans)>=2
        A(i)=chans(1); B(i)=chans(2); E(i)=sprintf('C%02d-C%02d',min(A(i),B(i)),max(A(i),B(i)));
    else
        A(i)=chans(1); E(i)=sprintf('C%02d',A(i));
    end
end
end

function x=lab(labels,ch)
if isKey(labels,double(ch)), x=string(labels(double(ch))); else, x="C"+compose("%02d",ch); end
end

function save_outputs(cfg,testName,contrast,runNum,phaseName,T,mode)
sig=isfinite(T.q_FDR)&T.q_FDR<cfg.alphaFDR; Ts=T(sig,:);
if isempty(Ts), return; end
if strcmp(mode,'signrank') && any(strcmp(Ts.Properties.VariableNames,'RankBiserial_paired'))
    Ts.AbsEffect=abs(Ts.RankBiserial_paired); Ts=sortrows(Ts,{'q_FDR','AbsEffect'},{'ascend','descend'});
elseif strcmp(mode,'friedman') && any(strcmp(Ts.Properties.VariableNames,'KendallW'))
    Ts=sortrows(Ts,{'q_FDR','KendallW'},{'ascend','descend'});
else
    Ts=sortrows(Ts,'q_FDR','ascend');
end
tag=regexprep(sprintf('%s_%s_run%d_%s',testName,contrast,runNum,phaseName),'[^\w]','_');
out=fullfile(cfg.outDir,tag); if ~exist(out,'dir'), mkdir(out); end
writetable(Ts,fullfile(out,sprintf('STEP83_%s_significant_features.csv',tag)));
writetable(channel_summary(Ts,mode),fullfile(out,sprintf('STEP83_%s_channel_summary.csv',tag)));
writetable(edge_summary(Ts,mode),fullfile(out,sprintf('STEP83_%s_edge_summary.csv',tag)));
end

function C=channel_summary(T,mode)
channels=[]; q=[]; p=[]; eff=[];
for i=1:height(T)
    chs=[]; if isfinite(T.ChannelA(i)), chs(end+1)=T.ChannelA(i); end; if isfinite(T.ChannelB(i)), chs(end+1)=T.ChannelB(i); end
    chs=unique(chs);
    for k=1:numel(chs), channels(end+1,1)=chs(k); q(end+1,1)=T.q_FDR(i); p(end+1,1)=getp(T(i,:),mode); eff(end+1,1)=abs(geteff(T(i,:),mode)); end %#ok<AGROW>
end
if isempty(channels), C=table(); return; end
[G,key]=findgroups(channels);
C=table(key,splitapply(@numel,q,G),splitapply(@(x)min(x,[],'omitnan'),p,G),splitapply(@(x)min(x,[],'omitnan'),q,G),splitapply(@(x)mean(x,'omitnan'),eff,G),splitapply(@(x)max(x,[],'omitnan'),eff,G), ...
    'VariableNames',{'Channel','N_sig_features_involving_channel','Min_p','Min_q','MeanAbsEffect','MaxAbsEffect'});
C=sortrows(C,{'N_sig_features_involving_channel','Min_q'},{'descend','ascend'});
end

function E=edge_summary(T,mode)
valid=isfinite(T.ChannelA)&isfinite(T.ChannelB)&strlength(string(T.EdgeLabel_FIXED))>0;
if ~any(valid), E=table(); return; end
T=T(valid,:); edge=string(T.EdgeLabel_FIXED); labs=string(T.Edge_ElectrodeLabels); [G,key]=findgroups(edge);
N=splitapply(@numel,T.q_FDR,G); Minp=splitapply(@(x)min(x,[],'omitnan'),getp(T,mode),G); Minq=splitapply(@(x)min(x,[],'omitnan'),T.q_FDR,G);
effect=geteff(T,mode); Me=splitapply(@(x)mean(x,'omitnan'),effect,G); Mae=splitapply(@(x)mean(abs(x),'omitnan'),effect,G); Mxe=splitapply(@(x)max(abs(x),[],'omitnan'),effect,G);
EL=strings(size(key)); Dir=strings(size(key));
for i=1:numel(key)
    idx=G==i; u=unique(labs(idx),'stable'); u=u(strlength(u)>0); if ~isempty(u), EL(i)=u(1); end
    if strcmp(mode,'signrank')
        m=mean(effect(idx),'omitnan'); if m>0, Dir(i)="A_greater_than_B"; elseif m<0, Dir(i)="B_greater_than_A"; else, Dir(i)="mixed_or_zero"; end
    else
        Dir(i)="three_class_no_binary_direction";
    end
end
E=table(key,EL,N,Minp,Minq,Me,Mae,Mxe,Dir,'VariableNames',{'Edge','Edge_ElectrodeLabels','N_sig_features_involving_edge','Min_p','Min_q','MeanEffect','MeanAbsEffect','MaxAbsEffect','Direction'});
E=sortrows(E,{'N_sig_features_involving_edge','Min_q'},{'descend','ascend'});
end

function p=getp(T,mode)
if strcmp(mode,'signrank'), p=T.p_signrank; else, p=T.p_friedman; end
end
function e=geteff(T,mode)
if strcmp(mode,'signrank'), e=T.RankBiserial_paired; else, e=T.KendallW; end
end
