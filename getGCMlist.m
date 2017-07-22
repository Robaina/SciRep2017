function Val=getGCMlist(GEM,GCAO,MAO,Metabolites,SolGCFVA,SolMFVA,PlotMetabolites,Directory)

%Arguments:
%GEM: genome-scale model in a COBRA-like structure
%GCAO,MAO: the matrix containing the sampled alternative optimal flux
%distributions for GC and M cells
%PoolMetabolites: a cell with strings containing the name of the
%metabolites from which the total pool is calculated
%PlotMetabolites: either 'T', which plots the flux-sum values of the
%the metabolites in PoolMetabolites, or 'F' which does not plot the
%flux-sum values.
%Directory: a string containing the directory to which the previous plots
%are saved, if empty then they are saved to current directory.


if nargin<5 || isempty(PlotMetabolites),
   PlotMetabolites='F';
end
if nargin<6,
    Directory=[];
end

RevRxns=find(GEM.rev==1);
Rxns=size(GEM.S,2);

for i=1:length(GEM.metNames),
    GEM.metNames{i}=strrep(GEM.metNames{i},'[c]','[cytoplasm]');
    GEM.metNames{i}=strrep(GEM.metNames{i},'[p]','[peroxisome]');
    GEM.metNames{i}=strrep(GEM.metNames{i},'[h]','[chloroplast]');
    GEM.metNames{i}=strrep(GEM.metNames{i},'[m]','[mitochondrion]');
    GEM.metNames{i}=strrep(GEM.metNames{i},'[i]','[intercellular]');
end

%Compare reaction flux values over the AO space
F=getDiffFluxes(GEM,GCAO,MAO);

%Compare metabolite flux-sums over the AO space
FS=getDiffMetabolites(GEM,GCAO,MAO,Metabolites);
FSrepGC=FS.MedFluxSumGC(:,1:2);
FSrepM=FS.MedFluxSumM(:,1:2);
FSrepGCM=FSrepGC(:,1)./FSrepM(:,1);FSrepGCM(isnan(FSrepGCM))=0;

%Treat forward and reverse direction as different reactions
FrepGC=F.MedFluxGC(:,1:2);FrepGC=[FrepGC;F.MedFluxGC(RevRxns,3:4)];
FrepM=F.MedFluxM(:,1:2);FrepM=[FrepM;F.MedFluxM(RevRxns,3:4)];
FrepGCM=FrepGC(:,1)./FrepM(:,1);FrepGCM(isnan(FrepGCM))=0;
for i=1:length(RevRxns),
    GEM.rxnNames=[GEM.rxnNames;GEM.rxnNames{RevRxns(i)},' ','[Backward]'];
    GEM.rxnNames{RevRxns(i)}=[GEM.rxnNames{RevRxns(i)},' ','[Forward]'];
end
GEM.subSystems=[GEM.subSystems;GEM.subSystems(RevRxns)];
GEM.S=[GEM.S,-GEM.S(:,RevRxns)];GEM.rev=[GEM.rev;ones(length(RevRxns),1)];
TestGCM=[F.TestGCM(:,1:2);F.TestGCM(RevRxns,3:4)];
TestMGC=[F.TestMGC(:,1:2);F.TestMGC(RevRxns,3:4)];
TestGCeqM=[F.TestGCeqM(:,1:2);F.TestGCeqM(RevRxns,3:4)];
GEM.rxnMechanisms=[GEM.rxnMechanisms;GEM.rxnMechanisms(RevRxns)];
GCrange=SolGCFVA.RangeSplit(:,2);Mrange=SolMFVA.RangeSplit(:,2);GCMRangeDiff=SolGCFVA.RangeSplit(:,2)-SolMFVA.RangeSplit(:,2);

%Restrict the number of decimal places to display
for i=1:length(FrepGC(:,1)),
    
    if (FrepGC(i,1)<1e-3 || FrepGC(i,1)>1e3) && FrepGC(i,1)~=0,
        A{i,1}=num2str(FrepGC(i,1),'%.3e');
    elseif (FrepGC(i,1)>=1e-3 || FrepGC(i,1)<=1e3) && FrepGC(i,1)~=0,
        A{i,1}=num2str(FrepGC(i,1),'%.3f');
    elseif FrepGC(i,1)==0,
        A{i,1}='0';
    end
    if (FrepGC(i,2)<1e-3 || FrepGC(i,2)>1e3) && FrepGC(i,2)~=0,
        A{i,2}=num2str(FrepGC(i,2),'%.3e');
    elseif (FrepGC(i,2)>=1e-3 || FrepGC(i,2)<=1e3) && FrepGC(i,2)~=0,
        A{i,2}=num2str(FrepGC(i,2),'%.3f');
    elseif FrepGC(i,2)==0,
        A{i,2}='0';
    end
    if (FrepM(i,1)<1e-3 || FrepM(i,1)>1e3) && FrepM(i,1)~=0,
        A{i,3}=num2str(FrepM(i,1),'%.3e');
    elseif (FrepM(i,1)>=1e-3 || FrepM(i,1)<=1e3) && FrepM(i,1)~=0,
        A{i,3}=num2str(FrepM(i,1),'%.3f');
    elseif FrepM(i,1)==0,
        A{i,3}='0';
    end
    if (FrepM(i,2)<1e-3 || FrepM(i,2)>1e3) && FrepM(i,2)~=0,
        A{i,4}=num2str(FrepM(i,2),'%.3e');
    elseif (FrepM(i,2)>=1e-3 || FrepM(i,2)<=1e3) && FrepM(i,2)~=0,
        A{i,4}=num2str(FrepM(i,2),'%.3f');
    elseif FrepM(i,2)==0,
        A{i,4}='0';    
    end
    if (FrepGCM(i,1)<1e-3 || FrepGCM(i,1)>1e3) && FrepGCM(i,1)~=0,
        A{i,5}=num2str(FrepGCM(i,1),'%.3e');
    elseif (FrepGCM(i,1)>=1e-3 || FrepGCM(i,1)<=1e3) && FrepGCM(i,1)~=0,
        A{i,5}=num2str(FrepGCM(i,1),'%.3f');
    elseif FrepGCM(i,1)==0,
        A{i,5}='0';    
    end
    if (TestGCeqM(i,1)<1e-3 || TestGCeqM(i,1)>1e3) && TestGCeqM(i,1)~=0,
        A{i,6}=num2str(TestGCeqM(i,1),'%.3e');
    elseif (TestGCeqM(i,1)>=1e-3 || TestGCeqM(i,1)<=1e3) && TestGCeqM(i,1)~=0,
        A{i,6}=num2str(TestGCeqM(i,1),'%.3f');
    elseif TestGCeqM(i,1)==0,
        A{i,6}='0';    
    end
    if (TestGCM(i,1)<1e-3 || TestGCM(i,1)>1e3) && TestGCM(i,1)~=0,
        A{i,7}=num2str(TestGCM(i,1),'%.3e');
    elseif (TestGCM(i,1)>=1e-3 || TestGCM(i,1)<=1e3) && TestGCM(i,1)~=0,
        A{i,7}=num2str(TestGCM(i,1),'%.3f');
    elseif TestGCM(i,1)==0,
        A{i,7}='0';        
    end
    if (TestMGC(i,1)<1e-3 || TestMGC(i,1)>1e3) && TestMGC(i,1)~=0,
        A{i,8}=num2str(TestMGC(i,1),'%.3e');
    elseif (TestMGC(i,1)>=1e-3 || TestMGC(i,1)<=1e3) && TestMGC(i,1)~=0,
        A{i,8}=num2str(TestMGC(i,1),'%.3f');
    elseif TestMGC(i,1)==0,
        A{i,8}='0';        
    end
    
    if (GCrange(i)<1e-3 || GCrange(i)>1e3) && GCrange(i)~=0,
        D{i,1}=num2str(GCrange(i),'%.3e');
    elseif (GCrange(i)>=1e-3 || GCrange(i)<=1e3) && GCrange(i)~=0,
        D{i,1}=num2str(GCrange(i),'%.3f');
    elseif GCrange(i)==0,
        D{i,1}='0';
    end
    
    if (Mrange(i)<1e-3 || Mrange(i)>1e3) && Mrange(i)~=0,
        D{i,2}=num2str(Mrange(i),'%.3e');
    elseif (Mrange(i)>=1e-3 || Mrange(i)<=1e3) && Mrange(i)~=0,
        D{i,2}=num2str(Mrange(i),'%.3f');
    elseif Mrange(i)==0,
        D{i,2}='0';
    end
    
    if (GCMRangeDiff(i)<1e-3 || GCMRangeDiff(i)>1e3) && GCMRangeDiff(i)~=0,
        D{i,3}=num2str(GCMRangeDiff(i),'%.3e');
    elseif (GCMRangeDiff(i)>=1e-3 || GCMRangeDiff(i)<=1e3) && GCMRangeDiff(i)~=0,
        D{i,3}=num2str(GCMRangeDiff(i),'%.3f');
    elseif GCMRangeDiff(i)==0,
        D{i,3}='0';
    end
end
for i=1:length(FSrepGC(:,1)),
    if (FSrepGC(i,1)<1e-3 || FSrepGC(i,1)>1e3) && FSrepGC(i,1)~=0,
        B{i,1}=num2str(FSrepGC(i,1),'%.3e');
    elseif (FSrepGC(i,1)>=1e-3 || FSrepGC(i,1)<=1e3) && FSrepGC(i,1)~=0,
        B{i,1}=num2str(FSrepGC(i,1),'%.3f');
    elseif FSrepGC(i,1)==0,
        B{i,1}='0';
    end
    if (FSrepGC(i,2)<1e-3 || FSrepGC(i,2)>1e3) && FSrepGC(i,2)~=0,
        B{i,2}=num2str(FSrepGC(i,2),'%.3e');
    elseif (FSrepGC(i,2)>=1e-3 || FSrepGC(i,2)<=1e3) && FSrepGC(i,2)~=0,
        B{i,2}=num2str(FSrepGC(i,2),'%.3f');
    elseif FSrepGC(i,2)==0,
        B{i,2}='0';
    end
    if (FSrepM(i,1)<1e-3 || FSrepM(i,1)>1e3) && FSrepM(i,1)~=0,
        B{i,3}=num2str(FSrepM(i,1),'%.3e');
    elseif (FSrepM(i,1)>=1e-3 || FSrepM(i,1)<=1e3) && FSrepM(i,1)~=0,
        B{i,3}=num2str(FSrepM(i,1),'%.3f');
    elseif FSrepM(i,1)==0,
        B{i,3}='0';
    end
    if (FSrepM(i,2)<1e-3 || FSrepM(i,2)>1e3) && FSrepM(i,2)~=0,
        B{i,4}=num2str(FSrepM(i,2),'%.3e');
    elseif (FSrepM(i,2)>=1e-3 || FSrepM(i,2)<=1e3) && FSrepM(i,2)~=0,
        B{i,4}=num2str(FSrepM(i,2),'%.3f');
    elseif FSrepM(i,2)==0,
        B{i,4}='0';    
    end
    if (FSrepGCM(i,1)<1e-3 || FSrepGCM(i,1)>1e3) && FSrepGCM(i,1)~=0,
        B{i,5}=num2str(FSrepGCM(i,1),'%.3e');
    elseif (FSrepGCM(i,1)>=1e-3 || FSrepGCM(i,1)<=1e3) && FSrepGCM(i,1)~=0,
        B{i,5}=num2str(FSrepGCM(i,1),'%.3f');
    elseif FSrepGCM(i,1)==0,
        B{i,5}='0';    
    end
    if (FS.TestGCeqM(i,1)<1e-3 || FS.TestGCeqM(i,1)>1e3) && FS.TestGCeqM(i,1)~=0,
        B{i,6}=num2str(FS.TestGCeqM(i,1),'%.3e');
    elseif (FS.TestGCeqM(i,1)>=1e-3 || FS.TestGCeqM(i,1)<=1e3) && FS.TestGCeqM(i,1)~=0,
        B{i,6}=num2str(FS.TestGCeqM(i,1),'%.3f');
    elseif FS.TestGCeqM(i,1)==0,
        B{i,6}='0';    
    end
    if (FS.TestGCM(i,1)<1e-3 || FS.TestGCM(i,1)>1e3) && FS.TestGCM(i,1)~=0,
        B{i,7}=num2str(FS.TestGCM(i,1),'%.3e');
    elseif (FS.TestGCM(i,1)>=1e-3 || FS.TestGCM(i,1)<=1e3) && FS.TestGCM(i,1)~=0,
        B{i,7}=num2str(FS.TestGCM(i,1),'%.3f');
    elseif FS.TestGCM(i,1)==0,
        B{i,7}='0';        
    end
    if (FS.TestMGC(i,1)<1e-3 || FS.TestMGC(i,1)>1e3) && FS.TestMGC(i,1)~=0,
        B{i,8}=num2str(FS.TestMGC(i,1),'%.3e');
    elseif (FS.TestMGC(i,1)>=1e-3 || FS.TestMGC(i,1)<=1e3) && FS.TestMGC(i,1)~=0,
        B{i,8}=num2str(FS.TestMGC(i,1),'%.3f');
    elseif FS.TestMGC(i,1)==0,
        B{i,8}='0';        
    end
end

Val.FluxTable=[{'Reaction Idx','Reaction [sense]','Mean Flux GC','Median','Mean Flux M','Median','Mean ratio (GC/M)','MWW (ho: GC=M) p-value','MWW (ho: M > GC) p-value','MWW (ho: GC > M) p-value','Mechanism','SubSystem','VmaxG','VmaxM','VmaxG - VmaxM'};[num2cell([(1:Rxns)';RevRxns]),GEM.rxnNames,A,GEM.rxnMechanisms,GEM.subSystems,D]];
Val.FluxSumTable=[{'Metabolite','Mean FluxSum GC','Median','Mean FluxSum M','Median','Mean ratio (GC/M)','MWW (ho: GC=M) p-value','MWW (ho: M > GC) p-value','MWW (ho: GC > M) p-value'};[GEM.metNames,B]];
Val.TestGCM=TestGCM;
Val.TestMGC=TestMGC;
Val.TestGCeqM=TestGCeqM;

%get Pool Metabolites list
TotFSGC=[];TotFSMets={};
n=0;
for i=1:length(Metabolites),
        GCD=FS.MedFluxSumGC(FS.CompMets{i},1)./FS.MedFluxSumM(FS.CompMets{i},1);
        PoolGCD=FS.PoolMedFluxSumGC(i,1)/FS.PoolMedFluxSumM(i,1);
        TotFSGC=[TotFSGC;[[FS.PoolMedFluxSumGC(i,1),FS.PoolMedFluxSumM(i,1),PoolGCD,FS.PoolTestGCeqM(i,1),FS.PoolTestGCM(i,1),FS.PoolTestMGC(i,1)];[FS.MedFluxSumGC(FS.CompMets{i},1),FS.MedFluxSumM(FS.CompMets{i},1),GCD,FS.TestGCeqM(FS.CompMets{i},1),FS.TestGCM(FS.CompMets{i},1),FS.TestMGC(FS.CompMets{i},1)]]];
        TotFSMets(n+1:(n+length([Metabolites{i};GEM.metNames(FS.CompMets{i})])))=[['Total',' ',Metabolites{i}];GEM.metNames(FS.CompMets{i})];
        n=length(TotFSMets);
end
TotFSMets=TotFSMets';

%Include sucrose flux-sum in futile cycle
SucRxns=find(GEM.S(patternfind(GEM.metNames,'sucrose[c]','T'),:)~=0);
SucRxns(ismember(SucRxns,[44,45]))=[];
SucFluxSum=zeros(size(GCAO,2),2);FutileSucFluxSum=SucFluxSum;
for i=1:size(GCAO,2),
    FutileSucFluxSum(i,1)=sum(abs(GCAO([44,45],i)));
    FutileSucFluxSum(i,2)=sum(abs(MAO([44,45],i)));
    SucFluxSum(i,1)=sum(abs(GCAO(SucRxns,i)));
    SucFluxSum(i,2)=sum(abs(MAO(SucRxns,i)));
end
[TestFutSucGCM]=ranksum(FutileSucFluxSum(:,1),FutileSucFluxSum(:,2),'tail','right');
[TestFutSucMGC]=ranksum(FutileSucFluxSum(:,1),FutileSucFluxSum(:,2),'tail','left');
[TestFutSucGCeqM]=ranksum(FutileSucFluxSum(:,1),FutileSucFluxSum(:,2));
[TestSucGCM]=ranksum(SucFluxSum(:,1),SucFluxSum(:,2),'tail','right');
[TestSucMGC]=ranksum(SucFluxSum(:,1),SucFluxSum(:,2),'tail','left');
[TestSucGCeqM]=ranksum(SucFluxSum(:,1),SucFluxSum(:,2));
TotFSMets2=TotFSMets;
SucIdx=patternfind(TotFSMets(:,1),'Sucrose[cytoplasm]');
TotFSMets{SucIdx}='Futile Cycle';
TotFSGC2=TotFSGC;
TotFSGC(SucIdx,:)=[mean(FutileSucFluxSum),mean(FutileSucFluxSum(:,1))/mean(FutileSucFluxSum(:,2)),TestFutSucGCM,TestFutSucMGC,TestFutSucGCeqM];
TotFSMets=[TotFSMets(1:SucIdx);'Non Futile Cycle';TotFSMets(SucIdx+1:end)];
TotFSGC=[TotFSGC(1:SucIdx,:);[mean(SucFluxSum),mean(SucFluxSum(:,1))/mean(SucFluxSum(:,2)),TestSucGCM,TestSucMGC,TestSucGCeqM];TotFSGC(SucIdx+1:end,:)];

%Restrict number of significant digits in table
for j=1:size(TotFSGC,2),
    for i=1:size(TotFSGC,1),
        if (TotFSGC(i,j)<1e-3 || TotFSGC(i,j)>1e3) && TotFSGC(i,j)~=0,
            C{i,j}=num2str(TotFSGC(i,j),'%.3e');
        elseif (TotFSGC(i,j)>=1e-3 || TotFSGC(i,j)<=1e3) && TotFSGC(i,j)~=0,
            C{i,j}=num2str(TotFSGC(i,j),'%.3f');
        elseif TotFSGC(i,j)==0,
            C{i,j}='0';
        end
    end
end
Val.PoolMetabolites=[{'Metabolite','Mean FluxSum GC','Mean FluxSum M','Ratio (GC/M)','MWW (ho: GC = M) p-value','MWW (ho: M > GC) p-value','MWW (ho: GC > M) p-value'};[TotFSMets,C]];

if PlotMetabolites=='T',
    %Plot Flux-Sums for selected metabolites
    numbins=20;
    for i=1:length(Metabolites),
        fr(i)=figure('Color','w','Visible','off');
        [Cr,Pr]=hist(FS.PoolFluxSumGC(i,:),numbins);
        bar(Pr,Cr,1,'FaceColor','r');
        title(['mean Flux-Sum (\sigma) =',' ',num2str(FS.PoolMedFluxSumGC(i,1)),'(',num2str(FS.PoolMedFluxSumGC(i,2)),')'],'FontSize',14);
        set(gca,'fontsize',18)
        xlabel([Metabolites{i},' ','Total Flux-Sum'],'FontSize',16,'FontWeight','bold')
        ylabel('counts','FontSize',18)
        print(fr(i),[Directory,Metabolites{i},'GC'],'-dtiff')

        gr(i)=figure('Color','w','Visible','off');
        [Cr,Pr]=hist(FS.PoolFluxSumM(i,:),numbins);
        h=bar(Pr,Cr,1,'FaceColor','k');
        set(h,'edgecolor','w')
        title(['mean Flux-Sum (\sigma) =',' ',num2str(FS.PoolMedFluxSumM(i,1)),'(',num2str(FS.PoolMedFluxSumM(i,2)),')'],'FontSize',14);
        set(gca,'fontsize',18)
        xlabel([Metabolites{i},' ','Total Flux-Sum'],'FontSize',16,'FontWeight','bold')
        ylabel('counts','FontSize',18)
        print(gr(i),[Directory,Metabolites{i},'M'],'-dtiff')    

       Lmets=patternfind(TotFSMets2,Metabolites{i});

       %GC-M comparison
       XTickLabels=TotFSMets2(Lmets);
       XTickLabels{1}='Total';
       for l=2:length(Lmets),
           if strfind(XTickLabels{l},'[cytoplasm]'),
               XTickLabels{l}='c';
           elseif strfind(XTickLabels{l},'[peroxisome]'),
               XTickLabels{l}='p';
           elseif strfind(XTickLabels{l},'[chloroplast]'),
               XTickLabels{l}='h';
           elseif strfind(XTickLabels{l},'[mitochondrion]'),
               XTickLabels{l}='m';
           elseif strfind(XTickLabels{l},'[intercellular]'),
               XTickLabels{l}='i';
           end
       end 
       for k=1:size(TotFSGC2,1),
        TotFSGC2(k,1:2)=TotFSGC2(k,1:2)/max(TotFSGC2(k,1:2));
       end
       for l=1:length(Lmets),
           XTickLabels{l}=[XTickLabels{l},'(',num2str(TotFSGC2(Lmets(l),5),3),')'];
       end
       hr(i)=figure('Color','w','Visible','off');
       g=bar(TotFSGC2(Lmets,[2,1]));
       if min(TotFSGC2(Lmets,[2,1]))>0.9,
           axis([0 length(Lmets)+1 0.9 1])
       elseif strcmpi(Metabolites{i},'Sucrose'),
           axis([0 length(Lmets)+1 0 0.02])
       end
       set(g(1),'FaceColor','k');set(g(2),'FaceColor','r')
           xlabel(Metabolites{i},'FontSize',18,'FontWeight','bold');
       ylabel('Normalized Flux-Sum','FontSize',18,'FontWeight','bold');
       axisHandle = gca;
       set(gca,'fontsize',22,'XTickLabel',XTickLabels,'FontWeight','bold')
       axisHandle.XTickLabelRotation = 45;
       print(hr(i),[Directory,Metabolites{i},'GC_M'],'-dtiff') 

    end
    close all 
end

FS.MedFluxSumGC(:,1:2)=FSrepGC;
FS.MedFluxSumM(:,1:2)=FSrepM;
Val.FluxSumGC=FS.FluxSumGC;
Val.FluxSumM=FS.FluxSumM;

end