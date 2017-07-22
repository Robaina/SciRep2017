function val=evalConstraintEffects(RESULTS,RESULTSCBCconstraints,SigThres)

%This functions computes the Hamming distance betwen the binary vectors
%containing the results of the Mann-Whitney tests comparing the altnerative
%optimal flux distributions of G and M cells

if nargin<4,
    SigThres=0.05;
end
 for i=1:525,
   if isempty(RESULTS.GCMlist.FluxTable{i+1,10}),
      GCMpv(i,1)=nan;
   elseif ~isempty(RESULTS.GCMlist.FluxTable{i+1,10}),
      GCMpv(i,1)=str2num(RESULTS.GCMlist.FluxTable{i+1,10});
   end
   if isempty(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,10}),
      GCMpv(i,2)=nan;
   elseif ~isempty(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,10}),
      GCMpv(i,2)=str2num(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,10});
   end

   if isempty(RESULTS.GCMlist.FluxTable{i+1,8}),
       GCeqMpv(i,1)=nan;
   elseif ~isempty(RESULTS.GCMlist.FluxTable{i+1,8}),
      GCeqMpv(i,1)=str2num(RESULTS.GCMlist.FluxTable{i+1,8});
   end
   if isempty(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,8}),
       GCeqMpv(i,2)=nan;
   elseif ~isempty(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,8}),
      GCeqMpv(i,2)=str2num(RESULTSCBCconstraints.GCMlist.FluxTable{i+1,8});
   end
 end
 
 GCMTests=zeros(525,2);GCeqMTests=zeros(525,2);
 GCMTests(GCMpv(:,1)<SigThres,1)=1;GCMTests(GCMpv(:,2)<SigThres,2)=1;
 GCeqMTests(GCeqMpv(:,1)<SigThres,1)=1;GCeqMTests(GCeqMpv(:,2)<SigThres,2)=1;
 GCMTests(isnan(GCeqMpv(:,1)),1)=2;GCMTests(isnan(GCeqMpv(:,2)),2)=2;
 
 %Get reactions in the CBC cycle
 NonCBCrxns=setdiff(1:size(GCMTests,1),patternfind(RESULTS.GCMlist.FluxTable(:,12),'Calvin-Benson cycle')-1);
 CBCrxns=patternfind(RESULTS.GCMlist.FluxTable(:,12),'Calvin-Benson cycle')-1;
 
 %Reactions depicted in Figure 1(main results of this study)
 TableRxns=cell2mat(RESULTS.GCMlist.FluxTable(2:end,1));
 Figure1Rxns=find(ismember(TableRxns,[1:3,5,54,115,339,113,374,510,413,152,179,4,420,9,328,346:348,498:500,473,117,35,36,478,39:40,327,491]));
 StarchRxns=find(ismember(TableRxns,[22:30,32,33]));
 
 %Obtain Hamming distance
 val.TotalDist=HammingMat(GCMTests)/525;
 val.NonCBCDist=HammingMat(GCMTests(NonCBCrxns,:))/525;
 val.Fig1Dist=HammingMat(GCMTests(Figure1Rxns,:))/525;
 
 %Identify mismatching reactions
 GCMTestsFig1=GCMTests(Figure1Rxns,:);
 GCMTestsCBC=GCMTests(CBCrxns,:);
 GCMTestsStarch=GCMTests(StarchRxns,:);
 
 A=find(GCMTestsFig1(:,1)~=GCMTestsFig1(:,2));
 C=find(GCMTestsCBC(:,1)~=GCMTestsCBC(:,2));
 E=find(GCMTestsStarch(:,1)~=GCMTestsStarch(:,2));
 
 val.DiffFig1_1_2=RESULTS.GCMlist.FluxTable(Figure1Rxns(A)+1,[1,2,11,12]);
 val.DiffCBC_1_2=RESULTS.GCMlist.FluxTable(CBCrxns(C)+1,[1,2,11,12]);
 val.DiffStarch_1_2=RESULTS.GCMlist.FluxTable(StarchRxns(E)+1,[1,2,11,12]);
end