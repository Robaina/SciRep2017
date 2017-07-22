function Sol=getDiffMetabolites(GEM,GCAO,MAO,PoolMetabolites)

if nargin<4,
    PoolMetabolites=[];
end
NPoints=size(GCAO,2);
NMets=size(GEM.S,1);

%Evaluate each metabolite in AraCOREred (different comparments)
for i=1:NMets,
    Rxns=[];
    Rxns=[Rxns,find(abs(GEM.S(i,:))>0)];
    for j=1:NPoints,
      FluxSumGC(i,j)=sum(abs(GCAO(Rxns,j)));
      FluxSumM(i,j)=sum(abs(MAO(Rxns,j)));
    end
    MedFluxSumGC(i,1)=mean(FluxSumGC(i,:));
    MedFluxSumGC(i,2)=median(FluxSumGC(i,:));
    MedFluxSumM(i,1)=mean(FluxSumM(i,:));
    MedFluxSumM(i,2)=median(FluxSumM(i,:));
    
    %Rank-Sum test on Flux-Sum distributions
    [TestGCM(i,1),TestGCM(i,2)]=ranksum(FluxSumGC(i,:),FluxSumM(i,:),'tail','right');
    [TestMGC(i,1),TestMGC(i,2)]=ranksum(FluxSumM(i,:),FluxSumGC(i,:),'tail','right');
    [TestGCeqM(i,1),TestGCeqM(i,2)]=ranksum(FluxSumGC(i,:),FluxSumM(i,:));
end
TestGCM(TestGCM<1e-8)=0;
TestMGC(TestMGC<1e-8)=0;
TestGCeqM(TestGCeqM<1e-8)=0;

%Evalute the total pool of each metabolite of interest (PoolMetabolites)
if ~isempty(PoolMetabolites),
    for i=1:length(PoolMetabolites),
        Rxns=[];
        CompMets{i}=patternfind(GEM.metNames,[PoolMetabolites{i},'['],'F');
        for k=1:length(CompMets{i}),
            Rxns=[Rxns,find(abs(GEM.S(CompMets{i}(k),:))>0)];
        end
        Rxns=unique(Rxns);
        for j=1:size(GCAO,2),
            PoolFluxSumGC(i,j)=sum(abs(GCAO(Rxns,j)));
            PoolFluxSumM(i,j)=sum(abs(MAO(Rxns,j)));
        end
        PoolMedFluxSumGC(i,1)=mean(PoolFluxSumGC(i,:));
        PoolMedFluxSumGC(i,2)=median(PoolFluxSumGC(i,:));
        PoolMedFluxSumM(i,1)=mean(PoolFluxSumM(i,:));
        PoolMedFluxSumM(i,2)=median(PoolFluxSumM(i,:));

        %Rank-Sum test on Flux-Sum distributions
        [PoolTestGCM(i,1),PoolTestGCM(i,2)]=ranksum(PoolFluxSumGC(i,:),PoolFluxSumM(i,:),'tail','right');
        [PoolTestMGC(i,1),PoolTestMGC(i,2)]=ranksum(PoolFluxSumM(i,:),PoolFluxSumGC(i,:),'tail','right');
        [PoolTestGCeqM(i,1),PoolTestGCeqM(i,2)]=ranksum(PoolFluxSumGC(i,:),PoolFluxSumM(i,:));
    end
    PoolTestGCM(PoolTestGCM<1e-8)=0;
    PoolTestMGC(PoolTestMGC<1e-8)=0;
    PoolTestGCeqM(PoolTestGCeqM<1e-8)=0;
    Sol.PoolMedFluxSumGC=PoolMedFluxSumGC;
    Sol.PoolMedFluxSumM=PoolMedFluxSumM;
    Sol.PoolTestGCM=PoolTestGCM;
    Sol.PoolTestMGC=PoolTestMGC;
    Sol.PoolTestGCeqM=PoolTestGCeqM;
    Sol.PoolFluxSumGC=PoolFluxSumGC;
    Sol.PoolFluxSumM=PoolFluxSumM;
    Sol.CompMets=CompMets;
end

Sol.MedFluxSumGC=MedFluxSumGC;
Sol.MedFluxSumM=MedFluxSumM;
Sol.TestGCM=TestGCM;
Sol.TestMGC=TestMGC;
Sol.TestGCeqM=TestGCeqM;
Sol.FluxSumGC=FluxSumGC;
Sol.FluxSumM=FluxSumM;

end
    
