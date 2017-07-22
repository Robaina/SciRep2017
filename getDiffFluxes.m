function Sol=getDiffFluxes(GEM,GCAO,MAO)

NRxns=size(GEM.S,2);
NPoints=size(GCAO,2);
%Evaluate (ranksum test) each reaction in AraCOREred
for i=1:NRxns,
    if ~GEM.rev(i),
    MedFluxGC(i,1)=mean(GCAO(i,:));
    MedFluxGC(i,2)=median(GCAO(i,:));
    MedFluxM(i,1)=mean(MAO(i,:));
    MedFluxM(i,2)=median(MAO(i,:));
    
    [TestGCM(i,1),TestGCM(i,2)]=ranksum(GCAO(i,:),MAO(i,:),'tail','right');
    [TestMGC(i,1),TestMGC(i,2)]=ranksum(MAO(i,:),GCAO(i,:),'tail','right');
    [TestGCeqM(i,1),TestGCeqM(i,2)]=ranksum(GCAO(i,:),MAO(i,:));
    
    elseif GEM.rev(i),
        GCAOfor=GCAO(i,GCAO(i,:)>0);GCAOfor=[GCAOfor,zeros(1,NPoints-length(GCAOfor))];
        GCAOrev=abs(GCAO(i,GCAO(i,:)<0));GCAOrev=[GCAOrev,zeros(1,NPoints-length(GCAOrev))];
        MAOfor=MAO(i,MAO(i,:)>0);MAOfor=[MAOfor,zeros(1,NPoints-length(MAOfor))];
        MAOrev=abs(MAO(i,MAO(i,:)<0));MAOrev=[MAOrev,zeros(1,NPoints-length(MAOrev))];        
        
        MedFluxGC(i,1)=mean(GCAOfor);
        MedFluxGC(i,2)=median(GCAOfor);
        MedFluxGC(i,3)=mean(GCAOrev);
        MedFluxGC(i,4)=median(GCAOrev);
        MedFluxM(i,1)=mean(MAOfor);
        MedFluxM(i,2)=median(MAOfor);
        MedFluxM(i,3)=mean(MAOrev);
        MedFluxM(i,4)=median(MAOrev);
                
        [TestGCM(i,1),TestGCM(i,2)]=ranksum(GCAOfor,MAOfor,'tail','right');
        [TestGCM(i,3),TestGCM(i,4)]=ranksum(GCAOrev,MAOrev,'tail','right');
        [TestMGC(i,1),TestMGC(i,2)]=ranksum(MAOfor,GCAOfor,'tail','right');
        [TestMGC(i,3),TestMGC(i,4)]=ranksum(MAOrev,GCAOrev,'tail','right');
        [TestGCeqM(i,1),TestGCeqM(i,2)]=ranksum(GCAOfor,MAOfor);
        [TestGCeqM(i,3),TestGCeqM(i,4)]=ranksum(GCAOrev,MAOrev);
    end
            
end
TestGCM(TestGCM<1e-8)=0;
TestMGC(TestMGC<1e-8)=0;
TestGCeqM(TestGCeqM<1e-8)=0;

Sol.MedFluxGC=MedFluxGC;
Sol.MedFluxM=MedFluxM;
Sol.TestGCM=TestGCM;
Sol.TestMGC=TestMGC;
Sol.TestGCeqM=TestGCeqM;


end
    
