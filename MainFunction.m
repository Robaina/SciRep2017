function Val=MainFunction(GEM,GeneNames,GeneValues,Metabolites)

    %This function generates the main computational results presented in this
    %study.
    %Arguments:
    %GEM: a genome-scale model in a COBRA-like structure
    %GeneNames: a cell array of strings containing the Arabidopsis gene names
    %in the same order as the data values in GeneValues
    %GeneValues: a matrix containing the columns of gene values for the GC and
    %M cells (three first columns correspond to M cells and the following to
    %GC).
    %Metabolites: a cell with strings containing the name of the
    %metabolites from which the total pool is calculated

    GCoptions.L_min=0;
    GCoptions.L_int=0;
    GCoptions.L_max=0;
    GCoptions.samplesize=1000;
    GCoptions.Nsamples=20;
    GCoptions.tasks=[454,448:450]; %Include constraints on minimum biomass production and maintenance
    %GCoptions.tasks=[454,11,13,14,448:450];GCoptions.Carb2Oxy='T'; %Include constraint on the ratio of carboxylation to oxygenation and CBC reactions
    GCoptions.time_limit=60;
    GCoptions.OutFlag=0;
    [FVArange(:,1),FVArange(:,2)]=FVA(GEM);
    GCoptions.FVArange=FVArange;
    Moptions=GCoptions;

    %Get expression data
    disp('Processing gene expression data...')
    disp('')
    MapGeneValues=zeros(length(GEM.rxns),6);
    for i=1:6,
      MapGeneValues(:,i)=AraCOREgene2rxn(GeneNames,GeneValues(:,i),GEM.genes,GEM.grRules,0);
    end
    GCdata=zeros(length(GEM.rxns),2);
    Mdata=zeros(length(GEM.rxns),2);
    MapGeneValues=MapGeneValues/max(max(MapGeneValues));
    for i=1:length(GEM.rxns),
        GCdata(i,1)=mean(MapGeneValues(i,4:6));
        GCdata(i,2)=std(MapGeneValues(i,4:6));
        Mdata(i,1)=mean(MapGeneValues(i,1:3));
        Mdata(i,2)=std(MapGeneValues(i,1:3));
    end
    GCoptions.Dstd=GCdata(:,2);
    Moptions.Dstd=Mdata(:,2);
    Val.GCdata=GCdata;Val.Mdata=Mdata;

    %Get RegrEx solution and AO space
    disp('Getting RegrExLAD solution and sampling AO space...')
    disp('')
    Val.SolGC=RegrExLAD(GEM,GCdata(:,1),GCoptions); 
    Val.SolM=RegrExLAD(GEM,Mdata(:,1),Moptions);
    Val.SolGCFVA=RegrExFVA(GEM,GCdata(:,1),Val.SolGC,GCoptions);
    Val.SolMFVA=RegrExFVA(GEM,Mdata(:,1),Val.SolM,Moptions);
    Val.SolGCAO=RegrExAOS(GEM,GCdata(:,1),Val.SolGC,GCoptions);
    Val.SolMAO=RegrExAOS(GEM,Mdata(:,1),Val.SolM,Moptions);
    
    %Get flux-sums of sucrose thought the futile cycle
    SucRxns=find(GEM.S(patternfind(GEM.metNames,'sucrose[c]','T'),:)~=0);
    TotSucFluxSum=zeros(size(Val.SolGCAO.Vsample,2),2);
    for i=1:size(Val.SolGCAO.Vsample,2),
        TotSucFluxSum(i,1)=sum(abs(Val.SolGCAO.Vsample(SucRxns,i)));
        TotSucFluxSum(i,2)=sum(abs(Val.SolMAO.Vsample(SucRxns,i))); 
        FutileSucFluxSum(i,1)=sum(abs(Val.SolGCAO.Vsample([44,45],i)));
        FutileSucFluxSum(i,2)=sum(abs(Val.SolMAO.Vsample([44,45],i)));
    end

    Val.meanTotSucFluxSum=mean(TotSucFluxSum);
    Val.meanFutileSucFluxSum=mean(FutileSucFluxSum);

    %Compare alternative optimal flux values of GC and M cell types
    GCMlist=getGCMlist(GEM,Val.SolGCAO.Vsample,Val.SolMAO.Vsample,Metabolites,Val.SolGCFVA,Val.SolMFVA,[]);Val.GCMlist=GCMlist;
     for i=2:size(GCMlist.FluxTable,1),  
         GCMlist.FluxTable{i,2}=[GCMlist.FluxTable{i,2},' ','(',num2str(GCMlist.FluxTable{i,1}),')'];
     end
     
     %Get Tables
    Figure1Rxns=[1:3,5,54,115,339,113,374,510,413,152,179,4,420,9,328,346:348,498:500,473,117,35,36,478,39:40,327,491];
    Tab3Rxns=[6,85,7:18,456:458,22:24,26,28:30,32];
    Val.Table1=GCMlist.FluxTable([1,1+Figure1Rxns],[1:3,5,7,8:10,12,11,13:15]);
    Val.Table2=GCMlist.PoolMetabolites;
    Val.Table3=GCMlist.FluxTable([1,1+Tab3Rxns],[1:3,5,7,8:10,12,11,13:15]);
    
    %**********************************************************************
    
    %This function evaluates the differences in Mann-Whitney test results
    %between the model predictions with and without the additional
    %constraints in the CBC. Requires to compute the solution structures
    %RESULTS and RESULTSCBC separately using MainFunction
    
    %CBCconstraintsDifferences=evalConstraintEffects(RESULTS,RESULTSCBC);

end