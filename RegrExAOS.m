function Val=RegrExAOS(GEM,D,RegrExSol,Options)
      
%************RegrEx Alternative Optima Sampling****************************

%Depends on Gurobi, which can be found in www.gurobi.com (free academic
%licenses available) FVA function (provided in supplementary material)
%
%Arguments:
%
% Required
%
% GEM: Metabolic model in COBRA structure. Blocked reactions should be
% removed first, for instance using the reduceModel function of the COBRA
% toolbox D: Experimental data vector (values already mapped to
% reactions)must be of same length as the number of reactions. It can be
% provided by the genetorxn function Vopt: flux vector obtained by RegrEx
%
% Optional
%
% Options must be a structure with the following fields (not all required)
%
%     Options.tasks: A vector indicating the indexes of the reactions that
%     should be active in the final model (i.e. flux value above eps)
%     
%     Options.blocked: A vector indicating the indexes of the reactions
%     that should not be active in the final model (i.e. zero flux value)
%     
%     Options.rxns_lb: a vector containing the lower bounds for each
%     reaction (default is set to 0, min value is 0, reversible reactions
%     are splitted) 
%     
%     Options.rxns_ub: a vector containing the upper bounds for each 
%     reaction (default is set to 1) 
%    
%     Options.time_limit: scalar indicating the time limit in seconds for 
%     the MIQP (default is set to 60) 
%     
%     Options.eps: scalar indicating the threshold value to consider a
%     reaction active, i.e. a reaction i with flux value vi is active if
%     vi>=eps Options.E_max: vector of length equal to the number of
%     reactions indicating the maximum allowed error value between data and
%     flux prediction, i.e. E=d-v, (default is set to 1000) 
%     
%     Options.k: diameter of the perturbation, k in [0,1], represents the 
%     fraction of the maximum and minimum absolute perturbation to V. 
%     Which corresponds to the maximum and minimum values calculated by FVA 
%     on the original GEM.
%     
%     Options.samplesize: sample size, default is 100 points
%     
%     Options.FVArange: a 2 column matrix containing the minimun and
%     maximum values for each reaction, calculated per FVA
%
%     Options.delta: deviation from the FVA interval to be taken into
%     account when sampling the flux space
%
%     Options.Carb2Oxy: If 'T' then the constraint on the ratio of
%     carboxlation to oxygenation is imposed.
%    
%Value:
%   Vsample: a matrix with sample alternative optimal flux distributions as
%   columns
%   Vsstotal: the first norm of each sampled flux distribution
%   Etotal: the total error between flux distrubion and data (weighted by
%   W), should be similar to the original ZE field in RegrExLAD
%
%**************************************************************************
%              Semidan(robaina@mpimp-golm.mpg.de), March, 2015
%**************************************************************************

%Argument evaluation

if ~exist('GEM','var'),
    error('The Genome-Scale model is missing...')
end

if ~exist('D','var'),
    error('The data vector is missing...')
end

if ~exist('RegrExSol','var'),
    error('The RegrEx solution structure is missing...')
end

if ~exist('Options','var'),
    Options=struct;
end

if ~isfield(Options,'tasks'), 
    tasks=[];
else
    tasks=Options.tasks;
end

if ~isfield(Options,'blocked'), 
    blocked=[];
else
    blocked=Options.blocked;
end

if ~isfield(Options,'rxns_lb'), 
    rxns_lb=0;
else
    rxns_lb=Options.rxns_lb;
end

if ~isfield(Options,'rxns_ub'), 
    rxns_ub=1;
else
    rxns_ub=Options.rxns_ub;
end

if ~isfield(Options,'time_limit'), 
    time_limit=60;
else
    time_limit=Options.time_limit;
end

if ~isfield(Options,'FeasTol'), 
    FeasTol=1e-9;
else
    FeasTol=Options.FeasTol;
end

if ~isfield(Options,'MaxCapacity'),
    MaxCapacity=1;
else
    MaxCapacity=Options.MaxCapacity;
end

if ~isfield(Options,'eps'), 
    eps=MaxCapacity*1e-6;
else
    eps=Options.eps;
end

if ~isfield(Options,'E_max'), 
    E_max=MaxCapacity*1e3;
else
    E_max=Options.E_max;
end
if ~isfield(Options,'samplesize') && ~isfield(Options,'Vp'),
    samplesize=1e2;
elseif isfield(Options,'Vp'),
    samplesize=size(Options.Vp,2);
elseif isfield(Options,'samplesize') && ~isfield(Options,'Vp'),
    samplesize=Options.samplesize;
end
if ~isfield(Options,'FVArange'), 
    [FVArange(:,1),FVArange(:,2)]=FVA(GEM,1,'gurobi','F');
else
    FVArange=Options.FVArange;
end

if ~isfield(Options,'OutFlag'), 
    OutFlag=0;
else
    OutFlag=Options.OutFlag;
end

if ~isfield(Options,'Nsamples'), 
    Nsamples=1;
else
    Nsamples=Options.Nsamples;
end

if ~isfield(RegrExSol,'W'),
    W=1;
end

if ~isfield(Options,'deltaZE'),
    deltaZE=1e-4;
elseif isfield(Options,'deltaZE'),
    deltaZE=Options.deltaZE;
end

tic
    %Parse Data    
    for i=1:length(D),
        if isnan(D(i)) || isinf(D(i)),
            D(i)=0;
            sprintf('Error: Data point D(%d) is NaN or Inf and has been converted to 0',i);
        end
    end
    if max(D)>MaxCapacity,
       D=D/max(D);
    end
    D=MaxCapacity*D;
    Vopt=RegrExSol.Flux;
    Vopt=round((1/eps)*Vopt)/(1/eps);
    S=GEM.S;
    SolV=zeros(size(S,2),samplesize);
    Vptotal=SolV;
    Etotal=zeros(samplesize,1);
    Vsstotal=Etotal;
    
    %Reaction partition: Irr_D, Irr_nD, Bfordat, Bfornodat, Rev_D, Rev_nD   
    Rev=find(GEM.rev==1);
    Irr=(setdiff(1:size(S,2),Rev))';
    Irr_D=Irr((D(Irr)~=0));
    Irr_nD=Irr((D(Irr)==0));
    Rev_D=Rev((D(Rev)~=0));
    Rev_nD=Rev((D(Rev)==0));
    NIrr_D=length(Irr_D);
    NIrr_nD=length(Irr_nD);
    NRev_D=length(Rev_D);
    NRev_nD=length(Rev_nD);
    NIrr=length(Irr);NRev=length(Rev);
    Rxns=size(S,2);Mets=size(S,1);
    Rxns_Or=[Irr_D;Irr_nD;Rev_D;Rev_nD];
    
    D=[D;D(Rev)];
    DIrr=D(Irr_D);
    DRev=D(Rev_D);
    DFor=DRev;
    ZV=RegrExSol.ZV;
    ZE=RegrExSol.ZE;
    
    %***************************Main Loop Starts***************************
    %**********************************************************************
    
    wbar = waitbar(0,'Sampling Alternative Optima...');
    counter=1;
    SolVT=[];
    EtotalT=[];
    VsstotalT=[];
    Ntotal=samplesize*Nsamples;
   for l=1:Nsamples,
       j=1;
    while j<=samplesize,
        waitbar(counter/Ntotal)
        
        %Generate random perturbation to Vopt and transform S matrix
        %accordingly
        
        clear flipped 
        St=S;
        if ~isfield(Options,'Vp'), 
            epsmin=-Vopt+FVArange(:,1);
            epsmax=FVArange(:,2)-Vopt;
            Vp=Vopt+((epsmax-epsmin).*rand(length(Vopt),1)+epsmin);
        elseif isfield(Options,'Vp'),
            Vp=Options.Vp(:,j);
        end
        p=1;
        for i=1:length(Vp),
            if Vp(i)<-eps,
                Vp(i)=abs(Vp(i));
                flipped(p)=i;
                p=p+1;
                St(:,i)=-St(:,i);
            end
        end
        
        %Stoichiometric Matrix reorganization
        
        Sam=[St(:,Irr_D),St(:,Irr_nD),St(:,Rev_D),St(:,Rev_nD)];  
        
        %Bounds
        
        %Construction of sense, c and b vectors, definition of tasks and 
        %blocked reactions
    
        Vss1min=zeros(Rxns,1);
        Vss1max=MaxCapacity*ones(Rxns,1);

        if length(rxns_lb)~=1,
          for h=1:size(S,2),
              Vss1min(ismember(Rxns_Or,h))=rxns_lb(h);
          end
        end
        if length(rxns_ub)~=1,
            for h=1:size(S,2),
                Vss1max(ismember(Rxns_Or,h))=rxns_ub(h);
            end
        end
        if ~isempty(tasks),
            for i=1:length(tasks),
                if tasks(i)==454,
                    Vss1min((ismember(Rxns_Or,tasks(i))))=eps;
                elseif tasks(i)~=454,
                    Vss1min((ismember(Rxns_Or,tasks(i))))=0.001;
                end
            end
        end
        if ~isempty(blocked),
           Vss1max((ismember(Rxns_Or,blocked)))=eps;
        end

        Vss2min=zeros(Rxns,1);
        Vss2min(NIrr+1:Rxns)=Vss1min(NIrr+1:Rxns);
        Vss2max=ones(Rxns,1);
        Vss2max(1:NIrr)=0;
        
        VpIrr=Vp([Irr_D;Irr_nD]);
        VpRev=Vp([Rev_D;Rev_nD]);

        vecsense=[repmat('=',Mets,1);repmat('=',NIrr+NRev,1);repmat('<',2*NRev,1);repmat('=',NRev+NIrr_D+2*NRev_D,1);'<';'>';repmat('=',1,1);repmat('>',2*NRev,1)];
        b=[zeros(Mets,1);VpIrr;zeros(NRev,1);zeros(NRev,1)+eps;Vss2max(NIrr+1:Rxns)+eps;VpRev;DIrr;zeros(NRev_D,1);DRev;ZE+deltaZE*ZE;ZE-deltaZE*ZE;ZV;zeros(NRev,1)-eps;Vss2min([Rev_D;Rev_nD])-eps];
        lb=[Vss1min(1:NIrr);zeros(NRev,1);Vss2min;zeros(2*Rxns+2*NRev+2*NIrr_D+4*NRev_D,1)];
        ub=[Vss1max;zeros(NIrr,1);MaxCapacity*ones(NRev,1);E_max*ones(2*Rxns+NRev+2*NIrr_D+4*NRev_D,1);ones(NRev,1)];

        %Integrates weighting based on standard deviation of each data point
   
        if isfield(RegrExSol,'W'),
            W=RegrExSol.W;            
            W=[W(Irr_D);W(Rev_D);W(Irr_D);W(Rev_D);W(Rev_D);W(Rev_D)]';
        elseif ~isfield(RegrExSol,'W'),
            W=1;
        end
        
        %Construction of A matrix:  

        %Variables: Vss1(Rxns),Vss2(Rxns),delta1plus(Rxns),delta1minus(Rxns),delta2(NRev),E1plus(NIrr_D+NRev_D),E1minus(NIrr_D+NRev_D),E2plus(NRev_D),E2minus(NRev_D), Y(NRev)
        A1=[Sam,-Sam,zeros(Mets,2*Rxns+2*NRev+2*NIrr_D+4*NRev_D)]; %S*(Vss1-Vss2)=0    
        A4=[diag(ones(NIrr,1)),zeros(NIrr,NRev+Rxns),-diag(ones(NIrr,1)),zeros(NIrr,NRev),diag(ones(NIrr,1)),zeros(NIrr,2*NRev+2*NIrr_D+4*NRev_D+NRev)]; % Vss1=V'+dplus-dminus (IrrRxns) 
        A5=[zeros(NRev,NIrr),diag(ones(NRev,1)),zeros(NRev,Rxns+NIrr),-diag(ones(NRev,1)),zeros(NRev,NIrr),diag(ones(NRev,1)),zeros(NRev,NRev+2*NIrr_D+2*NRev_D+2*NRev_D),-diag(VpRev.*ones(NRev,1))];% Vss1=y*V'+dplus-dminus (RevRxns)    
        A6=[zeros(NRev,NIrr),diag(ones(NRev,1)),zeros(NRev,3*Rxns+NRev+2*NIrr_D+2*NRev_D+2*NRev_D),-diag(Vss1max([Rev_D;Rev_nD]).*ones(NRev,1))];% Vss1 XOR Vss2 (RevRxns) 
        A7=[zeros(NRev,Rxns+NIrr),diag(ones(NRev,1)),zeros(NRev,2*Rxns+NRev+2*NIrr_D+2*NRev_D+2*NRev_D),diag(Vss2max(NIrr+1:Rxns).*ones(NRev,1))];% Vss1 XOR Vss2 (RevRxns)
        A8=[zeros(NRev,Rxns+NIrr),-diag(ones(NRev,1)),zeros(NRev,2*Rxns),diag(ones(NRev,1)),zeros(NRev,2*NIrr_D+2*NRev_D+2*NRev_D),diag(VpRev.*ones(NRev,1))];% K2=Vpf+Vss2
        A9=[diag(ones(NIrr_D,1)),zeros(NIrr_D,NIrr_nD+NRev),zeros(NIrr_D,3*Rxns+NRev),diag(ones(NIrr_D,1)),zeros(NIrr_D,NRev_D),-diag(ones(NIrr_D,1)),zeros(NIrr_D,NRev_D+2*NRev_D+NRev)];% Vss1+E1plus-E1minus=D (IrrRxns)
        A10=[zeros(NRev_D,NIrr),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD+3*Rxns+NRev),zeros(NRev_D,NIrr_D),diag(ones(NRev_D,1)),zeros(NRev_D,NIrr_D),-diag(ones(NRev_D,1)),zeros(NRev_D,2*NRev_D),-diag(DFor),zeros(NRev_D,NRev_nD)];% Vss1+E1plus-E1minus-y(DRev)=0
        A11=[zeros(NRev_D,Rxns+NIrr),diag(ones(NRev_D,1)),zeros(NRev_D,NRev_nD+2*Rxns+NRev+2*NIrr_D+2*NRev_D),diag(ones(NRev_D,1)),-diag(ones(NRev_D,1)),diag(DRev),zeros(NRev_D,NRev_nD)];% Vss2+E2plus-E2minus+y(DRev)=DRev 
        A16a=[zeros(1,4*Rxns+NRev),W.*ones(1,2*NIrr_D+4*NRev_D),zeros(1,NRev)]; % ||E1plus+E1minus+E2plus+E2minus||_1 < ZE + eps
        A16b=[zeros(1,4*Rxns+NRev),W.*ones(1,2*NIrr_D+4*NRev_D),zeros(1,NRev)]; % ||E1plus+E1minus+E2plus+E2minus||_1 > ZE - eps     
        A17=[ones(1,2*Rxns),zeros(1,2*Rxns+2*NRev+2*NIrr_D+4*NRev_D)]; % ||Vss1+Vss2||_1 = ZV
        A6b=[zeros(NRev,NIrr),diag(ones(NRev,1)),zeros(NRev,3*Rxns+NRev+2*NIrr_D+2*NRev_D+2*NRev_D),-diag(Vss1min([Rev_D;Rev_nD]).*ones(NRev,1))];% Vss1Rev>=Vss1Revmin
        A7b=[zeros(NRev,Rxns+NIrr),diag(ones(NRev,1)),zeros(NRev,2*Rxns+NRev+2*NIrr_D+2*NRev_D+2*NRev_D),diag(Vss2min([Rev_D;Rev_nD]).*ones(NRev,1))];% Vss2>=Vss2min
        
        if isfield(Options,'Carb2Oxy') && Options.Carb2Oxy=='T',
            %Include carboxylation to oxygenation ratio: 1.5 <= V6/V85 <= 4
            Rat1=zeros(1,Rxns);Rat1(ismember(Rxns_Or,6))=1;Rat1(ismember(Rxns_Or,85))=-4;
            Rat2=zeros(1,Rxns);Rat2(ismember(Rxns_Or,6))=1;Rat2(ismember(Rxns_Or,85))=-1.5;
            A18=[Rat1,zeros(1,3*Rxns+2*NRev+2*NIrr_D+4*NRev_D)];
            A19=[Rat2,zeros(1,3*Rxns+2*NRev+2*NIrr_D+4*NRev_D)];      
            Amat=[A1;A4;A5;A6;A7;A8;A9;A10;A11;A16a;A16b;A17;A6b;A7b;A18;A19];
            vecsense=[repmat('=',Mets,1);repmat('=',NIrr+NRev,1);repmat('<',2*NRev,1);repmat('=',NRev+NIrr_D+2*NRev_D,1);'<';'>';repmat('=',1,1);repmat('>',2*NRev,1);'<';'>'];
            b=[zeros(Mets,1);VpIrr;zeros(NRev,1);zeros(NRev,1)+eps;Vss2max(NIrr+1:Rxns)+eps;VpRev;DIrr;zeros(NRev_D,1);DRev;ZE+deltaZE*ZE;ZE-deltaZE*ZE;ZV;zeros(NRev,1)-eps;Vss2min([Rev_D;Rev_nD])-eps;0;0];
        elseif ~isfield(Options,'Carb2Oxy') || Options.Carb2Oxy=='F',
            Amat=[A1;A4;A5;A6;A7;A8;A9;A10;A11;A16a;A16b;A17;A6b;A7b];
        end
        
        %Construction of Q matrix:

        Q=zeros(size(Amat,2),size(Amat,2));
        cvec=zeros(size(Amat,2),1);cvec((2*Rxns+1):(4*Rxns+NRev))=1;
 
        %Solve MILP     

        m.Q=sparse(Q);
        m.obj=cvec;
        m.A=sparse(Amat);
        m.rhs=b;
        m.sense=vecsense;
        m.modelsense='min';
        m.vtype=[repmat('C',4*Rxns+NRev+2*NIrr_D+4*NRev_D,1);repmat('B',NRev,1)];
        m.lb=lb;
        m.ub=ub;
        params.OutputFlag=OutFlag;
        params.FeasibilityTol=FeasTol;
        params.Presolve=2;
        params.TimeLimit=time_limit;
        params.Threads=8;

        try
            gur=gurobi(m,params);
            X=gur.x; 
            
            VIrrD=X(1:NIrr_D);
            VIrrnD=X((NIrr_D+1):NIrr);
            VForD=X((NIrr+1):(NIrr+NRev_D));
            VFornD=X((NIrr+NRev_D+1):(NIrr+NRev));
            VRevD=X((Rxns+NIrr+1):(Rxns+NIrr+NRev_D));
            VRevnD=X((Rxns+NIrr+NRev_D+1):(Rxns+NIrr+NRev));
            SolV(Irr_D,j)=VIrrD;SolV(Irr_nD,j)=VIrrnD;SolV(Rev_D,j)=VForD-VRevD;SolV(Rev_nD,j)=VFornD-VRevnD;
            SolV(flipped,j)=-SolV(flipped,j);
            Etotal(j)=sum(W'.*abs(X((4*Rxns+NRev+1):(4*Rxns+NRev+2*NIrr_D+4*NRev_D))));
            Vsstotal(j)=sum(abs(SolV(:,j)));
            Vp(flipped)=-Vp(flipped);Vptotal(:,j)=Vp;
            counter=counter+1;
            j=j+1;
        catch
            if OutFlag==1,
               sprintf('Failed attempt to solve QP in n=%d',j);
            end
        end
    end
    
    SolVT=[SolVT,SolV];SolVT(abs(SolVT)<eps)=0;
    for i=1:length(Irr),
        SolVT(Irr(i),SolVT(Irr(i),:)<0)=0;
    end
    EtotalT=[EtotalT;Etotal];
    VsstotalT=[VsstotalT;Vsstotal];
   
   end
   
    close(wbar)
    Val.Vsample=SolVT;
    Val.Etotal=EtotalT;
    Val.Vsstotal=VsstotalT;
    Val.Time=toc;
    
end
