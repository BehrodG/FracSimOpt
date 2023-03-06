clc
clear

KMATRIX=5;
PORMAT=0.1;
SwcM=0.15;
SorM=0.1;
noM=2;
nwM=2;
Krw0M=0.5;
Kro0M=0.8;
%%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%       Fracture Characterization
zigma=0.01;
PORFRACT=0.005;
KFRACT=15*KMATRIX;
SwcF=0.1;
SorF=0.1;
noF=2;
nwF=2;
Krw0F=0.7;
Kro0F=0.9;

%%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%      Simulated Lab Results  or Importing of Original Lab Data

[Time,qo,qw]=QCalc(zigma, KMATRIX, KFRACT, PORMAT, PORFRACT, SwcM, SorM, SwcF, SorF, noM, nwM, noF, nwF, Kro0M, Krw0M, Kro0F, Krw0F);
[BreakThTime,Stablized_Flow]=Characterize(Time, qo, qw);

% save FracResult
% 
% load FracResult

% Min Values
%       Sigam   ,  Permeability  ,  Porosity  ,  Connate Water, Residual   OIL ,  no    ,  nw    ,  Kro Max Fracture,   Krw Max Fracture
MinVal=[0.005   ,       50       ,    0.001   ,     0.05      ,     0.02       ,  0.5   ,  0.5   ,          0.1     ,          0.1       ];

% Max Values
%       Sigam   ,  Permeability  ,  Porosity ,  Connate Water, Residual   OIL ,  no    ,  nw    ,  Kro Max Fracture,   Krw Max Fracture
MaxVal=[0.05    ,      150       ,      0.01 ,        0.7    ,      0.5       ,    4   ,    4   ,            1     ,          1       ];

%%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%      Initializing of Optimization

zigma_FIRST_GUESS=0.01;
K_FIRST_GUESS= 5*KMATRIX/(((1-(Stablized_Flow/max(qo)))^-1)-1);
K_FIRST_GUESS=60;
POR_FIRST_GUESS= max(qo)* BreakThTime/0.006/0.1/0.03;
POR_FIRST_GUESS=0.002;
SwcF_FIRST_GUESS=SwcM;
SorF_FIRST_GUESS=SorM;
noF_FIRST_GUESS=1;
nwF_FIRST_GUESS=1;
Kro0F_FIRST_GUESS=Kro0M;
Krw0F_FIRST_GUESS=Krw0M;

% [Time_ES0,qo_ES0,qw_ES0]=QCalc(zigma, KMATRIX, K_FIRST_GUESS, PORMAT, POR_FIRST_GUESS, SwcM, SorM, SwcF, SorF, noM, nwM, noF, nwF, Kro0M, Krw0M, Kro0F, Krw0F);
% 
% plot(Time,qo_ES-qo);
% figure (2);
% plot(Time,qo,Time,qw, Time_ES, qo_ES,Time_ES, qw_ES);





[Time_ES_0,QEstimated_0,qw_ES_0]=QCalc(zigma_FIRST_GUESS, KMATRIX, K_FIRST_GUESS, PORMAT, POR_FIRST_GUESS, SwcM, SorM, SwcF_FIRST_GUESS, SorF_FIRST_GUESS, noM, nwM, noF_FIRST_GUESS, nwF_FIRST_GUESS, Kro0M, Krw0M, Kro0F_FIRST_GUESS, Krw0F_FIRST_GUESS);
[BT_Time_0,Stablized_Flow_0]=Characterize(Time_ES_0, QEstimated_0, qw_ES_0);


WeightNo1=5;  % Weight of Break Through Time Value
WeightNo2=5;   % Weight of Value of Stablized Flow after First Break Through

Error=sum( ( (QEstimated_0-qo)/qo).^2 )+ ...
       WeightNo1*( ((BT_Time_0 - BreakThTime)/BreakThTime)^2  ) + ...
       WeightNo2*( ((Stablized_Flow_0 - Stablized_Flow)/Stablized_Flow) ^2);



UNKNOWN_0(1,1)=zigma_FIRST_GUESS;
UNKNOWN_0(2,1)=K_FIRST_GUESS;
UNKNOWN_0(3,1)=POR_FIRST_GUESS;
UNKNOWN_0(4,1)=SwcF_FIRST_GUESS;
UNKNOWN_0(5,1)=SorF_FIRST_GUESS;
UNKNOWN_0(6,1)=noF_FIRST_GUESS;
UNKNOWN_0(7,1)=nwF_FIRST_GUESS;
UNKNOWN_0(8,1)=Kro0F_FIRST_GUESS;
UNKNOWN_0(9,1)=Krw0F_FIRST_GUESS;

for i=1:9
    UNKNOWN_0(i)=(UNKNOWN_0(i)-MinVal(i))/(MaxVal(i)-MinVal(i));  %Normilize Parameters%
end
 
%%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%     Sensitivity of Gradient  

    for i=1:9
        UNKNOWN_1=UNKNOWN_0;

        UNKNOWN_1(i)=UNKNOWN_1(i)+.00001;
        
        for xxx=1:9
            UNKNOWN_1(xxx)=UNKNOWN_1(xxx)*(MaxVal(xxx)-MinVal(xxx))+MinVal(xxx);
        end
        
        [Time_ES_1,QEstimated_1,qw_ES_1]=QCalc(UNKNOWN_1(1), KMATRIX, UNKNOWN_1(2), PORMAT, UNKNOWN_1(3), SwcM, SorM, UNKNOWN_1(4), UNKNOWN_1(5), noM, nwM, UNKNOWN_1(6), UNKNOWN_1(7), Kro0M, Krw0M, UNKNOWN_1(8), UNKNOWN_1(9));
        [BT_Time_1,Stablized_Flow_1]=Characterize(Time_ES_1, QEstimated_1, qw_ES_1);
        
        for j=1:length(QEstimated_1)
            jac(j,i)=WeightNo1*(QEstimated_0(j)-QEstimated_1(j))/QEstimated_0(j)/.00001;
        end
        
        for j=length(QEstimated_1)+1:length(QEstimated_1)+1
            jac(j,i)=WeightNo2*(BT_Time_0-BT_Time_1)/BT_Time_0/0.00001;
        end
        
        
        for j=length(QEstimated_1)+2:length(QEstimated_1)+2
            jac(j,i)=(Stablized_Flow_0-Stablized_Flow_1)/Stablized_Flow_0/.00001;
        end        
        
    end
    
    OUTPUT(1,1:length(QEstimated_1))=1-QEstimated_0./qo;
    
    OUTPUT(1,length(QEstimated_1)+1)=WeightNo1*(BreakThTime-BT_Time_0)/BreakThTime;   
    OUTPUT(1,length(QEstimated_1)+2)=WeightNo2*(Stablized_Flow-Stablized_Flow_0)/Stablized_Flow;    

    g=jac'*OUTPUT';
    %=====================
    H=jac'*jac;
    
    
%%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%     VISE OPTIMIZATION    

[Gradinet_Sort , Gradient_ParamOrder]=sort(abs(g));  % sort parameters with low to high gradients%
NOAP=1;       %Number Of Active Parameters%

ActParam=zeros(9,1);
for i=9:-1:(9-NOAP+1)
    ActParam(Gradient_ParamOrder(i))=1;
end
   

    
count=0;
LAN=1;
COUNTER=1;
HISTORY(:,COUNTER)=UNKNOWN_0;
g=1;

while Error>10^-3 && max(abs(g))>0.00001
    jac(:,:)=0;
    for i=1:9
        UNKNOWN_1=UNKNOWN_0;
        
        if ActParam(i)==1
            UNKNOWN_1(i)=UNKNOWN_1(i)+.00001;
        end
        
        for xxx=1:9
            UNKNOWN_1(xxx)=UNKNOWN_1(xxx)*(MaxVal(xxx)-MinVal(xxx))+MinVal(xxx);
        end
        
        if ActParam(i)==1
            [Time_ES_1,QEstimated_1,qw_ES_1]=QCalc(UNKNOWN_1(1), KMATRIX, UNKNOWN_1(2), PORMAT, UNKNOWN_1(3), SwcM, SorM, UNKNOWN_1(4), UNKNOWN_1(5), noM, nwM, UNKNOWN_1(6), UNKNOWN_1(7), Kro0M, Krw0M, UNKNOWN_1(8), UNKNOWN_1(9));
            [BT_Time_1,Stablized_Flow_1]=Characterize(Time_ES_1, QEstimated_1, qw_ES_1);
        
            for j=1:length(QEstimated_1)
                jac(j,i)=WeightNo1*(QEstimated_0(j)-QEstimated_1(j))/QEstimated_0(j)/.00001;
            end
        
            for j=length(QEstimated_1)+1:length(QEstimated_1)+1
                jac(j,i)=WeightNo2*(BT_Time_0-BT_Time_1)/BT_Time_0/0.00001;
            end
        
        
            for j=length(QEstimated_1)+2:length(QEstimated_1)+2
                jac(j,i)=(Stablized_Flow_0-Stablized_Flow_1)/Stablized_Flow_0/.00001;
            end
        end
        
    end
    
    OUTPUT(1,1:length(QEstimated_1))=1-QEstimated_0./qo;
    
    OUTPUT(1,length(QEstimated_1)+1)=WeightNo1*(BreakThTime-BT_Time_0)/BreakThTime;   
    OUTPUT(1,length(QEstimated_1)+2)=WeightNo2*(Stablized_Flow-Stablized_Flow_0)/Stablized_Flow;    

    g=jac'*OUTPUT';
    %=====================
    H=jac'*jac;
    
    
    
    H_Active=zeros(NOAP,NOAP);  
    iii=1;
    jjj=1;
    
    for i=1:9
            if ActParam(i)==1
                for j=1:9
                    if ActParam(j)==1
                        H_Active(iii,jjj)=H(i,j);
                        jjj=jjj+1;
                    end
                end
                jjj=1;
                iii=iii+1;
            end
    end
    
    
    
    
    CovMat_Active=inv(H_Active+LAN*eye(NOAP));
    
    CovMat=zeros(9,9);
    
    iii=1;
    jjj=1;
       
    for i=1:9
            if ActParam(i)==1
                for j=1:9
                    if ActParam(j)==1
                        CovMat(i,j)=CovMat_Active(iii,jjj);
                        jjj=jjj+1;
                    end
                end
                jjj=1;
                iii=iii+1;
            end
    end    
    
    Delta=-CovMat*g;
    UNKNOWN_2_NOR= UNKNOWN_0+Delta;  %UNKNOWN_2 Normilized%
    
    for i=1:9
        if UNKNOWN_2_NOR(i)<0
            UNKNOWN_2_NOR(i)=0;
        end
        if UNKNOWN_2_NOR(i)>1
            UNKNOWN_2_NOR(i)=1;
        end
    end
    for xxx=1:9
        UNKNOWN_2(xxx)=UNKNOWN_2_NOR(xxx)*(MaxVal(xxx)-MinVal(xxx))+MinVal(xxx); %DeNormilize%
    end    
        
    
    [Time_ES_2,QEstimated_2,qw_ES_2]=QCalc(UNKNOWN_2(1), KMATRIX, UNKNOWN_2(2), PORMAT, UNKNOWN_2(3), SwcM, SorM, UNKNOWN_2(4), UNKNOWN_2(5), noM, nwM, UNKNOWN_2(6), UNKNOWN_2(7), Kro0M, Krw0M, UNKNOWN_2(8), UNKNOWN_2(9));
    [BT_Time_2,Stablized_Flow_2]=Characterize(Time_ES_2, QEstimated_2, qw_ES_2);
    
    Error2=sum( ( (QEstimated_2-qo)/qo).^2 )+ ...
           WeightNo1*( ((BT_Time_2 - BreakThTime)/BreakThTime)^2  ) + ...
           WeightNo2*( ((Stablized_Flow_2 - Stablized_Flow)/Stablized_Flow) ^2);
   
    if Error2<Error
        UNKNOWN_0=UNKNOWN_2_NOR;
        QEstimated_0=QEstimated_2;
        Error=Error2;
        LAN=LAN*0.25;
        
    else
        LAN=LAN*2;
        
    end
    RMS=sqrt(2*Error/(length(QEstimated_2)+2));
    Error
    Error2
    UNKNOWN_2
    errrror(COUNTER)=Error;
    
    
    COUNTER=COUNTER+1;
    HISTORY(:,COUNTER)=UNKNOWN_2_NOR;
end


errordlg('Finished');