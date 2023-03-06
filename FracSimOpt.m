clc
clear

zigma=0.01;

KMATRIX=5;  %% UNIT=========    SQUARE MICRO METER
PORMAT=0.1;
SwcM=0.15;
SorM=0.1;
noM=2;
nwM=2;
Krw0M=0.5;
Kro0M=0.8;


%%             %%%%%%%%%%%%%%%%%%%%%%%  Fracture Characterization
PORFRACT=0.005;
KFRACT=15*KMATRIX;

SwcF=0.05;
SorF=0.05;
noF=1.2;
nwF=1.2;
Krw0F=0.7;
Kro0F=0.9;

%%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%    Simulated Lab Results  or Importing of Original Lab Data

[Time,qo,qw]=QCalc(zigma, KMATRIX, KFRACT, PORMAT, PORFRACT, SwcM, SorM, SwcF, SorF, noM, nwM, noF, nwF, Kro0M, Krw0M, Kro0F, Krw0F);
[BreakThTime,Stablized_Flow]=Characterize(Time, qo, qw);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Initial Guess for Optimization

zigma_FIRST_GUESS=zigma;
%K_FIRST_GUESS= 5*KMATRIX/(((1-(Stablized_Flow/max(qo)))^-1)-1);
K_FIRST_GUESS=15*KMATRIX;
%POR_FIRST_GUESS= max(qo)* BreakThTime/0.006/0.1/0.03;
POR_FIRST_GUESS= PORFRACT;
SwcF_FIRST_GUESS=SwcF;
SorF_FIRST_GUESS=SorF;
noF_FIRST_GUESS=1.2;
nwF_FIRST_GUESS=1.2;
Krw0F_FIRST_GUESS=Krw0F;
Kro0F_FIRST_GUESS=Kro0F;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

UNKNOWN_0(1,1)=zigma_FIRST_GUESS;
UNKNOWN_0(2,1)=K_FIRST_GUESS;
UNKNOWN_0(3,1)=POR_FIRST_GUESS;
UNKNOWN_0(4,1)=SwcF_FIRST_GUESS;
UNKNOWN_0(5,1)=SorF_FIRST_GUESS;
UNKNOWN_0(6,1)=noF_FIRST_GUESS;
UNKNOWN_0(7,1)=nwF_FIRST_GUESS;
UNKNOWN_0(8,1)=Kro0F_FIRST_GUESS;
UNKNOWN_0(9,1)=Krw0F_FIRST_GUESS;


[Time_ES_0,QEstimated_0,qw_ES_0]=QCalc(zigma_FIRST_GUESS, KMATRIX, K_FIRST_GUESS, PORMAT, POR_FIRST_GUESS, SwcM, SorM, SwcF_FIRST_GUESS, SorF_FIRST_GUESS, noM, nwM, noF_FIRST_GUESS, nwF_FIRST_GUESS, Kro0M, Krw0M, Kro0F_FIRST_GUESS, Krw0F_FIRST_GUESS);
[BT_Time_0,Stablized_Flow_0]=Characterize(Time_ES_0, QEstimated_0, qw_ES_0);


WeightNo1=1;  % Weight of Break Through Time Value
WeightNo2=1;   % Weight of Value of Stablized Flow after First Break Through

Error=sum( ( (QEstimated_0-qo)/qo).^2 )+ ...
       WeightNo1*( ((BT_Time_0 - BreakThTime)/BreakThTime)^2  ) + ...
       WeightNo2*( ((Stablized_Flow_0 - Stablized_Flow)/Stablized_Flow) ^2);



%%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%    Optimization


count=0;
LAN=1;
COUNTER=1;
HISTORY(:,COUNTER)=UNKNOWN_0;
g=1;
Error=1;
while Error>10^-5 && abs(max(g))>0.000000001
    for i=1:9
        UNKNOWN_1=UNKNOWN_0;
        UNKNOWN__1=UNKNOWN_0;
        UNKNOWN_1(i)=UNKNOWN_1(i)+.0001*UNKNOWN_1(i);
        UNKNOWN__1(i)=UNKNOWN__1(i)-.0001*UNKNOWN__1(i);
        [Time_ES_1,QEstimated_1,qw_ES_1]=QCalc(UNKNOWN_1(1), KMATRIX, UNKNOWN_1(2), PORMAT, UNKNOWN_1(3), SwcM, SorM, UNKNOWN_1(4), UNKNOWN_1(5), noM, nwM, UNKNOWN_1(6), UNKNOWN_1(7), Kro0M, Krw0M, UNKNOWN_1(8), UNKNOWN_1(9));
        [BT_Time_1,Stablized_Flow_1]=Characterize(Time_ES_1, QEstimated_1, qw_ES_1);

        [Time_ES__1,QEstimated__1,qw_ES__1]=QCalc(UNKNOWN__1(1), KMATRIX, UNKNOWN__1(2), PORMAT, UNKNOWN__1(3), SwcM, SorM, UNKNOWN__1(4), UNKNOWN__1(5), noM, nwM, UNKNOWN__1(6), UNKNOWN__1(7), Kro0M, Krw0M, UNKNOWN__1(8), UNKNOWN__1(9));
        [BT_Time__1,Stablized_Flow__1]=Characterize(Time_ES__1, QEstimated__1, qw_ES__1);        

        for j=1:length(QEstimated_1)
            jac(j,i)=WeightNo1*(QEstimated_1(j)-QEstimated__1(j))/.0002/UNKNOWN_0(i);
        end
        
        for j=length(QEstimated_1)+1:length(QEstimated_1)+1
            jac(j,i)=WeightNo2*(BT_Time_1-BT_Time__1)/.0002/UNKNOWN_0(i);
        end
        
        
        for j=length(QEstimated_1)+2:length(QEstimated_1)+2
            jac(j,i)=(Stablized_Flow_1-Stablized_Flow__1)/.0002/UNKNOWN_0(i);
        end        
        
    end
    
    OUTPUT(1,1:length(QEstimated_1))=(qo-QEstimated_0);
    
    OUTPUT(1,length(QEstimated_1)+1)=WeightNo1*(BreakThTime-BT_Time_0);   
    OUTPUT(1,length(QEstimated_1)+2)=WeightNo2*(Stablized_Flow-Stablized_Flow_0);    

    g=jac'*OUTPUT';
    %=====================
    H=jac'*jac;
    UNKNOWN_2= UNKNOWN_0-(H+LAN* eye(9))\g;
    
    if UNKNOWN_2(4)>noM
        UNKNOWN_2(4)=noM;
    end
    
    if UNKNOWN_2(4)<1
        UNKNOWN_2(4)=1;
    end
    
    if UNKNOWN_2(5)>nwM
        UNKNOWN_2(5)=nwM;
    end   
    if UNKNOWN_2(5)<1
        UNKNOWN_2(5)=1;
    end    
    
    if UNKNOWN_2(3)>1
        UNKNOWN_2(3)=1;
        
    end

    [Time_ES_2,QEstimated_2,qw_ES_2]=QCalc(UNKNOWN_2(1), KMATRIX, UNKNOWN_2(2), PORMAT, UNKNOWN_2(3), SwcM, SorM, UNKNOWN_2(4), UNKNOWN_2(5), noM, nwM, UNKNOWN_2(6), UNKNOWN_2(7), Kro0M, Krw0M, UNKNOWN_2(8), UNKNOWN_2(9));
    [BT_Time_2,Stablized_Flow_2]=Characterize(Time_ES_2, QEstimated_2, qw_ES_2);
    
    Error2=sum( ( (QEstimated_2-qo)/qo).^2 )+ ...
       WeightNo1*( ((BT_Time_2 - BreakThTime)/BreakThTime)^2  ) + ...
       WeightNo2*( ((Stablized_Flow_2 - Stablized_Flow)/Stablized_Flow) ^2);
   
    if Error2<Error
        UNKNOWN_0=UNKNOWN_2;
        QEstimated_0=QEstimated_2;
        Error=Error2;
        LAN=LAN*.1;
        
    else
        LAN=LAN*10;
        
    end
    RMS=sqrt(2*Error/(length(QEstimated_2)+2))
    
    
    COUNTER=COUNTER+1;
    HISTORY(:,COUNTER)=UNKNOWN_2;
end


errordlg('Finished');