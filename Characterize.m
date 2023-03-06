function [BreakThTime,Stablized_Flow]=Characterize(Time, qo, qw)

Total_Itarations=length(qo);

for i=1:Total_Itarations

    if qw(i)>0.0001
        BreakThTime=Time(i);
        BreakThItaration=i;
        break;
    end
end
[dq]=derivative(Time, qo, qw);
[MAX_DQ MAX_DQ_INDICE] = max(dq,[],2); %Indice and value of maximum rate change
Stablized_Flow_Time=1.5*(Time(MAX_DQ_INDICE)-BreakThTime)+Time(MAX_DQ_INDICE);
Stablized_Flow=interp1(Time,qo,Stablized_Flow_Time);
    
end


