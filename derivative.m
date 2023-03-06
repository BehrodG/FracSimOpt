function [dq]=derivative(Time, qo, qw)

Total_Itarations=length(qo);
dq(1)=0;
for i=2:Total_Itarations

dq(i)=(qw(i)-qw(i-1))/(Time(i)-Time(i-1));
end

end


