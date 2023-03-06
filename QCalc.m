function [Time,qo,qw]=QCalc(zigma, KMATRIX, KFRACT, PORMAT, PORFRACT, SwcM, SorM, SwcF, SorF, noM, nwM, noF, nwF, Kro0M, Krw0M, Kro0F, Krw0F)


% Injection term is included in governing equation of fracture blocks
% adjacent to injection area


%  No well in end blocks. One hypothetical block considered after final
%  blocks in each raw and pressure of it is assumed equal to atmospheric
%  pressure. Moreover, there is flow from end fracture blocks to
%  atmospheric pressure

Tyfo=zeros(18,1); %Transmissibility in Y, Forward direction for Oil
Tyfw=zeros(18,1); %Transmissibility in Y, Forward direction for Water
Txfo=zeros(18,1); %Transmissibility in X, Forward direction for Oil
Txfw=zeros(18,1); %Transmissibility in X, Forward direction for Water
Tybo=zeros(18,1); %Transmissibility in Y, Backward direction for Oil
Tybw=zeros(18,1); %Transmissibility in Y, Backward direction for Water
Txbo=zeros(18,1); %Transmissibility in X, Backward direction for Oil
Txbw=zeros(18,1); %Transmissibility in X, Backward direction for Water




cpoo=zeros(18,1);
cswo=zeros(18,1);
cpow=zeros(18,1);
csww=zeros(18,1);

Porosity=zeros(18,1);
KD=zeros(18,1);
LandaWater=zeros(18,1);
LandaOil=zeros(18,1);

Bw=zeros(18,1);
DiffBw=zeros(18,1);

Bo=zeros(18,1);
DiffBo=zeros(18,1);

Visw=zeros(18,1);
Viso=zeros(18,1);

krw=zeros(18,1);
kro=zeros(18,1);

DiffPc=zeros(18,1);
Pcow=zeros(18,1);

p=zeros(18,1);
so=zeros(18,1);

Alpha=zeros(18,1);

R=zeros(10,1);
RForward=zeros(18,1);
RBackward=zeros(18,1);
RCentral=zeros(18,1);

BlockGIP=zeros(18,1);
BlockVol=zeros(18,1);

ZZ=zeros(20,1);
ZForward=zeros(18,1);
ZBackward=zeros(18,1);
ZCentral=zeros(18,1);


RAxis=zeros(20,10);
ZAxis=zeros(20,10);
S=zeros(6,3);

sw=zeros(18,1);
swnew=zeros(18,1);

ACTNUM=zeros(18,1);

BlockOIP=zeros(18,1);
BlockWIP=zeros(18,1);


WC=zeros(10000,1);

%% Reservoir Properties


dX(1:18)=0.1/3;     % Unit =======      M E T E R



dY(1:18)=0.01;
dZ=0.03; 

ACTNUM(1:9)=1;


% ACTNUM(101:200)=[1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1 ...
%                  1  1   1   1   1   1   1   1   1   1  ];
             
% %              
ACTNUM(10:18)=  [0  0   0 ...
                 1  1   1 ...
                0   0   0  ];

% % 
% ACTNUM(26:50)=  [0  0   0   0   0 ...
%                  0  0   0   0   0 ...
%                  0  0   0   0   0 ...
%                  0  0   0   0   0 ...
%                  0  0   0   0   0  ];

TOTACTFRCELL=sum(ACTNUM(10:18)); 
ACTFRCELLNO=1;


for i=10:18
    if ACTNUM(i)==1 
        ACTFRCELL(ACTFRCELLNO)=i;
        ACTFRCELLNO=ACTFRCELLNO+1;
    end
end








%Determination of Boundary of each block
% 


porosity(1:9,1)=PORMAT; %               UNIT==========   UNITLESS
porosity(10:18,1)=PORFRACT;

KD(1:9,1)=KMATRIX; %% Permeability Assignment         UNIT=========    SQUARE MICRO METER
KD(10:18,1)=KFRACT;

kX=KD;  % PERMEABILITY IN X DIRECTION
kY=kX;  % PERMEABILITY IN Y DIRECTION


Cr=3*0.1450377 * 1e-6; %rock compressibility            ( UNIT =====   KPA^-1 )     
% =   3E-6 psi^-1



%% PVT Table

	
	PVT=[101.325	1        0.001;...
        6894.757	1.03     0.0013;...
        20684.271   1.1      0.0014;];

% Fluids Properties

Reg3= Polyfit(PVT(:,1),PVT(:,2),2); % Bo

Reg5= Polyfit(PVT(:,1),PVT(:,3),2); % Viso


VisoFunc=@(x) Reg5(1)*x.^2 + Reg5(2)*x + Reg5(3);


% BwFunc=@(x) Reg2(1)*x.^2 + Reg2(2)*x + Reg2(3);

BoFunc=@(x) Reg3(1)*x.^2 + Reg3(2)*x + Reg3(3);
DiffBoFunc=@(x) Reg3(1)*x*2 + Reg3(2);


%% Rock Properties

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rock Properties
SwM=[SwcM,1-SorM ];
PcM=[0.9383 0];
p3=polyfit(SwM,PcM,1);
%Determining of Capillary Pressure & Diff of Pc
Pc=@(x) p3(1)*x + p3(2);
diffPc=p3(1);


%Determining of Relative Permabilities
KrwFuncM=@(x) Krw0M*(x^nwM);
KroFuncM=@(x) Kro0M*((1-x)^noM);

KrwFuncF=@(x) Krw0F*(x^nwF);
KroFuncF=@(x) Kro0F*((1-x)^noF);

%% INITIAL CONDITION 

for i=1:9
    p(i)=101.325; %psi
    sw(i)=SwcM;
end

for i=10:18
    p(i)=101.325; %psi
    sw(i)=SwcF;
end
%%   Well Properties

Pbh=101.325; %KPa

qi=0.1;      %UNIT ============    CUBIC-MTER/DAY for each well

KHINJECTION=0;
for i=1:18
    if mod(i,3)==1 && ACTNUM(i)==1
        KHINJECTION=KHINJECTION+KD(i);
    end
end

for i=1:18
    if mod(i,3)==1 && ACTNUM(i)==1
        qii(i)=qi*KD(i)/KHINJECTION;
    end
end
%% DATE & Time STEPS

dt=zeros(100000,1);
dt(:,1)=0.00000001;
dt(1:500,1)=0.00000001;
dt(10000:1000000)=0.00000001;

%% ****************************** BLOCK VOLUMES & INITIAL Oil IN PLACE ********************************
for i=1:18
BlockVol(i)=dZ*dX(i)*dY(i);
end

BlockVol=BlockVol; %%%CUBIC M E T E R

        

for i=1:18         
  Bo(i)=BoFunc(p(i)); 
  Bw(i)=1;
end



%% ****************************** MAIN LOOP ********************************

time=0;
counter=1;
pold=p;
swold=sw;
WC(1)=0;

while (counter<=1000 ) && (max(WC) <= 0.95)

      
    time=time+dt(counter);
    Time(counter)=time;
    error=100;
    %IMPES
    while error>=5
        %***************** relative permeabiliTz data *********************
        for i=1:9
        NormilizedSw(i)=((swold(i)-SwcM)/(1-SorM-SwcM));
        end        
        for i=10:18
        NormilizedSw(i)=((swold(i)-SwcF)/(1-SorF-SwcF));
        end 
        for i=1:18
            if NormilizedSw(i)<0
                NormilizedSw(i)=0;
            end
        end
        
        for i=1:9
            kro(i)=KroFuncM(NormilizedSw(i));
            krw(i)=KrwFuncM(NormilizedSw(i));
        end
        
        for i=10:18
            kro(i)=KroFuncF(NormilizedSw(i));
            krw(i)=KrwFuncF(NormilizedSw(i));
        end   
        
         for i=1:9
            Pcow(i)=Pc(swold(i));
            DiffPc(i)=diffPc;
        end                   
        
        for i=1:18
            if kro(i)>1
                kro(i)=1;
            end
            
            if kro(i)<0
                kro(i)=0;
            end
            
            if krw(i)>1
                krw(i)=1;
            end
            
            if krw(i)<0
                krw(i)=0;
            end

            
        end
   
        
         for i=1:18
             
%           Bw(i)=BwFunc(pold(i));
            Bw(i)=1;      
%           DiffBw(i)=((1/BwFunc(1.00001*pold(i)))-(1/BwFunc(0.99999*pold(i))))/...
%           (0.00002*pold(i));
            DiffBw(i)=0.00004;

            Bo(i)=BoFunc(pold(i));
            DiffBo(i)=DiffBoFunc(pold(i));
%           Visw(i)=ViswFunc(pold(i));
            Visw(i)=0.001;
            Viso(i)=VisoFunc(pold(i));
            
         end
        
        for i=1:18

           LandaWater(i)=krw(i)/(Visw(i)*Bw(i));   %%%%%%%%%%%%  ????? Bo and Bw
           LandaOil(i)=kro(i)/(Viso(i)*Bo(i)); 
        
        end
        


        %**************** Calculation of Trb, Trf, Tzb, Tzf ***********************
        for i=1:18
            if mod(i,3)==1
                Gxb(i)=0;
            else
                Gxb(i)=ACTNUM(i)*ACTNUM(i-1)*( (dZ*dY(i) )/ ( ( (dX(i)/2)/kX(i))  + ( (dX(i-1)/2)/kX(i-1))  ) );
            end
  
            if mod(i,3)==0
                Gxf(i)=ACTNUM(i)*ACTNUM(i)*( (dZ*dY(i) )/ ( ( (dX(i)/2)/kX(i)) + ( (dX(i)/2)/kX(i))  ) );
            else
                Gxf(i)=ACTNUM(i)*ACTNUM(i+1)*( (dZ*dY(i) )/ ( ( (dX(i)/2)/kX(i))  + ( (dX(i+1)/2)/kX(i+1))  ) );
            end

            if (1<=i && i<=3) || (10<=i && i<=12) 
                Gyb(i)=0;
            else
                Gyb(i)=ACTNUM(i)*ACTNUM(i-3)*( (dZ*dX(i) )/ ( ( (dY(i)/2)/kY(i))  + ( (dY(i-3)/2)/kY(i-3))  ) );
            end

            if (7<=i && i<=9) || (16<=i && i<=18)
                Gyf(i)=0;
            else
                Gyf(i)=ACTNUM(i)*ACTNUM(i+3)*( (dZ*dX(i) )/ ( ( (dY(i)/2)/kY(i))  + ( (dY(i+3)/2)/kY(i+3))  ) );
            end
        end



        
        
          for i=1:18
            if mod(i,3)==1
                Landa_Water_Back_R(i)=0;
                Landa_Oil_Back_R(i)=0;
            else
                if p(i) >= p(i-1)
                    Landa_Water_Back_R(i)=LandaWater(i);
                    Landa_Oil_Back_R(i)=LandaOil(i);
                else
                    Landa_Water_Back_R(i)=LandaWater(i-1);
                    Landa_Oil_Back_R(i)=LandaOil(i-1);
                end
            end
            

            if mod(i,3)==0
                if p(i) >= Pbh
                    Landa_Water_Forward_R(i)=LandaWater(i);
                    Landa_Oil_Forward_R(i)=LandaOil(i);
                else
                    Landa_Water_Forward_R(i)=LandaWater(i);
                    Landa_Oil_Forward_R(i)=LandaOil(i);
                end
            else
                if p(i) >= p(i+1)
                    Landa_Water_Forward_R(i)=LandaWater(i);
                    Landa_Oil_Forward_R(i)=LandaOil(i);
                else
                    Landa_Water_Forward_R(i)=LandaWater(i+1);
                    Landa_Oil_Forward_R(i)=LandaOil(i+1);
                end
            end

            if (1<=i && i<=3) || (10<=i && i<=12) 
                Landa_Water_Back_Z(i)=0;
                Landa_Oil_Back_Z(i)=0;
            else
                if p(i) >= p(i-3)
                Landa_Water_Back_Z(i)=LandaWater(i);
                Landa_Oil_Back_Z(i)=LandaOil(i);
                else
                Landa_Water_Back_Z(i)=LandaWater(i-3);
                Landa_Oil_Back_Z(i)=LandaOil(i-3);
                end
           end

            if (7<=i && i<=9) || (16<=i && i<=18 ) 
                Landa_Water_Forward_Z(i)=0;
                Landa_Oil_Forward_Z(i)=0;
            else
                if p(i) >= p(i+3)
                Landa_Water_Forward_Z(i)=LandaWater(i);
                Landa_Oil_Forward_Z(i)=LandaOil(i);
                else
                Landa_Water_Forward_Z(i)=LandaWater(i+3);
                Landa_Oil_Forward_Z(i)=LandaOil(i+3);
                end
            end
            
          end

       for i=1:18
            
                Txbw(i)=86.4E-6* Landa_Water_Back_R(i)*Gxb(i);
                Txbo(i)=86.4E-6* Landa_Oil_Back_R(i)*Gxb(i);
            
                Txfw(i)=86.4E-6* Landa_Water_Forward_R(i)*Gxf(i);
                Txfo(i)=86.4E-6* Landa_Oil_Forward_R(i)*Gxf(i);

                Tybw(i)=86.4E-6* Landa_Water_Back_Z(i)*Gyb(i);
                Tybo(i)=86.4E-6* Landa_Oil_Back_Z(i)*Gyb(i);

                Tyfw(i)=86.4E-6* Landa_Water_Forward_Z(i)*Gyf(i);
                Tyfo(i)=86.4E-6* Landa_Oil_Forward_Z(i)*Gyf(i);
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrix - Fracture Transmissilibilty  %%%
 
        for i=1:9
                if p(i) >= p(i+9)
                    Landa_Water_FM(i)=LandaWater(i);
                    Landa_Oil_FM(i)=LandaOil(i);
                else
                    Landa_Water_FM(i)=LandaWater(i+9);
                    Landa_Oil_FM(i)=LandaOil(i+9);
                end
            
                TFMw(i)=ACTNUM(i)*ACTNUM(i+9)*86.4E-6* Landa_Water_FM(i);
                TFMo(i)=ACTNUM(i)*ACTNUM(i+9)*86.4E-6* Landa_Oil_FM(i);
            
        end  
        
        for i=10:18
            
                if p(i) > p(i-9)
                    Landa_Water_FM(i)=LandaWater(i);
                    Landa_Oil_FM(i)=LandaOil(i);
                else
                    Landa_Water_FM(i)=LandaWater(i-9);
                    Landa_Oil_FM(i)=LandaOil(i-9);
                end
            
                TFMw(i)=ACTNUM(i)*ACTNUM(i-9)*86.4E-6* Landa_Water_FM(i);
                TFMo(i)=ACTNUM(i)*ACTNUM(i-9)*86.4E-6* Landa_Oil_FM(i);
            
        end     
        
        
        Sigma=zigma;
        
        TFMw=TFMw*Sigma;
        TFMo=TFMo*Sigma;
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
        
 

        for i=1:18
          cpow(i)=(BlockVol(i)*porosity(i)*(swold(i))/dt(counter))*(Cr/Bw(i) + DiffBw(i));

          cswo(i)=-(BlockVol(i)*porosity(i))/(Bo(i)*dt(counter)); 
            
          cpoo(i)=(BlockVol(i)*porosity(i)/dt(counter))*  (1-swold(i))*((Cr/Bo(i) )+ DiffBo(i) );
            
          csww(i)=(BlockVol(i)*porosity(i)/(Bw(i)*dt(counter)))- DiffPc(i)*cpow(i);
        end


        %*************************** A*X=C => A & C ?? ****************************
        for i=1:18

            
                Alpha(i)= -cswo(i)/csww(i); 
                

            if i==1
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                

                A(i,i+1) = Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)= Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))...
                        - (Tyfo(i)+ Alpha(i)* (Tyfw(i)))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i));  


            elseif i==10
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                

                A(i,i+1) = Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)= Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))...
                        - (Tyfo(i)+ Alpha(i)* (Tyfw(i)))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION;  
                  

            elseif i==3  
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i) ); 
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i)  );
                A(i,i)=-(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...   %because of hypothetical grid
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i)); 


            elseif i==12      
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i) ); 
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i)  );
                A(i,i)=-(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...   %because of hypothetical grid
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i); 


            elseif i==7
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );               
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i) );
                A(i,i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))...
                        - (Tybo(i)+ Alpha(i)* (Tybw(i)))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i));  

           
            elseif i==16
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );               
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i) );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i) ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i) ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION;
              
            elseif i==9
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i)  );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );  
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...   %because of hypothetical grid
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i));


                 
            elseif i==18
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i)  );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );  
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...   %because of hypothetical grid
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i));
                   
                c(i)=  - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i); 


                  
            elseif (1<i && i<3)
                
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-1)= Txbo(i)+ Alpha(i)* Txbw(i);                
                A(i,i+1)= Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)=Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)=-(Txbo(i)+ Alpha(i)*  Txbw(i))...
                       -(Txfo(i)+ Alpha(i)*  Txfw(i))...
                       -(Tyfo(i)+ Alpha(i)*  Tyfw(i))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i));
                
                
                
            elseif (10<i && i<12)
  
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-1)= Txbo(i)+ Alpha(i)* Txbw(i);                
                A(i,i+1)= Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)=Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)=-(Txbo(i)+ Alpha(i)*  Txbw(i))...
                       -(Txfo(i)+ Alpha(i)*  Txfw(i))...
                       -(Tyfo(i)+ Alpha(i)*  Tyfw(i))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i);
                
                
                
                
            elseif (7<i && i<9)
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i)  );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );                
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i)  );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i));
            
                
            elseif (16<i && i<18)
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i)  );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );                
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i)  );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i);                
                
                  
            elseif mod(i,3)==1 && i<9
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                

                A(i,i-3)= Tybo(i)+ Alpha(i)* Tybw(i);
                A(i,i+1) = Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)= Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)= - (Tybo(i)+ Alpha(i)* (Tybw(i)))...
                        - (Txfo(i)+ Alpha(i)* (Txfw(i)))...
                        - (Tyfo(i)+ Alpha(i)* (Tyfw(i)))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i));
                   
             elseif mod(i,3)==1 && (i>9)
                 
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);
                 

                A(i,i-3)= Tybo(i)+ Alpha(i)* Tybw(i);
                A(i,i+1) = Txfo(i)+ Alpha(i)* Txfw(i);
                A(i,i+3)= Tyfo(i)+ Alpha(i)* Tyfw(i);
                A(i,i)= - (Tybo(i)+ Alpha(i)* (Tybw(i)))...
                        - (Txfo(i)+ Alpha(i)* (Txfw(i)))...
                        - (Tyfo(i)+ Alpha(i)* (Tyfw(i)))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                        - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       - Alpha(i)*qi*ACTNUM(i)*KD(i)/KHINJECTION;                   
                  
    
            elseif mod(i,3)==0  && i<9
                
                
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );                
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i)  );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...   %because of hypothetical grid
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   

                c(i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i)); 

                  
                  
            elseif mod(i,3)==0  && i>9
                
                
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i)  );                
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i)  );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i)  ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i)  ))...   %because of hypothetical grid
                       -(Txfo(i)+ Alpha(i)* (Txfw(i)  ))...
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i)  ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
                   

                c(i)= - (Txfo(i)+ Alpha(i)* (Txfw(i)))* Pbh ...
                      - (cpoo(i)+ Alpha(i)*cpow(i))*p(i);

                
            elseif  i<9
               
                A(i,i+9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i));                
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i));
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i) );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i) ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i) ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i) ))...
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i) ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i)...
                       + Alpha(i)* Txfw(i)*(Pcow(i+1)-Pcow(i))...
                       + Alpha(i)* Txbw(i)*(Pcow(i-1)-Pcow(i))...
                       + Alpha(i)* Tyfw(i)*(Pcow(i+3)-Pcow(i))...
                       + Alpha(i)* Tybw(i)*(Pcow(i-3)-Pcow(i)); 
                
                
            elseif  i>9
          
                A(i,i-9)= TFMo(i)+ Alpha(i)* TFMw(i);                
                
                A(i,i-3)=Tybo(i)+ Alpha(i)* (Tybw(i) );
                A(i,i-1)= Txbo(i)+ Alpha(i)* (Txbw(i));                
                A(i,i+1)= Txfo(i)+ Alpha(i)* (Txfw(i));
                A(i,i+3)=Tyfo(i)+ Alpha(i)* (Tyfw(i) );
                A(i,i)=-(Tybo(i)+ Alpha(i)* (Tybw(i) ))...
                       -(Txbo(i)+ Alpha(i)* (Txbw(i) ))...
                       -(Txfo(i)+ Alpha(i)* (Txfw(i) ))...
                       -(Tyfo(i)+ Alpha(i)* (Tyfw(i) ))...
                        - (TFMo(i)+ Alpha(i)* TFMw(i))...
                       - (cpoo(i)+ Alpha(i)*cpow(i)) ;
   
                c(i)= - (cpoo(i)+ Alpha(i)*cpow(i))*p(i);                 
                
                
                  
            end

        end
        
   for i=1:9      
     CoefMat2(i,:)=A(i,:);
     ConstMat(i,:)=c(:,i);
   end
       
   
   cou=1;
   for i=10:18
       if ACTNUM(i)==1
           CoefMat2(9+cou,:)=A(i,:);
           ConstMat(9+cou,:)=c(:,i);
           cou=cou+1;
       end
   end
   
   
   for i=1:9      
     CoefMat(:,i)=CoefMat2(:,i);
   end
       
   cou=1;
   for i=10:18
       if ACTNUM(i)==1
           CoefMat(:,9+cou)=CoefMat2(:,i);
           cou=cou+1;
       end
   end
           
pnewActive=inv(CoefMat)*ConstMat;
  
%     for i=1:9
%      Xi=(floor((i-1)/3)+1);
%      Yi=(mod(i-1,3)+1);
%      Valuei=pnewActive(i);
%      S(Xi,Yi)=Valuei;
%     end
%     
    
% error;     
% imagesc(R,ZZ,S);
% contour(R,ZZ,S)
% colorbar;
% pause(0.0001);


   for i=1:9      
     pnew(i)=pnewActive(i);
   end

   cou=1;
   for i=10:18
       if ACTNUM(i)==1
           pnew(i)=pnewActive(9+cou);
           cou=cou+1;
       else
           
           pnew(i)=p(i);
           
       end
   end 
   
    for i=1:18
     Xi=(floor((i-1)/3)+1);
     Yi=(mod(i-1,3)+1);
     Valuei=pnew(i);
     S(Xi,Yi)=Valuei;
    end 
        
        for i=1:18
             
            if i==1
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          + Tyfw(i)*(pnew(i+3)-pnew(i))...
                          + TFMw(i)*(pnew(i+9)-pnew(i))...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION...
                          - cpow(i)*(pnew(i)-p(i)) ...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i)));
                      
          
            elseif i==3
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh - pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          + TFMw(i)*(pnew(i+9)-pnew(i))...   
                          -cpow(i)*(pnew(i)-p(i)) ...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i)));
                      
                      
            elseif i==7
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION...
                          -cpow(i)*(pnew(i)-p(i)) ...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i)) );
                      
                      
            elseif i==9
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh -pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) ...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i)) );                      
        
            elseif i==10
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) ...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION);                       

            elseif i==12
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh - pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );   
                      
                    
            elseif i==16
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) ...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION);
                      
                     
                     
            elseif i==18
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh - pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );
                     
                      
                      
                      
            elseif (1<i && i<3) 
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) )...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i)); 
                      
                     
                     
            elseif  (7<i && i<9)
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) )...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i)); 
                      
                      
            elseif (10<i && i<12)
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );  
                      
            elseif (16<i && i<18)
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );                      
                     

            elseif mod(i,3)==1 && (i<9)

                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION...
                          -cpow(i)*(pnew(i)-p(i)) )...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i));
                   
             elseif mod(i,3)==1 && (i>9)

                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...                       
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                           -cpow(i)*(pnew(i)-p(i)) ...
                          + qi*ACTNUM(i)*KD(i)/KHINJECTION);
                      
            elseif mod(i,3)==0 && (i<9)
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh -pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) )...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i));
                      
                      
            elseif mod(i,3)==0 && (i>9)
                 swnew(i)=sw(i)+(1/csww(i))*(Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Txfw(i)*(Pbh - pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );
                      
            elseif i<9
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i+9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) )...
                          + Txfw(i)*(-Pcow(i+1)+Pcow(i))...
                          + Txbw(i)*(-Pcow(i-1)+Pcow(i))...
                          + Tyfw(i)*(-Pcow(i+3)+Pcow(i))...
                          + Tybw(i)*(-Pcow(i-3)+Pcow(i));                     

            else
                 swnew(i)=sw(i)+(1/csww(i))*(Txfw(i)*(pnew(i+1)-pnew(i))...
                          +Txbw(i)*(pnew(i-1)-pnew(i))...
                          +Tyfw(i)*(pnew(i+3)-pnew(i))...
                          +Tybw(i)*(pnew(i-3)-pnew(i))...
                          +TFMw(i)*(pnew(i-9)-pnew(i))...
                          -cpow(i)*(pnew(i)-p(i)) );
            end
        end
        


    for i=1:18
     Xi=(floor((i-1)/3)+1);
     Yi=(mod(i-1,3)+1);
     Valuei=swnew(i);
     S1(Xi,Yi)=Valuei;
    end
S;

       
        error=0;
        for i=1:18
            error1=abs(pnew(i)-pold(i));
            if error1>error
                error=error1;
            end
        end
        %                 error
        pold=pnew;
        swold=swnew;
    end
    
    
    p=pold;
    sw=swnew;
    

qw1=Txfw(3)*(p(3)-Pbh);
qw2=Txfw(6)*(p(6)-Pbh);
qw3=Txfw(9)*(p(9)-Pbh);


qw8=Txfw(15)*(p(15)-Pbh);


qo1=Txfo(3)*(p(3)-Pbh);
qo2=Txfo(6)*(p(6)-Pbh);
qo3=Txfo(9)*(p(9)-Pbh);

qo8=Txfo(15)*(p(15)-Pbh);


qw(counter)=qw1+qw2+qw3+qw8;
qo(counter)=qo1+qo2+qo3+qo8;

WC(counter)=qw(counter)/(qo(counter)+qw(counter));
% WGR1(counter)=LandaWater(1)/LandaOil(1);
    

% if ( mod(counter,100)==0)  
% subplot(2,3,1)    
% 
% pcolor(RAxis(1,:), ZAxis(:,1),S1);figure(gcf)
% colorbar;
% 
% subplot(2,3,2)    
% grifs=[p(1),p(200)];
% pcolor(RAxis, ZAxis, S); figure(gcf)
% colorbar;
% 
% subplot(2,3,3)    
% plot(Time,qg);
% 
% subplot(2,3,4)    
% plot(Time,qw);
% 
% 
% subplot(2,3,5)    
% plot(Time,GIP);
% 
% 
% subplot(2,3,6)    
% plot(Time,WGR1);
% axis([0,max(Time),0,0.0002])
% 
% 
% pause(0.01);
% end

% if ( mod(counter,50)==1)  
% subplot(2,2,1)    
% 
% imagesc(S);
% colorbar;
% 
% subplot(2,2,2)    
% 
% imagesc(S1);
% 
% colorbar;
% 
% subplot(2,2,3)
% plot(Time,qo);
% 
% subplot(2,2,4)
% plot(Time,qw);
% 
% 
% % subplot(2,3,3)    
% % plot(Time,qg);
% % 
% % subplot(2,3,4)    
% % plot(Time,qw);
% % 
% % 
% % subplot(2,3,5)    
% % plot(Time,GIP);
% % 
% % 
% % subplot(2,3,6)    
% % plot(Time,WGR1);
% % axis([0,max(Time),0,0.0002])
% 
% 
% pause(0.001);
% end
% pause(0.000000001);

counter=counter+1;


end

% 
% 
% subplot(2,2,1)    
% 
% imagesc(S);
% colorbar;
% 
% subplot(2,2,2)    
% 
% imagesc(S1);
% 
% colorbar;
% 
% subplot(2,2,3)
% plot(Time,qo);
% 
% subplot(2,2,4)
% plot(Time,qw);
% 
% % time    ,    WGR
%  

end