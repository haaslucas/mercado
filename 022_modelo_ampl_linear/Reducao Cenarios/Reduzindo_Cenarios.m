clc
clear all
close all

Reducted_Scenarios = 5;
Original_Scenarios = 100;

NominalSolar = [0 0 0 0 0 0 0 171.1034  416.4702  652.4412  839.1873  957.9717  998.6655  957.9717  839.1873  652.4412  416.4702  171.1034 0 0 0 0 0 0 0];
NominalWind =[0.29902	0.310306667	0.31622	0.316853333	0.315206667	0.3127	0.34776	0.414133333	0.38066	0.2872	0.200526667	0.145333333	0.114113333	0.097993333	0.092666667	0.096326667	0.10932	0.12884	0.175366667	0.217853333	0.245026667	0.2628	0.274993333	0.28666];

for h =1:24
    for k=1:Original_Scenarios
        wind(h,k)=(1+normrnd(0,0.1))*NominalWind(h);
        wind2(h,k)=(1+normrnd(0,0.1))*NominalWind(h);
        load1(h,k) = 1+normrnd(0,0.15);
        load2(h,k) = 1+normrnd(0,0.05);
        solar(h,k) = ((1+normrnd(0,0.1))*NominalSolar(h))/1000;
    end
end
variables={wind, wind2, load1, load2, solar};
T = size(variables{1},1);
N = Reducted_Scenarios;
[S,P,J,L] = scenred(variables, 'cityblock','nodes',round(linspace(1,N,T)));

 clear wind wind2 load1 load2 solar;
 wind=S(:,:,1);
 wind2=S(:,:,2);
 load1=S(:,:,3);
 load2=S(:,:,4);
 solar=S(:,:,5);
 
%  figure()
%  hold on
%  for k=1:5
%      
%      plot(wind(:,k));
%  end


fileMatrixA = fopen('Scenarios.dat','w');
fprintf(fileMatrixA,'set S := ');

formatSpec1 = '%d ';
for k =1:Reducted_Scenarios
    fprintf(fileMatrixA,formatSpec1,k);
end
fprintf(fileMatrixA,';\n\n\n');




formatSpec = '%d %d %6.2f \n';

fprintf(fileMatrixA,'param wind1 := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,wind(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');

fprintf(fileMatrixA,'param wind2 := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,wind2(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');


fprintf(fileMatrixA,'param load1 := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,load1(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');


fprintf(fileMatrixA,'param load2 := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,load2(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');


fprintf(fileMatrixA,'param solar := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,solar(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');



for t=1:T
    qtd = 0;
    pos=[];
    for k =1:Reducted_Scenarios
        if P(t,k)==0
            qtd = qtd + 1;
            pos(qtd)=k;
        end
    end
    
    for k=1:size(pos,2)
        P(t,pos(k)) = max(P(t,:))/(qtd+1);
    end
    [value, position]=max(P(t,:));
    P(t,position)=value/(qtd+1);
end


a=[55*30/195
55*60/195
55*45/195
55*40/195
55*20/195];

formatSpec = '%d %d %6.5f \n';
fprintf(fileMatrixA,'param probability := \n');

for k =1:Reducted_Scenarios
    for t=1:T
        fprintf(fileMatrixA,formatSpec,k,t,P(t,k));
    end
end
fprintf(fileMatrixA,';\n\n\n');
fclose(fileMatrixA);

