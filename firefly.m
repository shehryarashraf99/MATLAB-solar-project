function[tilt,power]=firefly(x0,objective,num_fireflies,alpha,beta0,gamma,theta,Max_Iterations,NVARS,lb,ub)

%% hyperparamters
% num_fireflies=10; %num fireflies
% alpha=1;
% beta0=1;
% gamma=0.01;
% theta=0.97;

% Max_Iterations=5;

fbest=zeros(1,num_fireflies);
xbest=cell(1,num_fireflies);


%% aditional function to ensure boundary conditions are met

function [ns]=findlimits(n,ns,Lb,Ub)
for l=1:n
  nsol_tmp=ns(l,:);
  % Apply the lower bound
  I=nsol_tmp<Lb;  nsol_tmp(I)=Lb(I);
  % Apply the upper bounds
  J=nsol_tmp>Ub;  nsol_tmp(J)=Ub(J);
  % Update this new move
  ns(l,:)=nsol_tmp;
end
end

%% calculate initial function values here

parfor i=1:num_fireflies
    %ns(i,:)=lb+(ub-lb).*rand(1,NVARS); % random initialization within bounds (for now)
    lightning(i)=objective([x0(i,:)]); %function evaluation at initial point
    
end    
%% main function loop here.    

for k=1:Max_Iterations
    alpha=alpha*theta;
    scale=abs(ub-lb);
    for i=1:num_fireflies
        for j=1:num_fireflies
            lightning(i)=objective([x0(i,:)]);
            if lightning(i)<=lightning(j) % we are doing minimization, hence update
                r=sqrt(sum((x0(i,:)-x0(j,:)).^2));
                beta=beta0*exp(-gamma*r.^2);
                steps=alpha.*(rand(1,NVARS)-0.5).*scale;
                x0(i,:)=x0(i,:)+beta*(x0(j,:)-x0(i,:))+steps;
            end
        end
    end
  %% I want to see num_function evaluations
 x=findlimits(num_fireflies,x0,lb,ub);
 [Lightn,Index]=sort(lightning,'descend');
temp=x;
for u=1:num_fireflies
 x(u,:)=temp(Index(u),:);
end
fbest(k)=Lightn(1); 
xbest{k}=x(1,:);
    
    
end
%% sort fbest and xbest here

[best,index]=sort(fbest,'descend');
power=best;
tilt=xbest{index};





    
end       
                


