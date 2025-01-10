
function [tilt,power] = surrogate_opt(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt)


  NVARS=num_trackers;


%% hyper parameters
 modifier=3.9;
a=0;

runs=1; %default is 50

num_particles=50;
% std=0.5/3;
% mean=0;



%%

surrogateoptions = optimoptions('surrogateopt', ...
                      'PlotFcn','surrogateoptplot',... 
                       'UseParallel',true,'MaxFunctionEvaluations',3000);
                    
intcon=[];



                                     
  

fminconoptions=optimoptions('fmincon', ...
                               'FunctionTolerance',10^-24,...
                               'MaxFunctionEvaluations',200000,...
                               'StepTolerance',0.0005,...
                               'OptimalityTolerance',0,...
                               'ConstraintTolerance',0,...
                               'Display','iter',...
                                'UseParallel',true,...
                               'PlotFcn','optimplotfval');
  
                                      
                                                     
 %% initializing tracker tilt angles at the start of optimization. generally best if initialized to 0                                      
                                    
if rem(num_trackers,2)==0
    x0=[-a,a];
    x0=repmat(x0,1,num_trackers/2);
else
    x0=[a,-a];
    x0=repmat(x0,1,floor(num_trackers/2));
    x0=[x0,a];
end

%% optimizating during the morning

 for k=1:length(morning_zenith)
     
     if morning_zenith(k)>50
         fminconoptions.MaxIterations=10;
     else
         fminconoptions.MaxIterations=1;
     end
     
     
   [~,shading_profile]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x0,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);          
   shading_profile=transpose(shading_profile)+shading_profile;  
   temp=max(shading_profile);
 
counter=1;
is_shaded=zeros(1,length(temp)/2);
for i=1:2:length(temp)
    if temp(i)==1 || temp(i+1)==1
        is_shaded(counter)=1;
    else
        is_shaded(counter)=0;
    end
    counter=counter+1;
end
is_shaded=is_shaded*modifier;

   % is_shaded returns an array that tells us if tracker number (i) is either
   % shading or is being shaded by any other tracker. if so, then that
   % tracker will not be allowed to move up than (swing -modifier) degrees. in our
   % case the value is 0.1 degrees. This extra constraint helps us in
   % improving overall energy output, foregoing power output in the short
   % term to get more power as the day progresses.  
     
     
     
   
lb=x0- (swing-is_shaded); % lower bounds of shaded/shading trackers now modified 
parfor i=1:length(lb)
    if lb(i)<min_tilt(i)
        lb(i)=min_tilt(i);
    end
end

ub=x0+ (swing);
parfor i=1:length(ub)
    if ub(i)>max_tilt(i)
        ub(i)=max_tilt(i);
    end
end 

parfor i=1:length(x0)
    if x0(i)>4
        ub(i)=x0(i);
    end
end

%% calculating initial points

A=-4:0.1:4;
[m,n]=ndgrid(A,A);
Z=[m(:),n(:)];
init=cell(1,length(Z));
for i=1:length(Z)
    init{i}=[Z(i,1),Z(i,2)];
 if rem(num_trackers,2)==0
    init{i}=repmat(init{i},1,num_trackers/2);
 else
    init{i}=repmat(init{i},1,floor(num_trackers/2));
    init{i}=[init{i},Z(i,1)];
end
init{i}=x0+init{i};

for j=1:num_trackers
    if init{i}(j)>ub(j)
        init{i}(j)=ub(j);
    elseif init{i}(j)<lb(j)
        init{i}(j)=lb(j);
    end
end

init{i}=repmat(init{i},num_particles,1);

end

% 2461,4101,1,6561
        
     
        
          

First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized.                                     


%is_not_smooth tells us if we can expect any shading when the trackers have
%been tilted up by value swing. this tells us if we can expect any non_smoothness (discountinuities) in the function,
% if is_not_smooth==0 we can speed up the optimization by reducing the
% number ofdatpoints we are sampling from 6561 to 4. if we are using a non
% smooth optimizer, we no longer need to use it and can use a gradient
% based optimizer instead to speed up the process without a loss in power
% output.

[~,is_not_smooth]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,lb,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
 
is_not_smooth=max(is_not_smooth(:));
    
  if is_not_smooth==0
        init={init{2461},init{4101},init{1},init{6561},init{4019},init{3937},init{3855},init{3773}};
  
    



parfor i=1:length(init)
  %options.InitialSwarmMatrix=init{i};
  [fval(i),~]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,init{i}(1,:),NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
%[Z,z]=particleswarm(objective,NVARS,lb,ub,options);
%X{i}=Z;
%fval(i)=z;
 end
 [~,y]=sort(fval,'descend');
 index_list=ones(1,4);
 for i=1:4
     index_list(i)=y(i);
 end
 if is_not_smooth==1
      index_list=index_list(1);
     index_list=repmat(index_list,1,runs);
      else
     index_list=index_list(1);
 end
 

 
 for n=1:length(index_list)
 Y{n}=init{index_list(n)};
 end
 end
 
 if is_not_smooth==1 
index_list=1;
[potential_tilt_angle{1},~]=surrogateopt(objective,lb,ub,intcon,surrogateoptions);

 
else 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=fmincon(objective,Y{m}(1,:),[],[],[],[],lb-0.04,ub+0.04,[],fminconoptions);
end
end


 


%% rounding to the nearest 10th

% many algorithms such as particlewarm, patternsearch, fmincon do not have
% integer contraints. as the accuracy of tracker tilt angles is 0.1 degrees
% we round the answer to the nearest 10th of a degree. also to avoid
% shading we round down. 3.99 degrees is rounded to 3.9 degrees (not 4)
% because of the possibility that shading occurs at 4.0 but if there is no
% shading at 3.99 degrees then there is no shading at 3.9 degrees.

 for i=1:length(index_list)
     rounded_potential_tilt_angle{i}=round(potential_tilt_angle{i},1);
     for j=1:length(ub)
    if  rounded_potential_tilt_angle{i}(j)<potential_tilt_angle{i}(j)
        rounded_potential_tilt_angle{i}(j)= rounded_potential_tilt_angle{i}(j)+0.1;
    end
     end
 end
 


%% rounding to the nearest degree
% for i=1:length(index_list)
% for j=1:length(ub)
% rounded_potential_tilt_angle{i}(j)=ceil(potential_tilt_angle{i}(j));
% end
% end
    
 %% no rounding
% rounded_potential_tilt_angle=potential_tilt_angle;


for i=1:length(index_list)
 [rounded_potential_power(i),~]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,rounded_potential_tilt_angle{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
rounded_potential_power(i)=rounded_potential_power(i)/1000;
end

 [~,v2]=max(rounded_potential_power);
 tilt{k}=rounded_potential_tilt_angle{v2};
 power(k)=rounded_potential_power(v2);
 x0=tilt{k};
 
%  if runs>10
%     runs=runs-2;
%  end
 
 tilt{k}
  morning_zenith(k)
  length(init)
  clear init
 end
  
parfor i=1:length(morning_zenith)
output1(i,:)=[morning_timestamp(i),tilt{i},power(i)];
end

clear x0 year init lb ub Y is_Shaded objective First index_list potential_tilt_angle rounded_potential_tilt_angle rounded_potential_power power tilt m n i j k y z Z fval X





%% repeat the process for evening time period

if rem(num_trackers,2)==0
    x0=[-a,a];
    x0=repmat(x0,1,num_trackers/2);
else
    x0=[a,-a];
    x0=repmat(x0,1,floor(num_trackers/2));
    x0=[x0,a];
end


% runs=50;

for k=1:length(evening_zenith)
    
    if evening_zenith(k)>50
         fminconoptions.MaxIterations=10;
     else
         fminconoptions.MaxIterations=1;
     end
    
    [~,shading_profile]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x0,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);         
    shading_profile=transpose(shading_profile)+shading_profile;  
    temp=max(shading_profile);
counter=1;
is_shaded=zeros(1,length(temp)/2);
for i=1:2:length(temp)
    if temp(i)==1 || temp(i+1)==1
        is_shaded(counter)=1;
    else
        is_shaded(counter)=0;
    end
    counter=counter+1;
end
is_shaded=is_shaded*modifier;

  lb=x0- (swing);
parfor i=1:length(lb)
    if lb(i)<min_tilt(i)
        lb(i)=min_tilt(i);
    end
end

ub=x0+ (swing-is_shaded);
parfor i=1:length(ub)
    if ub(i)>max_tilt(i)
        ub(i)=max_tilt(i);
    end
end 
  
parfor i=1:length(x0)
    if x0(i)<-4
        lb(i)=x0(i);
    end
end



 A=-4:0.1:4;
[m,n]=ndgrid(A,A);
Z=[m(:),n(:)];
init=cell(1,length(Z));
for i=1:length(Z)
    init{i}=[Z(i,1),Z(i,2)];
 if rem(num_trackers,2)==0
    init{i}=repmat(init{i},1,num_trackers/2);
 else
    init{i}=repmat(init{i},1,floor(num_trackers/2));
    init{i}=[init{i},Z(i,1)];
end
init{i}=x0+init{i};

for j=1:num_trackers
    if init{i}(j)>ub(j)
        init{i}(j)=ub(j);
    elseif init{i}(j)<lb(j)
        init{i}(j)=lb(j);
    end
end

init{i}=repmat(init{i},num_particles,1);

end

% 2461,4101,1,6561   



First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized. 

[~,is_not_smooth]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,ub,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       

is_not_smooth=max(is_not_smooth(:));
  
                          
   if is_not_smooth==0
          init={init{2461},init{4101},init{1},init{6561},init{2543},init{2625},init{2707},init{2789}};
  
                                      
    


parfor i=1:length(init)
  %options.InitialSwarmMatrix=init{i};
%[Z,z]=particleswarm(objective,NVARS,lb,ub,options);

[fval(i),~]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,init{i}(1,:),NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);         
%X{i}=Z;
%fval(i)=z;
 end
[~,y]=sort(fval,'descend');

 index_list=ones(1,4);
 for i=1:4
     index_list(i)=y(i);
 end
 if is_not_smooth==1
     index_list=index_list(1);
     index_list=repmat(index_list,1,runs);
 else
     index_list=index_list(1);
 end
 
 
     
 for n=1:length(index_list)
 Y{n}=init{index_list(n)};
 end
  end 

 if is_not_smooth==1 
index_list=1;
[potential_tilt_angle{1},~]=surrogateopt(objective,lb,ub,intcon,surrogateoptions);
 
else 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=fmincon(objective,Y{m}(1,:),[],[],[],[],lb-0.04,ub+0.04,[],fminconoptions);
end
end
 


%% applying roudning to the nearest 10th
for i=1:length(index_list)
    rounded_potential_tilt_angle{i}=round(potential_tilt_angle{i},1);
     for j=1:length(ub)
    if  rounded_potential_tilt_angle{i}(j)>potential_tilt_angle{i}(j)
        rounded_potential_tilt_angle{i}(j)= rounded_potential_tilt_angle{i}(j)-0.1;
    end
     end
 end
%% apply rounding to the nearest degree
% for i=1:length(index_list)
% for j=1:length(ub)
% rounded_potential_tilt_angle{i}(j)=floor(potential_tilt_angle{i}(j));
% end
% end

%% no rounding
% rounded_potential_tilt_angle=potential_tilt_angle;

 
for i=1:length(index_list)
 [rounded_potential_power(i),~]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,rounded_potential_tilt_angle{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
rounded_potential_power(i)=rounded_potential_power(i)/1000;
end

 [~,v2]=max(rounded_potential_power);
 tilt{k}=rounded_potential_tilt_angle{v2};
 power(k)=rounded_potential_power(v2);
 x0=tilt{k};
 
%  if runs>10
%     runs=runs-2;
%  end
length(init)
clear init
tilt{k}
  evening_zenith(k)
end
parfor i=1:length(evening_zenith)
output2(i,:)=[evening_timestamp(i),tilt{i},power(i)];
end


final=sortrows([output1;output2]);



header=('local time');
     t=('power(kW)');
   formatSpec = "Tracker: %d";
 parfor i=1:num_trackers
     temporary_cell{i}=sprintf(formatSpec,i');
      array{i}=convertStringsToChars(temporary_cell{i});
 end
 headings=[header,array,t];
 values=num2cell(final);
% final=num2cell(final);
% mirrorfinal=num2cell(mirrorfinal);
% temp=sortrows([final;mirrorfinal]);
matrix=[headings;values];
%xlswrite('fminimum.xls',matrix);
writecell(matrix,'surrogateopt.xlsx');
end
