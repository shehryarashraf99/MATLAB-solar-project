function [tilt,power] = TLBO_opt(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt)


  NVARS=num_trackers;


%% hyper parameters
 modifier=swing-0.1;
a=0;
%num_iter=71;
runs=1; %default is 1
Initial_Mesh_Size=8;
%Function_Tolerance=10^-12; %10^-6
Mesh_Tolerance=0.05; %10^-12
%Self_adjustment_weight=0.25; % default is 0.25
%Social_Adjustment_weight=0.205; %default is 0.205
%Max_Stall_iterations=1000;%default is NVARS
%Max_iterations=100;
StepTolerance=0.05;
%InertiaRange=[0.00004 0.00009]; % previous is [0.04 0.09]
%min_neighbours=0.9;
population=35;
num_iterations=1000;


%num_function_evalutions =num_fireflies^2 *num_iterations so keep
%num_fireflies at max 10


% std=0.5/3;
% mean=0;



%%







  patternoptions =optimoptions('patternsearch',   ...
                                          'UseParallel', true, ...
                                          'StepTolerance',StepTolerance,...
                                           'ScaleMesh','on',...'
                                            'AccelerateMesh' , false,   ...
                                           'InitialMeshSize',Initial_Mesh_Size,...
                                            'UseCompletePoll',true,... 
                                            'MeshRotate', 'on', ...
                                            'TolMesh'        , Mesh_Tolerance,    ...
                                             'Cache', 'on',...
                                             'PollOrderAlgorithm','Success',...
                                          'PollMethod', 'GSSPositiveBasis2N',...
                                            'Display'        , 'Iter',  ...
                                          'PlotFcn'        ,{@psplotbestf,@psplotmeshsize});
                                          
                                     
  

fminconoptions=optimoptions('fmincon', ...
                               'FunctionTolerance',10^-24,...
                               'MaxFunctionEvaluations',10000,...
                               'StepTolerance',0.05,...
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
      
%      if morning_zenith(k)>85
%          particleoptions.MaxIterations=10;
%          particleoptions.MaxStallIterations=10;
%      else
%          particleoptions.MaxIterations=Max_iterations;
%           particleoptions.MaxStallIterations=Max_Stall_iterations;
%      end
      
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

A=-swing:0.1:swing;
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

%init{i}=repmat(init{i},population,1);

end

% 2461,4101,1,6561
        
     
if morning_zenith(k)>88
    iter=10;
else
    iter=num_iterations;
end
          

First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized.                                     

%is_not_smooth tells us if we can expect any shading when the trackers have
%been tilted up by value swing. this tells us if we can expect any non_smoothness (discountinuities) in the function,
% if is_not_smooth==0 we can speed up the optimization by reducing the
% number of datapoints we are sampling from 6561 to 4. if we are using a non
% smooth optimizer, we no longer need to use it and can use a gradient
% based optimizer instead to speed up the process without a loss in power
% output.



[~,is_not_smooth]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,lb,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
 
is_not_smooth=max(is_not_smooth(:));
    
  
   if is_not_smooth==0
clear init;
init{1}=x0+lb;
init{2}=x0-1;
init{3}=x0+1;
init{4}=x0+ub;
init{5}=x0+0.9;
init{6}=x0+0.8;
init{7}=x0+0.7;
init{8}=x0+0.6;
init{9}=x0+0.5;
init{10}=x0+0.4;
init{11}=x0+0.3;
init{12}=x0+0.2;
init{13}=x0+0.1;
      % init={init{2461},init{4101},init{1},init{6561},init{4019},init{3937},init{3855},init{3773}};
   end
    



parfor i=1:length(init)
  %options.InitialSwarmMatrix=init{i};
  [fval(i),~]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,init{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
%[Z,z]=particleswarm(objective,NVARS,lb,ub,options);
%X{i}=Z;
%fval(i)=z;
 end
 [~,y]=sort(fval,'descend');
 

   
 Y{1}=init{y(1)}; %best power output
 
 
 
 
  
  
  
 if is_not_smooth>=0.5 
for m=1:runs
    
   for i=1:population
     TLBO_matrix(i,:)=init{y(i)};
  end
 
[TLBO_tilt_angle{m},~,~,~]=TLBO(TLBO_matrix,objective,lb,ub,iter,population);
[potential_tilt_angle{m},~]=patternsearch(objective,TLBO_tilt_angle{m}(1,:),[],[],[],[],lb,ub,[],patternoptions);
end 
 
 else 

%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{1},~]=fmincon(objective,Y{1}(1,:),[],[],[],[],lb-0.01,ub+0.01,[],fminconoptions);

end


 


%% rounding to the nearest 10th

% many algorithms such as particlewarm, patternsearch, fmincon do not have
% integer contraints. as the accuracy of tracker tilt angles is 0.1 degrees
% we round the answer to the nearest 10th of a degree. also to avoid
% shading we round down. 3.99 degrees is rounded to 3.9 degrees (not 4)
% because of the possibility that shading occurs at 4.0 but if there is no
% shading at 3.99 degrees then there is no shading at 3.9 degrees.

 for i=1:length(potential_tilt_angle)
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

 %% getting the best answer from the the rounded candidates
for i=1:length(rounded_potential_tilt_angle)
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


%% repeat the process for the evening times


 %% initializing tracker tilt angles at the start of optimization. generally best if initialized to 0  

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
     
%      if evening_zenith(k)>85
%          particleoptions.MaxIterations=10;
%          particleoptions.MaxStallIterations=10;
%      else
%          particleoptions.MaxIterations=Max_iterations;
%           particleoptions.MaxStallIterations=Max_Stall_iterations;
%      end
    
    
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

ub=x0+ (swing-is_shaded); % lower bounds of shaded/shading trackers now modified 
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


 A=-swing:0.1:swing;
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

%init{i}=repmat(init{i},population,1);

end

% 2461,4101,1,6561   


if evening_zenith(k)>88
    iter=10;
else
    iter=num_iterations;
end


First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized. 

[~,is_not_smooth]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,ub,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       

is_not_smooth=max(is_not_smooth(:));
  
                          
   if is_not_smooth==0
init{1}=x0-1;
init{2}=x0+1;
init{3}=x0+lb;
init{4}=x0-0.9;
init{5}=x0-0.8;
init{6}=x0-0.7;
init{7}=x0-0.6;
init{8}=x0+ub;
init{9}=x0-0.5;
init{10}=x0-0.4;
init{11}=x0-0.3;
init{12}=x0-0.2;
init{13}=x0-0.1;
       
          %init={init{2461},init{4101},init{1},init{6561},init{2543},init{2625},init{2707},init{2789}};
   end 
                                      
    


parfor i=1:length(init)
  %options.InitialSwarmMatrix=init{i};
%[Z,z]=particleswarm(objective,NVARS,lb,ub,options);

[fval(i),~]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,init{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);         
%X{i}=Z;
%fval(i)=z;
 end
[~,y]=sort(fval,'descend');




 Y{1}=init{y(1)}; %best initial power output
 
 if is_not_smooth>=0.5 
for m=1:runs
    
 for i=1:population
     TLBO_matrix(i,:)=init{y(i)};
 end
[TLBO_tilt_angle{m},~,~,~]=TLBO(TLBO_matrix,objective,lb,ub,iter,population);
[potential_tilt_angle{m},~]=patternsearch(objective,TLBO_tilt_angle{m}(1,:),[],[],[],[],lb,ub,[],patternoptions);
end 
 
 else 

%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{1},~]=fmincon(objective,Y{1},[],[],[],[],lb-0.01,ub+0.01,[],fminconoptions);

end




%% applying rounding to the nearest 10th
 for i=1:length(potential_tilt_angle)
     rounded_potential_tilt_angle{i}=round(potential_tilt_angle{i},1);
     for j=1:length(ub)
    if  rounded_potential_tilt_angle{i}(j)<potential_tilt_angle{i}(j)
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

 %% getting the best answer from the the rounded candidates
for i=1:length(rounded_potential_tilt_angle)
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
writecell(matrix,'TLBO.xlsx');
end

