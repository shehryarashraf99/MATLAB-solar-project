function [tilt,power] = simulated_annealing(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,slope,max_tilt,min_tilt)
%% hyper parameters

a=0;
num_iter=35;
MaxStallIterations=1000;
runs=1; %default is 50
Initial_Mesh_Size=8;
Function_Tolerance=10^-24;
Mesh_Tolerance=10^-8;
Temperature=10;   
%temperature function= temperaturefast or temperatureboltz or temperatureexp
ReannealInterval=100;

%%




slope_adjustment=atand(slope/100); % not used for now.

if slope_adjustment>0
slope_adjustment=ceil(slope_adjustment);
else
    slope_adjustment=floor(slope_adjustment);
end



  NVARS=num_trackers;

                                      
options= optimoptions('particleswarm', ... 
                       'FunctionTolerance',1,...
                       'MaxStallIterations',1,...
                       'InitialSwarmSpan',1,...
                      'UseParallel',true,...
                      'MinNeighborsFraction',1,...    %default is 0.25
                      'SelfAdjustmentWeight',1,... %default 1.49
                      'SocialAdjustmentWeight',1,... % default 1.49                    
                      'PlotFcn','pswplotbestf',...
                       'SwarmSize',30); %default is 100 or 10*NVARS                  
  
                   
                   

                               
                               
                   
                   
                                    
if rem(num_trackers,2)==0
    x0=[-a,a];
    x0=repmat(x0,1,num_trackers/2);
else
    x0=[a,-a];
    x0=repmat(x0,1,floor(num_trackers/2));
    x0=[x0,a];
end


 for k=1:length(morning_zenith)
     
     if abs(morning_zenith(k))>=80
        modifier=3;
        
      else
         modifier=0;
     end   
     
     
     
   
lb=x0- (swing-modifier);
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

 for i=1:length(lb)
    if rem(i,2)==1
        
        P1(i)=lb(i);
        P2(i)=x0(i);
        P3(i)=lb(i)+1;
        P4(i)=x0(i);
        P5(i)=lb(i)+1;
        P6(i)=lb(i)+2;
        P7(i)=lb(i);
        P8(i)=lb(i)+2;
        P9(i)=lb(i)+2;
        P10(i)=x0(i)-1;
        P11(i)=x0(i)-1;
        P12(i)=lb(i)+1;
        P13(i)=lb(i);
        P14(i)=x0(i)-1;
        P15(i)=x0(i)-1;
        P16(i)=x0(i);
        P17(i)=x0(i);
        P18(i)=lb(i)+2;
        P19(i)=x0(i);
        P20(i)=lb(i)+1;
        P21(i)=ub(i);
        P22(i)=x0(i);
        P23(i)=lb(i);
        P24(i)=x0(i)+1;
        P25(i)=lb(i);
        P26(i)=x0(i)+2;
        
        
  
    else
        P1(i)=x0(i);
        P2(i)=lb(i);
        P3(i)=x0(i);
        P4(i)=lb(i)+1;
        P5(i)=lb(i)+2;
        P6(i)=lb(i)+1;
        P7(i)=lb(i)+2;
        P8(i)=lb(i);
        P9(i)=x0(i)-1;
        P10(i)=lb(i)+2;
        P11(i)=lb(i)+1;
        P12(i)=x0(i)-1;
        P13(i)=x0(i)-1;
        P14(i)=lb(i);
        P15(i)=x0(i);
        P16(i)=x0(i)-1;
        P17(i)=lb(i)+2;
        P18(i)=x0(i);
        P19(i)=lb(i)+1;
        P20(i)=x0(i);
        P21(i)=x0(i);
        P22(i)=ub(i);
        P23(i)=x0(i)+1;
        P24(i)=lb(i);
        P25(i)=x0(i)+2;
        P26(i)=lb(i);
        
        
        
    end

init{1}=repmat(P1,30,1); %(4,0)
init{2}=repmat(P2,30,1);%(0,4)
init{3}=repmat(P3,30,1);%(3,4)
init{4}=repmat(P4,30,1);%(4,3);
init{5}=repmat(P5,30,1);%(2,3);
init{6}=repmat(P6,30,1);%(3,2)
init{7}=repmat(P7,30,1);%(4,2)
init{8}=repmat(P8,30,1);%(2,4)
init{9}=repmat(lb,30,1);%(4,4)
init{10}=repmat(x0,30,1);%(0,0);
init{11}=repmat(ub,30,1);%(-4,-4);
init{12}=repmat(P9,30,1);%(2,1)
init{13}=repmat(P10,30,1);%(1,2)
init{14}=repmat(P11,30,1);%(1,3)
init{15}=repmat(P12,30,1);%(3,1)
init{16}=repmat(lb+1,30,1);%(3,3)
init{17}=repmat(lb+2,30,1);%(2,2);
init{18}=repmat(x0-1,30,1);%(1,1);
init{19}=repmat(P13,30,1);%(4,1);
init{20}=repmat(P14,30,1);%(1,4);
init{21}=repmat(P15,30,1);%(0,1);
init{22}=repmat(P16,30,1);%(1,0);
init{23}=repmat(P17,30,1);%(2,0);
init{24}=repmat(P18,30,1);%(0,2);
init{25}=repmat(P19,30,1);%(0,3);
init{26}=repmat(P20,30,1);%(3,0)
init{27}=repmat(P21,30,1);
init{28}=repmat(P22,30,1);
init{29}=repmat(x0+1,30,1);
init{30}=repmat(P23,30,1);
init{31}=repmat(P24,30,1);
init{32}=repmat(x0+2,30,1);
init{33}=repmat(x0+3,30,1);
init{34}=repmat(P25,30,1);
init{35}=repmat(P26,30,1);

 end   

First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized.                                     

[~,is_Shaded]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,lb,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
 
if is_Shaded==1
      patternoptions =optimoptions('patternsearch',   ...
                                          'UseParallel', true, ...
                                           'ScaleMesh','on',...'
                                            'AccelerateMesh' , false,   ...
                                           'InitialMeshSize',Initial_Mesh_Size,...
                                            'UseCompletePoll',true,... 
                                            'MeshRotate', 'on', ...
                                            'TolMesh'        , Mesh_Tolerance,    ...
                                             'Cache', 'on',...
                                             'PollOrderAlgorithm','Consecutive',...
                                          'PollMethod', 'GSSPositiveBasisNp1');
                                      
 
                                      
                                      
else
 patternoptions =optimoptions('patternsearch',   ...
                                          'UseParallel', true, ...
                                           'ScaleMesh','on',...'
                                            'AccelerateMesh' , true,   ...
                                           'InitialMeshSize',8,...
                                            'UseCompletePoll',true,... 
                                            'MeshRotate', 'off', ...
                                            'TolMesh'        , 1e-2,    ...
                                             'Cache', 'on',...
                                          'PollMethod', 'GSSPositiveBasisNp1'); 
                                      
      
                                      
    
end


 annealingoptions=optimoptions('simulannealbnd', ...
                               'FunctionTolerance',Function_Tolerance,...
                               'MaxStallIterations',MaxStallIterations,...
                               'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping},...
                               'InitialTemperature',Temperature,...
                               'TemperatureFcn','temperatureboltz',...
                               'Hybridfcn',{@patternsearch,patternoptions},...
                               'ReannealInterval',ReannealInterval);

fminconoptions=optimoptions('fmincon', ...
                               'FunctionTolerance',Function_Tolerance,...
                               'MaxFunctionEvaluations',1000,...
                               'MaxIterations',500,...
                               'StepTolerance',0.1,...
                               'OptimalityTolerance',0.05,...
                               'ConstraintTolerance',0.05,...
                               'Display','iter',...
                                'UseParallel',true,...
                               'PlotFcn','optimplotfval');
                                      




for i=1:num_iter
  options.InitialSwarmMatrix=init{i};
[Z,z]=particleswarm(objective,NVARS,lb,ub,options);
X{i}=Z;
fval(i)=z;
 end
 [~,y]=sort(fval);
 index_list=ones(1,3);
 for i=1:3
     index_list(i)=y(i);
 end
 if is_Shaded==1
      index_list=index_list(1);
     index_list=repmat(index_list,1,runs);
      else
     index_list=index_list(1);
 end
 
 if modifier~=0
index_list=index_list(1);
 end
 
 for n=1:length(index_list)
 Y{n}=init{index_list(n)};
 end
 

 
 
if is_Shaded==1 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=simulannealbnd(objective,Y{m}(1,:),lb,ub,annealingoptions);
end 
 
else 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=fmincon(objective,Y{m}(1,:),[],[],[],[],lb,ub,[],fminconoptions);
end
end
%% rounding to the nearest 10th

% for i=1:length(index_list)
%     rounded_potential_tilt_angle{i}=round(potential_tilt_angle{i},1);
%      for j=1:length(ub)
%     if  rounded_potential_tilt_angle{i}(j)<potential_tilt_angle{i}(j)
%         rounded_potential_tilt_angle{i}(j)= rounded_potential_tilt_angle{i}(j)+0.1;
%     end
%      end
% end

%% rounding to the nearest degree
% for i=1:length(index_list)
% for j=1:length(ub)
% rounded_potential_tilt_angle{i}(j)=ceil(potential_tilt_angle{i}(j));
% end
% end
    
 %% no rounding
rounded_potential_tilt_angle=potential_tilt_angle;


parfor i=1:length(index_list)
 [rounded_potential_power(i),~]=calctotalpower(morning_GHI(k),morning_DNI(k),morning_DHI(k),morning_zenith(k),morning_azimuth(k),morning_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,rounded_potential_tilt_angle{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
rounded_potential_power(i)=rounded_potential_power(i)/1000;
end

 [~,v2]=max(rounded_potential_power);
 tilt{k}=rounded_potential_tilt_angle{v2};
 power(k)=rounded_potential_power(v2);
 x0=tilt{k};
 
 
 
 tilt{k}
  morning_zenith(k)
 end
  
parfor i=1:length(morning_zenith)
output1(i,:)=[morning_timestamp(i),tilt{i},power(i)];
end

clear x0 year init lb ub Y P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16 P17 P18 P19 P20 P21 P22 P23 P24 P25 P26 is_Shaded objective First particleoptions patternoptions index_list potential_tilt_angle rounded_potential_tilt_angle rounded_potential_power power tilt m n i j k y z Z fval X







if rem(num_trackers,2)==0
    x0=[-a,a];
    x0=repmat(x0,1,num_trackers/2);
else
    x0=[a,-a];
    x0=repmat(x0,1,floor(num_trackers/2));
    x0=[x0,a];
end




for k=1:length(evening_zenith)
    if evening_zenith(k)>=80
      
    modifier=3;
    else
      modifier=0;
    end
  lb=x0- (swing);
parfor i=1:length(lb)
    if lb(i)<min_tilt(i)
        lb(i)=min_tilt(i);
    end
end

ub=x0+ (swing-modifier);
parfor i=1:length(ub)
    if ub(i)>max_tilt(i)
        ub(i)=max_tilt(i);
    end
end 
  
    

for i=1:length(ub) 

    
    if rem(i,2)==1 
        P1(i)=ub(i); 
        P2(i)=x0(i); 
        P3(i)=ub(i)-1; 
        P4(i)=x0(i); 
        P5(i)=ub(i)-1; 
        P6(i)=ub(i)-2; 
        P7(i)=ub(i); 
        P8(i)=ub(i)-2; 
        P9(i)=ub(i)-2; 
        P10(i)=x0(i)+1; 
        P11(i)=x0(i)+1; 
        P12(i)=ub(i)-1; 
        P13(i)=ub(i); 
        P14(i)=x0(i)+1; 
        P15(i)=x0(i)+1; 
        P16(i)=x0(i); 
        P17(i)=x0(i); 
        P18(i)=ub(i)-2; 
        P19(i)=x0(i); 
        P20(i)=ub(i)-1; 
        P21(i)=x0(i);
        P22(i)=lb(i);
        P23(i)=ub(i);
        P24(i)=x0(i)-1;
        P25(i)=x0(i)-2;
        P26(i)=ub(i);
        
    else 
        P1(i)=x0(i); 
        P2(i)=ub(i); 
        P3(i)=x0(i); 
        P4(i)=ub(i)-1; 
        P5(i)=ub(i)-2; 
        P6(i)=ub(i)-1; 
        P7(i)=ub(i)-2; 
        P8(i)=ub(i); 
        P9(i)=x0(i)+1; 
        P10(i)=ub(i)-2; 
        P11(i)=ub(i)-1; 
        P12(i)=x0(i)+1; 
        P13(i)=x0(i)+1; 
        P14(i)=ub(i); 
        P15(i)=x0(i); 
        P16(i)=x0(i)+1; 
        P17(i)=ub(i)-2; 
        P18(i)=x0(i); 
        P19(i)=ub(i)-1; 
        P20(i)=x0(i);
        P21(i)=lb(i);
        P22(i)=x0(i);
        P23(i)=x0(i)-1;
        P24(i)=ub(i);
        P25(i)=ub(i);
        P26(i)=x0(i)-2;
        
    end 
end
init{1}=repmat(P1,30,1); %(4,0)
init{2}=repmat(P2,30,1);%(0,4)
init{3}=repmat(P3,30,1);%(3,4)
init{4}=repmat(P4,30,1);%(4,3);
init{5}=repmat(P5,30,1);%(2,3);
init{6}=repmat(P6,30,1);%(3,2)
init{7}=repmat(P7,30,1);%(4,2)
init{8}=repmat(P8,30,1);%(2,4)
init{9}=repmat(ub,30,1);%(4,4)
init{10}=repmat(x0,30,1);%(0,0);
init{11}=repmat(lb,30,1);%(-4,-4);
init{12}=repmat(P9,30,1);%(2,1)
init{13}=repmat(P10,30,1);%(1,2)
init{14}=repmat(P11,30,1);%(1,3)
init{15}=repmat(P12,30,1);%(3,1)
init{16}=repmat(ub-1,30,1);%(3,3)
init{17}=repmat(ub-2,30,1);%(2,2);
init{18}=repmat(x0+1,30,1);%(1,1);
init{19}=repmat(P13,30,1);%(4,1);
init{20}=repmat(P14,30,1);%(1,4);
init{21}=repmat(P15,30,1);%(0,1);
init{22}=repmat(P16,30,1);%(1,0);
init{23}=repmat(P17,30,1);%(2,0);
init{24}=repmat(P18,30,1);%(0,2);
init{25}=repmat(P19,30,1);%(0,3);
init{26}=repmat(P20,30,1);%(3,0)% 
init{27}=repmat(P21,30,1);%(-4,0)
init{28}=repmat(P22,30,1);%(0,-4)
init{29}=repmat(x0-1,30,1);
init{30}=repmat(P23,30,1);%(-4,1)
init{31}=repmat(P24,30,1);%(1,-4)
init{32}=repmat(x0-2,30,1);%(-2,-2)
init{33}=repmat(x0-3,30,1);%(-3,-3)
init{34}=repmat(P25,30,1);%(-2,4)
init{35}=repmat(P26,30,1);%(4,-2)

First=@(V)V(1); % I want to optimize the first output of calctotalpower
objective=@(x) -1*First(calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,x,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)); % objective function to be minimized. 

[~,is_Shaded]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,ub,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       

if is_Shaded==1
      patternoptions =optimoptions('patternsearch',   ...
                                          'UseParallel', true, ...
                                           'ScaleMesh','on',...'
                                            'AccelerateMesh' , false,   ...
                                           'InitialMeshSize',Initial_Mesh_Size,...
                                            'UseCompletePoll',true,... 
                                            'MeshRotate', 'on', ...
                                            'TolMesh'        ,Mesh_Tolerance,    ...
                                             'Cache', 'on',...
                                             'PollOrderAlgorithm','Consecutive',...
                                          'PollMethod', 'GSSPositiveBasisNp1');
                                      
   
                                      
                                      
else
 patternoptions =optimoptions('patternsearch',   ...
                                          'UseParallel', true, ...
                                           'ScaleMesh','on',...'
                                            'AccelerateMesh' , true,   ...
                                           'InitialMeshSize',8,...
                                            'UseCompletePoll',true,... 
                                            'MeshRotate', 'off', ...
                                            'TolMesh'        , 1e-2,    ...
                                             'Cache', 'on',...
                                          'PollMethod', 'GSSPositiveBasisNp1'); 
                                      
     
 
    
end

 for i=1:num_iter
  options.InitialSwarmMatrix=init{i};
[Z,z]=particleswarm(objective,NVARS,lb,ub,options);


X{i}=Z;
fval(i)=z;
 end
[~,y]=sort(fval);

 index_list=ones(1,3);
 for i=1:3
     index_list(i)=y(i);
 end
 if is_Shaded==1
     index_list=index_list(1);
     index_list=repmat(index_list,1,runs);
 else
     index_list=index_list(1);
 end
 
  if modifier~=0
index_list=index_list(1);
 end
     
 for n=1:length(index_list)
 Y{n}=init{index_list(n)};
 end
 

 
 
if is_Shaded==1 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=simulannealbnd(objective,Y{m}(1,:),lb,ub,annealingoptions);
end 
 
else 
for m=1:length(index_list) 
%particleoptions.InitialSwarmMatrix=Y{m};
[potential_tilt_angle{m},~]=fmincon(objective,Y{m}(1,:),[],[],[],[],lb,ub,[],fminconoptions);
end
end

%% applying roudning to the nearest 10th
% for i=1:length(index_list)
%     rounded_potential_tilt_angle{i}=round(potential_tilt_angle{i},1);
%      for j=1:length(ub)
%     if  rounded_potential_tilt_angle{i}(j)>potential_tilt_angle{i}(j)
%         rounded_potential_tilt_angle{i}(j)= rounded_potential_tilt_angle{i}(j)-0.1;
%     end
%      end
% end
%% apply rounding to the nearest degree
% for i=1:length(index_list)
% for j=1:length(ub)
% rounded_potential_tilt_angle{i}(j)=floor(potential_tilt_angle{i}(j));
% end
% end

%% no rounding
rounded_potential_tilt_angle=potential_tilt_angle;

 
for i=1:length(index_list)
 [rounded_potential_power(i),~]=calctotalpower(evening_GHI(k),evening_DNI(k),evening_DHI(k),evening_zenith(k),evening_azimuth(k),evening_timestamp(k),num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day,temperature,Pressure,reflect,site_specs,rounded_potential_tilt_angle{i},NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area);       
rounded_potential_power(i)=rounded_potential_power(i)/1000;
end

 [~,v2]=max(rounded_potential_power);
 tilt{k}=rounded_potential_tilt_angle{v2};
 power(k)=rounded_potential_power(v2);
 x0=tilt{k};


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
     temp{i}=sprintf(formatSpec,i');
      array{i}=convertStringsToChars(temp{i});
 end
 ok=[header,array,t];
 temp=num2cell(final);
% final=num2cell(final);
% mirrorfinal=num2cell(mirrorfinal);
% temp=sortrows([final;mirrorfinal]);
matrix=[ok;temp];
xlswrite('simualatedannealing.xls',matrix);
    


end


% the output of the file is in pattern_output.xls  the output is in the
% form
% solar zenith, tilt1,tilt2,tilt3 power output for 3 trackers