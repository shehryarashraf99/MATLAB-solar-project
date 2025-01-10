clc
clear all

%% Choose the Site
[file,path] = uigetfile('*.xlsx');
site=xlsread(file); % reads the file chosen
site_specs=[site(:,1),site(:,2),site(:,3)]; % first three columns of excel file gives us tracker coordinates
 num_trackers=length(site_specs)/2; % number of trackers present in the site.
 A=site(:,4);
 B=site(:,5);
 C=site(:,6);
 D=site(:,7);
 E=site(:,8);
chord_length = A(~isnan(A))';
panels_per_tracker = B(~isnan(B))';
panel_area=C(~isnan(C))';
max_tilt=D(~isnan(D))';
min_tilt=E(~isnan(E))';

%% calculate distances and  tracker id lists for irradiance corrections

EAST_WEST=site_specs(:,1);
NORTH_SOUTH=site_specs(:,2);

 k=1;
for i=1:length(EAST_WEST)
   
    if rem(i,2)==0
        
    updated(k)=EAST_WEST(i);
     ok(k)=NORTH_SOUTH(i);
    k=k+1;
    end
    
end
east_distance_matrix=(distmat(transpose(updated),transpose(updated))).^2;% distance matrix
north_distance_matrix=(distmat(transpose(ok),transpose(ok))).^2;
distance_matrix=sqrt(east_distance_matrix+north_distance_matrix);
morning_tracker_distance=zeros(1,num_trackers);
morning_tracker_id=zeros(1,num_trackers);
evening_tracker_distance=zeros(1,num_trackers);
evening_tracker_id=zeros(1,num_trackers);



for i=1:num_trackers
        temp=distance_matrix(i,:);
        temp2=sort(temp);
        out1=temp2(2);
        out2=temp2(3);
       idx1=find(temp==out1);
       idx2=find(temp==out2);
       if length(idx1)>1
           idx1(2)=[];
       end
       if length(idx2)>1
           idx2(1)=[];
       end
        
        disp1=updated(i)-updated(idx1);
        disp2=updated(i)-updated(idx2);
        if disp1 && disp2 <0
            morning_tracker_distance(i)=out1;
            morning_tracker_id(i)=idx1;
        end
        if disp1 && disp2 >0
            evening_tracker_distance(i)=out1;
            evening_tracker_id(i)=idx1;
        end
        
        if disp1>0 && disp2<0
            evening_tracker_distance(i)=out1;
            evening_tracker_id(i)=idx1;
            morning_tracker_distance(i)=out2;
            morning_tracker_id(i)=idx2;
        end
        if disp1<0 && disp2>0
             evening_tracker_distance(i)=out2;
            evening_tracker_id(i)=idx2;
          morning_tracker_distance(i)=out1;
           morning_tracker_id(i)=idx1;
        end
        
end






%% Choose the Site Characteristics
prompt={'Enter Temperature (in degrees celsius)', 'Enter Altitude in feet', 'Enter Latitude', 'Enter Longitude', 'Enter UTC offset', 'Enter Year', 'Enter Month number', 'Enter day of the month'};
dlgtitle='Enter Site Characteristics';
dims=[1 70];
definput={'20','5180.8','35.04','-106.62','-7','2012','1','1'};
site_characteristics = inputdlg(prompt,dlgtitle,dims,definput);

temperature=str2double(site_characteristics{1});
altitude=str2double(site_characteristics{2});
latitude=str2double(site_characteristics{3});
longitude=str2double(site_characteristics{4});
UTC=str2double(site_characteristics{5});
year=str2double(site_characteristics{6});
month=str2double(site_characteristics{7});
day_of_month=str2double(site_characteristics{8});
Pressure=pvl_alt2pres(altitude);

%% enter panel Characteristics
prompt={' STC Efficiency (in decimal)', ' STC cell Temperature  (in degrees Celsius)', 'Nominal Operating Cell Temperature (in degrees Celsius)', ' STC open circuit voltage (in volts)','Power Coefficient (in negative decimal)'};
dlgtitle='Enter Panel characteristics';
dims=[1 70];
definput={'0.192','25','45','44.4','-0.36'};
panel_characteristics=inputdlg(prompt,dlgtitle,dims,definput);

STC_eff=str2double(panel_characteristics{1});
STC_cell_temp=str2double(panel_characteristics{2});
NOCT=str2double(panel_characteristics{3});
STC_OC=str2double(panel_characteristics{4});
Power_Coeff=str2double(panel_characteristics{5});

%% choose optimization strategy

list={'Genetic Algorithm', 'Pattern Search' ,'Surrogate Optimization','particleswarm','simulated anealing','fminimum','firefly','TLBO'}; % three optimizers to choose from

[choice,tf] = listdlg('PromptString',{'Select an Optimization Strategy.',...
    ''},...
    'SelectionMode','single','ListString',list);
         
 
          
 %% Choose whether to include reflection as part of the optimization         
          
third_list={'Yes', 'No'};          % check if reflection irradiance be part of the calculation.

[second_choice,~] = listdlg('PromptString',{'Include Reflection Irradiance?',...
    ''},...
    'SelectionMode','single','ListString',third_list);
if second_choice==1
    reflect=1;
else
    reflect=0;
end



 
 
 
 %% enter optimization constraints
 
 prompt={'Enter Swing limit in degrees', 'Enter Slope adjustment in percent' };
 dlgtitle= 'Enter Optimization constriants';
 dims=[1 70];
definput={'4','0'};
swing_char=inputdlg(prompt,dlgtitle,dims,definput);

swing=str2double(swing_char(1));
slope=str2double(swing_char(2));



%% caulculating irradiances, and sun coordinates through the day

Location.latitude = latitude;
Location.longitude =longitude;
Location.altitude = altitude/3.2;
DN = datenum(year, month,day_of_month):1/(24*60):datenum(year, month, day_of_month, 23, 59, 59);
Time = pvl_maketimestruct(DN,UTC);
[SunAz, SunEl, ~,~]=pvl_ephemeris(Time, Location);
[ClearSkyGHI, ClearSkyDNI, ClearSkyDHI]= pvl_clearsky_ineichen(Time, Location);
counter=1;
for i=1:length(SunEl)
    if SunEl(i)>0
    elevation(counter)=SunEl(i);
    minute(counter)=Time.minute(i);
    hour(counter)=Time.hour(i);
    azimuth(counter)=SunAz(i);
    GHI(counter)=ClearSkyGHI(i);
    DNI(counter)=ClearSkyDNI(i);
    DHI(counter)=ClearSkyDHI(i);
    counter=counter+1;
    end

end

zenith=90-elevation;




sampler=2; %update every 2 minutes
swing=swing*sampler/4;

zenith=zenith(1:sampler:end);
minute=minute(1:sampler:end);
hour=hour(1:sampler:end);
azimuth=azimuth(1:sampler:end);
GHI=GHI(1:sampler:end);
DNI=DNI(1:sampler:end);
DHI=DHI(1:sampler:end);
for i=1:length(hour)
    timestamp(i)=hour(i)*100 +minute(i);
end

if rem(length(azimuth),2)==0
    partition=length(azimuth)/2;
else
    partition=ceil(length(azimuth)/2);
end

morning_azimuth=azimuth(1:partition);
morning_zenith=zenith(1:partition);
morning_timestamp=timestamp(1:partition);
morning_GHI=GHI(1:partition);
morning_DNI=DNI(1:partition);
morning_DHI=DHI(1:partition);


evening_azimuth=azimuth(partition+1:end);
evening_zenith=zenith(partition+1:end);
evening_timestamp=timestamp(partition+1:end);
evening_GHI=GHI(partition+1:end);
evening_DNI=DNI(partition+1:end);
evening_DHI=DHI(partition+1:end);

evening_azimuth=flip(evening_azimuth);
evening_zenith=flip(evening_zenith);
evening_timestamp=flip(evening_timestamp);
evening_GHI=flip(evening_GHI);
evening_DNI=flip(evening_DNI);
evening_DHI=flip(evening_DHI);

t=datetime(year,month,day_of_month);
day_of_year=day(t,'dayofyear');



morning_slope=zeros(1,num_trackers);
for i=1:1:length(morning_slope)
     if morning_tracker_id(i)==0
         morning_slope(i)=0;
    else
    morning_slope(i)=atand((site_specs(2*morning_tracker_id(i),3)-site_specs(2*i,3))/morning_tracker_distance(i));
    end
end






evening_slope=zeros(1,num_trackers);
for i=1:1:length(evening_slope)
     if evening_tracker_id(i)==0
         evening_slope(i)=0;
    else
    evening_slope(i)=atand((site_specs(2*evening_tracker_id(i),3)-site_specs(2*i,3))/evening_tracker_distance(i));
    end
end
%% giving the front tracker the average value of the slope
avg=mean(morning_slope);
for i=1:length(morning_slope)
    if morning_tracker_id(i)==0
        morning_slope(i)=avg;
    end
end

avg=mean(evening_slope);
for i=1:length(evening_slope)
    if evening_tracker_id(i)==0
        evening_slope(i)=avg;
    end
end








 
%% running optimization

tic;

  if choice==1
% 
 [tilt,fval]=ga_opt(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);
% %the answer will be in a xls file called ga_opt.xls
% 
 elseif choice==2
 [tilt,fval]=pattern_search(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);
% %the answer is in a xls file called pattern_search.xls
%  
 elseif choice ==3
 [tilt,fval]=surrogate_opt(num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing);
% 
% 
 elseif choice ==4
     [tilt,fval]=particle(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);
       
  elseif choice ==5
     [tilt,fval]=simulated_annealing(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);
 
  elseif choice ==6
   [tilt,fval]=fminimum(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);
  
  elseif choice==7
   [tilt,fval]=firefly_opt(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);   
  
   
  elseif choice==8
   [tilt,fval]=TLBO_opt(morning_azimuth,morning_zenith,morning_timestamp,morning_GHI,morning_DNI,morning_DHI,evening_azimuth,evening_zenith,evening_timestamp,evening_GHI,evening_DNI,evening_DHI,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,day_of_month,temperature,Pressure,reflect,site_specs,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area,swing,morning_slope,evening_slope,max_tilt,min_tilt);   
   end
  
  runtime=toc;
      
 runtime=runtime/3600;    


h = msgbox(['Task Completed']);




    










%%





% output in a csv file not a text file
% make a template input file.
