function [output,shading_matrix] = calctotalpower(GHI,DNI,DHI,zenith,azimuth,timestamp,num_trackers,panels_per_tracker,chord_length,morning_tracker_id,evening_tracker_id,morning_tracker_distance,evening_tracker_distance,morning_slope,evening_slope,day,temperature,Pressure,reflect,site_specs,collector_tilt,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)
% output is the power output of the tracker half

%zenith is the solar zenith in degrees

% azimuth is the solar azimuth in degrees

% num_trackers is the number of trackers in the site

% temperature is the ambient temperature of the site in celsius

% collector_tilt is the an array containing the tilt angle of each tracker
% in degrees

% site specs is the .xls file of the site given as an input to set
%up the trackers

%chord_length is the length of each tracker. in feet

% panels per tracker is the number of panels per tracker


% reflect is binary, 1 means reflected irradiance is included in the
% calculation, 0 means that it is not included in the calculation


% calctotalpower is main objective function that is to be maximized in
%  with an output in Watts



%% check if it is morning or evening

if timestamp>1200
    front_tracker_number=evening_tracker_id; % to be read from excel file in future iteration
    front_tracker_distance=evening_tracker_distance;
    slope=evening_slope;
   
                                                                       
else
    front_tracker_number=morning_tracker_id; % to be read from excel file in future iteration
    front_tracker_distance=morning_tracker_distance;
    slope=morning_slope;
end




%% first compute shading matrix
% if zenith<0
%     collector_tilt=flip(collector_tilt);
% end
%     


[x,X]=calculate_coordinates(num_trackers,panels_per_tracker,chord_length,zenith,azimuth,collector_tilt,site_specs); % returns us all panel coordinates and its shadow coordinates
shading_matrix=return_shading_matrix(x,X,num_trackers,timestamp,slope,zenith); % returns us the tracker-tracker shading matrix


%% now compute irradiance and power output of each tracker half, and proceed to sum it.


 % front tracker number gives us the order in terms of west-east during the evening, and east-west during the morning the positions of the
    % trackers and tells us the tracker number of the tracker that reduces
    % the reflected and  diffuse irradiance of the tracker that is behind it due to its tilt angle. this is
    % necessary for irradiance correction as done in Helioscope documentation
    % this is done to ensure irradiance calculations are as accurate as it
    % should be in a practical scenario in a commercial pv system.
    
    


parfor i=1:num_trackers*2
    collector_azimuth_angle(i)=get_azimuth(X{i},collector_tilt(round(i/2)));
    irr(i)=irradiance(GHI,DNI,DHI,day,azimuth,zenith,Pressure,collector_tilt(round(i/2)),collector_azimuth_angle(i),chord_length,shading_matrix,collector_tilt,front_tracker_number,front_tracker_distance,i,reflect);
    %power(i)=getpower(irr(i),temperature,panels_per_tracker(round(i/2)),chord_length(round(i/2)));
    power(i)=getpower(irr(i),temperature,panels_per_tracker(round(i/2)),NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area(round(i/2)));
end
output=sum(power);


end

