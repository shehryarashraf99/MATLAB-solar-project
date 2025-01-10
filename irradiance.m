function [TotalIrradiance] = irradiance(GHI,DNI,DHI,day_number,azimuth,zenith,Pressure,collector_tilt,collector_azimuth_angle,chord_length,shading_matrix,all_tilts,front_tracker_number,front_tracker_distance,tracker_number,reflection)
% day_number is the the day of the year, i.e 1st of January is the first
% day of the year, 1st of february is the32nd day of the year

% Pressure is atmopheric pressure in Pascals
% azimuth is the solar azimuth

% zenith is the solar zenith

% collector tilt is the array that contains the tilt angles of every
% tracker

% shading matrix is the matrix N X N matrix where N is the total number of
% tracker halves it tells us which tracker half shades which tracker half,
% if shading matrix (1,5) =1, it means tracker half 1 is shading tracker
% half 5.
% if shading_matrix(2,4) =0, then tracker half 2 is NOT shading tracker half
% 4

% tracker number is the tracker half number in question. This function runs
% in a loop N times, where N is the number of tracker halves

% reflection is the user choice. if reflection==1, reflection irradiance is
% part of the optimization, else it is not.






%%ensure all angles are greater than zero  formula for calculation of beam, diffuse and reflected irradiance require positive angles. will not work

if zenith <0
    zenith=-zenith;
end
elevation=90-zenith;


% a tracker/collector facing west with negative tilt is the same as a tracker with positive tilt but facing east.

% if collector_tilt>0
%     collector_azimuth_angle=270; % if tracker tilt angle is positive, it is facing west
% else
%     collector_azimuth_angle=90; % else it is pointing eastwards
% end

if collector_tilt<0
    collector_tilt=-collector_tilt;      
end

if collector_tilt==0
    collector_azimuth_angle=180;
end

%%


%% we assume standard pressure, that can  be adjusted
%Pressure= 1.01325*10^5; % Pascals
%%

%%  calculate beam and diffuse irradiance


 AM=pvl_relativeairmass(zenith);

 AM=pvl_absoluteairmass(AM, Pressure);
Ibh=DNI;
Idh=DHI;
if GHI<0
    GHI=0;
end
if Ibh<0
    Ibh=0;
end
if Idh<0
    Idh=0;
end

 %%

Beam_Irradiance=Ibh*(cosd(elevation)*sind(collector_tilt)*cosd(collector_azimuth_angle-azimuth)+(sind(elevation)*cosd(collector_tilt)));

if Beam_Irradiance<0  % irradiance cannot be a negative quantity.
    Beam_Irradiance=0;
end




%% incorporate shading

shade=max(shading_matrix(:,tracker_number)); % receiving shading matrix from compute_shading matrix function

if shade>0
    Beam_Irradiance=0;
end
%%


%%  calculate diffuse irradiance.
Hextra=pvl_extraradiation(day_number); % calculate extraterrestial irradiation in watts per meter square

model='1990';
%model='phoenix1988';
%model='usacomposite1988';
%model='sandiacomposite1988';
%model= 'france1988';

perez=pvl_perez(collector_tilt,collector_azimuth_angle,Idh,Ibh,Hextra,zenith,azimuth,AM,model);


%%


%% incorporate reflected irradiance

if reflection==1
ref=pvl_grounddiffuse(collector_tilt,GHI,0.2);
else
ref=0;
end

%%
%perez=0; %setting diffuse irradiance to zero because why not?


%% adjusting diffuse irradiance, beam irradiance and reflected irradiance based on tracker number and the number of the tracker in front of it
 temp=front_tracker_number(round(tracker_number/2));
  b0=0.05;
   rho=1- (b0*(secd(pvl_getaoi(collector_tilt, collector_azimuth_angle, zenith,azimuth))-1));
   if rho<0 
      rho=0;
   end
   if rho>1
       rho=0;
   end
   
 if temp==0
    diffuse_correction=0.5*(1+cosd(collector_tilt));
    reflected_correction=0.5*(1-cosd(collector_tilt));
    perez=perez*diffuse_correction;
    ref=ref*reflected_correction;
 else
 front_tilt=all_tilts(temp);
 if front_tilt<0
     front_tilt=-front_tilt;
 end

 L=chord_length(round(tracker_number/2))/12; % it is 2*chord length of the front tracker and should be read from tracker id. hard coded for testing
 R=front_tracker_distance(round(tracker_number/2)); % it is distance between tracker and the front tracker, hard coded as a constant
 Y=sqrt(L^2+R^2 -2*L*R*cosd(front_tilt));
 gamma=asind(L*sind(front_tilt)/Y);
 
 %theta=pvl_getaoi(collector_tilt, collector_azimuth_angle, zenith,azimuth);
 D0=(1+cosd(collector_tilt+gamma))/2;
 %gamma=gamma*pi/180; %% convert to radians for integration matlab cannot do integration in degrees
%  syms phi theta
%  f(phi,theta)=(sin(phi))^2 *sin(theta+collector_tilt*pi/180);
%  p1=int(f,phi,[0 3*pi/2]);
%  p2=int(f,phi,[pi/2 pi]);
%  p3=vpaintegral(p1+p2,theta,[0 gamma],'AbsTol',1e-5,'RelTol',1e-4);
%  p3=double(p3);
%  P=p3/(2*pi);
%  D=D0+P;

D=0.5*(cosd(collector_tilt)-cosd(gamma+collector_tilt)); 
D=D+D0;
 diffuse_correction=D0/D;
 
 
 
 
  %% adjusting  irradiance due to diffuse irradiance
 R=1/4 *(1- cosd(180-front_tilt-gamma+collector_tilt));
 R0=(1-cosd(collector_tilt))/2;
 reflected_correction=R0/R;
 
 
 
 
 perez=perez*diffuse_correction;
 ref=ref*reflected_correction;
end
  %% adjusting irradiances in presence of reflection


% Diam=rho*D;
%Riam=rho*R;
%ref=ref*Riam/R0

%
   

TotalIrradiance=Beam_Irradiance+perez+ ref; % output of this function

if reflection==1
TotalIrradiance=TotalIrradiance*rho;
end


end






