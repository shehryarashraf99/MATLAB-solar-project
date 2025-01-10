function [x,X] = calculate_coordinates(num_trackers,panels_per_tracker,chord_length,zenith,azimuth,collector_tilt,site_specs)
% sun_angle is the solar zenith
% tilt is an array containing all the tilt angles of the trackers
% input is the .xls file that contains the coordinates of the trackers


% the function outputs us the panel coordinates for each tracker and its corresponding
% shadow coordinates

num_tracker_halves=num_trackers*2;

% if zenith<0
%     collector_tilt=flip(collector_tilt);
% end



%% INITIALIZE TRACKER COORDINATES




%input=xlsread('shadow validation.xlsx'); instead of reading the file.

%input=[0 8 80; 0 10 -80; 30 12 80; 30 10 -80; 60 11 80; 60 11 -80];

% transform=zeros(size(site_specs));
% transform(:,1)=site_specs(:,1);
% transform(:,2)=site_specs(:,3);
% transform(:,3)=site_specs(:,2);

transform=site_specs;

%c=80/12;


coord=cell(1,num_tracker_halves);
parfor i=1:num_tracker_halves
    coord{i}=zeros(4,3); % initializing cell array contents for faster runtime
end



for i=1:2:num_tracker_halves
    c=chord_length(round(i/2))/2; % will be an input
c=c/12;
    
    p1=transform(i,:);                      % given xls file, reading two lines, using tracker length c,I am calculating the four corner co ordinates of the tracker
    p2=transform(i+1,:);                     % given the four corner coordinates, p3 p4 are the upper corners, p5 p6 are lower corners, p1, p2 are the center corners and form the axis of rotation    
    p3=[p1(1)+c,p1(2),p1(3)];
    p4=[p2(1)+c,p2(2),p2(3)];
    p5=[p1(1)-c,p1(2),p1(3)];
    p6=[p2(1)-c,p2(2),p2(3)];
    
    coord{i}=[p5;p6;p2;p1];  % p5 p6 p2 p1 are corner co ordinates of tracker half i, and are stored in coord{i}
    coord{i+1}=[p3;p4;p2;p1]; % p3 p4 p2 p1 are corner co ordinates of tracker half i+1 and are stored in coord{i+1}
end
%%


%% Applying rotation 


%create axis of rotation for each tracker
% initializing axis cell, to preallocate memory and increase processing speed
axis =cell(1,num_trackers);
parfor i=1:num_trackers
    axis{i}=zeros(3,3); 
end
 %for each tracker I am storing axis start point, its end point, and the vector (endpoint -start point), I then store it in the cell axis
for i=1:2:num_tracker_halves
    axis{round(i/2)}(1,:)=transform(i,:); %
    axis{round(i/2)}(2,:)=transform(i+1,:);
    axis{round(i/2)}(3,:)=transform(i,:)-transform(i+1,:);
end

%rotate each tracker around its axis of rotation

%rotation of a point along an axis is done in three stages
%1) subtract it with the point of the axis closest to it
%2) perform rotation using rotation function rot
%3) add it with the same amount you subtracted

% there are  4 corners of each tracker half and there are N trackers or 2N
% tracker half, so the function runs 2N times 4 iterations




parfor i=1:num_tracker_halves
    for j=1:4
        if j==1 || j==4
            coord{i}(j,:)=coord{i}(j,:)- axis{round(i/2)}(1,:);
            coord{i}(j,:)=rot(coord{i}(j,:),axis{round(i/2)}(3,:),collector_tilt(round(i/2)));
            coord{i}(j,:)=coord{i}(j,:)+ axis{round(i/2)}(1,:);
        else
             coord{i}(j,:)=coord{i}(j,:)- axis{round(i/2)}(2,:);
            coord{i}(j,:)=rot(coord{i}(j,:),axis{round(i/2)}(3,:),collector_tilt(round(i/2)));
            coord{i}(j,:)=coord{i}(j,:)+ axis{round(i/2)}(2,:);
        end
    end
end
%% calculate the coordinates of N panels for each tracker



% we want to find the coordinates of each panel on each tracker half

X=cell(1,num_tracker_halves);

parfor i=1:num_tracker_halves
   X{i}=break_tracker(coord{i},panels_per_tracker(round(i/2)));  %break tracker will have variable numbers of panels
end

%X{i}{j} gives us the 4 corner coordinates of the jth panel on the ith
%tracker half
%%


%% initialzing sun_position,sun_angle & ground plane

p1=[-10000,-10000,-50000];
p2=[10000,-10000,-50000];
p3=[-10000,10000,-50000];
%p4=[1000,1000,-5];

normal=cross(p1-p2,p1-p3);

  
height=492126000000;
start_pos(3)=height*cosd(zenith);
start_pos(2)=height*sind(zenith)*cosd(azimuth);
start_pos(1)=height*sind(zenith)*sind(azimuth);

% [start_pos(1),start_pos(2),start_pos(3)]=sph2cart(deg2rad(azimuth),deg2rad(zenith),height);

  
  % we have initialised the sun position, such that the sun angle is of
  % angle( sun_angle) degrees
  
%% finding intersection of all shadowing lines with ground plane



   
x=cell(1,num_tracker_halves);
                                                                                                            % parfor i=1:num_tracker_halves
                                                                                                            %     x{i}=zeros(3,2*panels_per_tracker(i)+2);
                                                                                                            % end


% we have rays hitting each panel corner of each tracker half. those two
% points constitutue a line segment., we want to know the coordinates of each line
% segment, when extrapolated,
% intersects the ground plane beneath the trackers, defined by normal vector normal, and point p1

% these coordinates are stored in x.  x{i}{j} tells the coordinates that
% rays from the sun, hit the panel corners X{i}{j} and hit the ground at these 4 co ordinates 

 parfor i=1:num_tracker_halves


   x{i}=compute_intersections(normal,p1,start_pos,X{i});
    
    
    
end

 

end

