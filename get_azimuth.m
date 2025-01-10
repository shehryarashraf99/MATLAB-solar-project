function [azimuth] = get_azimuth(X,collector_tilt)
% we only need 1 panel of a tracker to get the collector azimuth angle. it
% can be any panel, but we take first one in our calculations
p1=X(:,1);
p2=X(:,2);
p3=X(:,length(X));

p1=transpose(p1);
p2=transpose(p2);
p3=transpose(p3);
% we need 3 co ordinates to calculate normal vector of the panel. the
% normal vector of the panel is the same as the normal vector of the
% tracker

n1=cross(p1-p2,p1-p3);
n2=[0,0,1]; % normal vector of the earth's surface
n3=dot(n1,n2);
n4=n2.*n3;
n5=n1-n4;
n6=[0 1 0];
angle=atan2(norm(cross(n5,n6)),dot(n5,n6));
azimuth=angle*180/pi;

if collector_tilt>0
    azimuth=360-azimuth;
end

end

