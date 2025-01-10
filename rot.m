function output=rot(input,axis,angle)
% input is the input coordinates of a point
% axis is the axis of rotation, it is a vector
% angle is the angle at whoch the coordinates are to be rotated by in
% degrees

% this function rotates coordinates along an axis by an angle in degrees

% the output is the new coordinates after rotation

angle=deg2rad(angle);
axis=axis/norm(axis);
output=input.*cos(angle) + cross(input,axis)*sin(angle) + axis.*(dot(axis,input)).*(1-cos(angle));

end

