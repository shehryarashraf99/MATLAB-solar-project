function [K] = break_tracker(tracker,N)   
% this function returns us the co ordinates of N panels, belonging to
% tracker called tracker



%% compute the four corners of tracker
temp1=tracker(1,:);
temp1=transpose(temp1);

temp2=tracker(2,:);
temp2=transpose(temp2);

temp3=tracker(3,:);
temp3 =transpose(temp3);

temp4=tracker(4,:);
temp4=transpose(temp4);

%% compute all unique points on the tracker, these points form panel coordinates
 

 
 A1=linspace(temp1(1),temp2(1),N+1);
 A2=linspace(temp1(2),temp2(2),N+1);
 A3=linspace(temp1(3),temp2(3),N+1);
 A=[A1;A2;A3];
 
  B1=linspace(temp3(1),temp4(1),N+1);
 B2=linspace(temp3(2),temp4(2),N+1);
 B3=linspace(temp3(3),temp4(3),N+1);
 B=[B1;B2;B3];
 
K=[A,B];


     % output corner coordinates
     
  
     

end

