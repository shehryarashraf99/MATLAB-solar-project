function output = return_shading_matrix(x,X,num_trackers,timestamp,slope,zenith)
% this function takes in coordinates of each panel of each tracker, as well
% as its shadow coordinates, and outputs the resultant N x N shading
% matrix, where N is the number of tracker halves

elevation=90-zenith;

num_tracker_halves=num_trackers*2;

output=zeros(num_tracker_halves,num_tracker_halves); 

 % we have all the panel coordinates, and we know the co ordinates where
 % sun rays hit the panel and hits the ground. this gives us 4 lines  all
 % we have to do now is
 % check if that line intersects any other panel of any other tracker and
 % that is done below
 

%% check if tracker half i is shading tracker half k


parfor i=1:num_tracker_halves
    for k=1:num_tracker_halves
        east_west=X{i}(1,1)-X{k}(1,1);
        height=max(X{i}(3,:))-min(X{k}(3,:));
         if i==k || round(i/2)==round(k/2) 
             output(i,k)=0;
        elseif east_west<0 && timestamp<1200
             output(i,k)=0;
       elseif east_west>0 && timestamp>1200
             output(i,k)=0;
       elseif height<=0 
           output(i,k)=0;
        else
            
            output(i,k)=compute(X{k},x{i},X{i});
        end
    end
end    

 for i=1:num_tracker_halves
     if slope(round(i/2))>elevation
         output(i,i)=0.5;
     end
 end


end

