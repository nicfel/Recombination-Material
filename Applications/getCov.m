function [cov] = getCov(Sequence, name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

cov = sum(Sequence=='N')/length(Sequence);

%% check how many small blocks of non N's there are to evaluate alignment quality
vals = false(length(Sequence),1);
vals(Sequence=='N') = true;

if cov>0.25

    blocks = 0;
    curr = vals(1);
    curr_length = 0;
    for i = 2 : length(vals)
        if curr~=vals(i)            
            curr=vals(i);
            if curr==1 && curr_length<100
                blocks = blocks+1;
            end
            curr_length = 0;
        end
        curr_length = curr_length+1;
    end
    if blocks > 5
        cov = 1;
    end
end

end

