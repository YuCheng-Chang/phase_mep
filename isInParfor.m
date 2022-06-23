function [flag]=isInParfor()
    flag=~ isempty(getCurrentTask());
end