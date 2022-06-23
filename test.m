function [s]=test()
    if isInParfor
        s='is in parfor';
    else
        s='not in parfor';
    end
end