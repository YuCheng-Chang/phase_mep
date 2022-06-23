function [flag]=isPositiveInt(x)
    flag=isequal(floor(x),x)&&isscalar(x)&&all(x>0);
end