function [subject]=file2subject(filename)
    subject=split(filename,"_");
    subject=subject{1};
end