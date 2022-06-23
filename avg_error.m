function [err]=avg_error(est_phases,true_phases)
%     fprintf('est_phases=');
%     disp(size(est_phases));
%     fprintf('true_phases=');
%     disp(size(true_phases));
    assert(isequal(size(est_phases),size(true_phases)),'input arrays must have same size...');
    err=abs(angle(exp(1i*true_phases)./exp(1i*est_phases)));
    err=mean(err);
end