function [t_sample, x_sample] = sample(t, xt, fs)
% Samples xt (defined on the dense grid t) at rate fs.
if iscolumn(t),  t  = t';  end
if iscolumn(xt), xt = xt'; end

t_sample = min(t):1/fs:max(t);
x_sample = interp1(t, xt, t_sample, 'spline');   % spline: accurate between grid points
end