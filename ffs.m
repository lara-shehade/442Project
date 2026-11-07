function [xhat, ck] = ffs(xt, t, n, T)
% Finite Fourier series of the chunk xt(t) with period T and 2n+1 terms:
%   ck   = (1/T) * integral over one period of x(t) e^{-j 2 pi k t / T} dt
%   xhat = sum_{k=-n}^{n} ck e^{j 2 pi k t / T}
% The period starts at t(1). If T is longer than the chunk, x(t) = 0 on
% the rest of the period; if T is shorter, only the first T seconds are used.
if iscolumn(t),  t  = t';  end
if iscolumn(xt), xt = xt'; end

tau = t - t(1);                       % time measured from the start of the period
in  = tau <= T;                       % samples inside one period
k   = (-n:n).';                       % column of harmonic indices

E  = exp(-1j*2*pi*k*tau(in)/T);       % (2n+1) x (samples in period)
ck = trapz(tau(in), E .* xt(in), 2).' / T;

xhat = real(ck * exp(1j*2*pi*k*tau/T));   % x(t) is real
end