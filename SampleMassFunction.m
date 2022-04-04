function [m] = SampleMassFunction(M0,sigma,N)
% sample masses out of Gaussian in log-space, centered around mass M0, with
% STD sigma (non-dimensional), total N GCs. 
% Returning masses list m

m = exp(normrnd(log(M0),sigma,[N, 1]));

end

