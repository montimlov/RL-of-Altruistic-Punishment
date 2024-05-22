function [alphaT]=bias2norm(alphaS)
% given values between 0-1 make them normally distributed
alphaT=-log(0.5./alphaS - 1);