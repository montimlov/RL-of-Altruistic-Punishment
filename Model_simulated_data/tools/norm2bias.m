function [ alpha_out ] = norm2bias(alpha)
% transformation from gaussian space to alpha space, as suggested in Daw 2009 Tutorial
% MKW October 2017

alpha_out = 0.5*logsig(alpha);

end

