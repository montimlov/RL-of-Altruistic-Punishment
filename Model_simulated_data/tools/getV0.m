function [Va0, Vb0] = getV0()

Va0 = round((randi(6)-1)/5,1);
Vb0 = 1-Va0;
