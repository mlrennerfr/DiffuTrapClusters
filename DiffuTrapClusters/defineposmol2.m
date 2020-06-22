function [possynx,possyny,newposx,newposy]=defineposmol2(newposx, newposy, maxsize)
% function [possynx,possyny,newposx,newposy]=defineposmol2(newposx, newposy, maxsize)
% called by DoSimulate.m
% provides a pixel localization for comparison with synaptic image
% (possynx,possyny)
% corrects newposx,newposy considering bouncing back at the borders.
%
% Marianne Renner 11/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if newposy<0  %transposed!
    newposy=abs(newposy);
else
    if newposy>maxsize
        newposy=maxsize-(newposy-maxsize);
    end
end

possynx=floor(newposy);
if floor(newposy)<1
    possynx=1;
end
if floor(newposy)>maxsize
    possynx=maxsize; %size image;
end

if newposx<0  %transposed!
    newposx=abs(newposx);
else
    if newposx>maxsize
        newposx=maxsize-(newposx-maxsize);
    end
end

possyny=floor(newposx);
if floor(newposx)<1
    possyny=1;
end
if floor(newposx)>maxsize
    possyny=maxsize; %size image;
end

%eof%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


