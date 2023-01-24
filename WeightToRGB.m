function [RGB]=WeightToRGB(wR, wG, wB)
% Weight to RGB color

        maxRGB=max([wR wG wB]);
        RGB=[wR wG wB]/maxRGB;

end