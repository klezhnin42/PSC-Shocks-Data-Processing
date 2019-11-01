% bins data into histogram h
% similar to histc, but operating along both dimensions.
%
% datax, datay is a pair of vectors representing x and y data.  They
%  must be the same length (they are treated as 1-d data no matter the
%  shape)
% xedge and yedge are edges of histogram bins.  See histc for definition

function hh = histc2d(xdata, ydata, xedge, yedge)

if ( numel(xdata) ~= numel(ydata) )
    error('Expect length datax = length datay')
end

if numel(yedge) > numel(xedge)
    hh = histc2d(ydata, xdata, yedge, xedge)';
    return 
end


[sxdata,ind] = sort(xdata(:));
sydata = ydata(ind);

hx = histc(sxdata, [-inf, xedge]);
binx = 1+cumsum(hx);

hh = zeros(length(xedge), length(yedge));
for m=1:length(xedge)
    hh(m,:) = histc(sydata(binx(m):(binx(m+1)-1)), yedge);
end