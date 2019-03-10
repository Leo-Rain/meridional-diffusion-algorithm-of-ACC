function [fax,fay] = dist_gradient(f,lat,lon)

% only for  2 dims  [lon*lat];
ndim = ndims(f);
siz = size(f);
nn=size(f,3);

loc = cell(1, ndims(f));
for k = 1:ndims(f)
    loc(k) = {1:siz(k)};
end

dlat=lat(2)-lat(1);
dlon=lon(2)-lon(1);

dist_lon=sw_dist([0,dlat],[180,180],'km')*1000;
for j=1:siz(2)
    dist_lat(j)=sw_dist([lat(j),lat(j)],[0,dlon],'km')*1000;
end

% first dimension

g  = zeros(size(f),class(f)); % case of singleton dimension
n = siz(1);
for i=1:siz(2)
    h = loc{1}(:);
    % Take forward differences on left and right edges
    if n > 1
        g(1,i) = (f(2,i) - f(1,i))/(h(2)-h(1))/dist_lat(i);
        g(n,i) = (f(n,i) - f(n-1,i))/(h(end)-h(end-1))/dist_lat(i);
    end
    
    % Take centered differences on interior points
    if n > 2
        h = (h(3:n) - h(1:n-2))*dist_lat(i);
        g(2:n-1,i) = bsxfun(@rdivide,(f(3:n,i)-f(1:n-2,i)),h);
    end    
end

fax = g;

% second dimensions and beyond
% special case 2-D matrices to support sparse matrices,
% which lack support for N-D operations including reshape
% and indexing
n = siz(2);
h = reshape(loc{2},1,[]);
g = zeros(size(f),class(f));

% Take forward differences on left and right edges
if n > 1
    g(:,1) = (f(:,2) - f(:,1))/(h(2)-h(1))/dist_lon;
    g(:,n) = (f(:,n) - f(:,n-1))/(h(end)-h(end-1))/dist_lon;
end

% Take centered differences on interior points
if n > 2
    h = (h(3:n) - h(1:n-2))*dist_lon;
    g(:,2:n-1) = bsxfun(@rdivide,(f(:,3:n)-f(:,1:n-2)),h);
end
fay= g;