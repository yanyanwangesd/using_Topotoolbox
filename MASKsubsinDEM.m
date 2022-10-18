function Sref = MASKsubsinDEM(DEM,MASK)
% this function find the subscripts of non-nans in the MASK(logical) in DEM.
% MASK is logical with size the same with DEM.size, i.e., size(MASK) =
% DEM.size.

% Input: DEM, GRIDobj from Topotoolbox
%        MASK, a logical matrix that size(MASK) = DEM.size
% outputs: Sref is a struct where the subscripts are stored in the Sref.subs


% Usage: in situations related with or similar to using crop(DEM)

% Author: Yanyan Wang (wangyanyan0607@hotmail.com)
% Last modification date: Oct. 17, 2022


IX = find(MASK);
siz  = DEM.size;
k    = [1 cumprod(siz(1:end-1))];
% preallocate subsref structure
Sref    = substruct('()',cell(1,2));
% subset size
sizout = zeros(1,2);
% loop through dimensions (see ind2sub)
% and get subscripts of minimum bounding rectangle/box/...
for r = 2:-1:1
    IX2       = rem(IX-1,k(r))+1;
    subdim    = (IX-IX2)/k(r)+1;
    Sref.subs{r} = min(subdim):max(subdim);
    sizout(r) = numel(Sref.subs{r});
    IX        = IX2;
end


