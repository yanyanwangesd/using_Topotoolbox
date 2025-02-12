

function [DEM, FD, A] = use_hydrosheds_dem_flowdir(dem_WGS, fd_WGS, shptxt)
% this function uses hydrosheds product directly downloaded (no interpolation at all)
% to derive the CORRECT flowpath, flow accumulation, for stream
% extraction. This function can precisely extract the flow path, and
% connections, especially in flat area. See documentation from HydroSHEDS
% website. 
%
% Inputs: dem_WGS, GRIDobj, DEM
%         fd_WGS, GRIDobj, flowdirection
%         shptxt, string, shapefile name including the file path
%
% Syntax:
%       flowdirtxt = '/Users/yanywang/Downloads/hyd_as_dir_30s/hyd_as_dir_30s.tif';
%       demtxt ='/Users/yanywang/Downloads/hyd_as_dem_30s/hyd_as_dem_30s.tif';
%       shptxt = '/Users/yanywang/basin_polygons/east_asia_basins_buffer2coverALLsamples.shp';
%       fd_WGS = GRIDobj(flowdirtxt);
%       dem_WGS = GRIDobj(demtxt);
%       [DEM, FD, A ] = use_hydrosheds_dem_flowdir(dem_WGS, fd_WGS, shptxt);
%       crita = 10000e6; % in [m2]
%       S = STREAMobj(FD,A>crita); % Create stream object
%
% 
% Outputs:
%       DEM, GRIDobj with georeference in WGS84
%       FD, FLOWojb with  georeference in WGS84, size aligned with DEM
%       A, GRIDobj with georeference in WGS84, the scalar A.Z is in unit of [m2]


% Author: Yanyan Wang (wangyanyan0607@hotmail.com)
% Last update on Feb. 12, 2025



%% step 2)  crop GRIDobj to the study area polygon
basinshed = shaperead(shptxt);
P = polygon2GRIDobj(fd_WGS,basinshed);

% flow direction
hydroFD = fd_WGS;
hydroFD.Z(~P.Z) = nan;
hydroFD = crop(hydroFD);

% dem
DEM = dem_WGS;
DEM.Z(~P.Z) = nan;
DEM = crop(DEM);


%% step 3) contruct flow diretion receiver and donor arrays from hydrosheds flow direction raster
% contruct flow diretion receiver and donor arrays
% 1) find the receivers from HydroSHEDS flowdirection
ic = 1:1:hydroFD.size(1)*hydroFD.size(2); % linear index of all data points, donors
[ici,icj] = ind2sub(hydroFD.size,ic);  % convert the linear index to subscript for donors
fdir8 = hydroFD.Z(:); % the D8 flow direction, find the documentation of definition of flowdirection from HydroSHEDS website

% initialise the reveiver array
icdi = ici;
icdj = icj;

% case 1, flow to right
id = fdir8==1;
icdi(id) = ici(id);
icdj(id) = icj(id)+1;

% case 2, flow to right bottom corner
id = fdir8==2;
icdi(id) = ici(id)+1;
icdj(id) = icj(id)+1;

% case 3,flow to bottom
id = fdir8==4;
icdi(id) = ici(id)+1;
icdj(id) = icj(id);

% case 4,flow to left bottom corner
id = fdir8==8;
icdi(id) = ici(id)+1;
icdj(id) = icj(id)-1;

% case 5, flow to left
id = fdir8==16;
icdi(id) = ici(id);
icdj(id) = icj(id)-1;

% case 6, flow to top left corner
id = fdir8==32;
icdi(id) = ici(id)-1;
icdj(id) = icj(id)-1;

% case 7, flow to top
id = fdir8==64;
icdi(id) = ici(id)-1;
icdj(id) = icj(id);

% case 8, flow to top right corner
id = fdir8==128;
icdi(id) = ici(id)-1;
icdj(id) = icj(id)+1;

% case 9, flow to itselft
id = fdir8==0;
icdi(id) = ici(id);
icdj(id) = icj(id);

% case 10, make the four boundaries flow to itself
ID = ici~=1&ici~= hydroFD.size(1)&icj~=1&icj~= hydroFD.size(2); % the non-boundaries
icdi(~ID) = ici(~ID);
icdj(~ID) = icj(~ID);

% convert to linear index
icd = sub2ind(hydroFD.size, icdi, icdj);% receiver
icc = ic;

%remove nan points
idnan = isnan(fdir8);
icc = icc(~idnan); % the donor
icd = icd(~idnan); % the receiver

clear ic

%% step 4) construct sparse matrix M to feed into FLOWobj function
% first get rid of cyclic connections
ix = icc;
ixc = icd;
cycleid = ix==ixc;
ix = ix(~cycleid);
ixc = ixc(~cycleid);
% construct M
nrc = numel(hydroFD.Z);
M = sparse(ix,ixc,1,nrc,nrc);

%% step 5) construct FLOWobj with M
cs = hydroFD.cellsize;
siz = hydroFD.size;
FD = FLOWobj(M,'cellsize',cs,'size',siz,'refmat',hydroFD.refmat);
FD.refmat = hydroFD.refmat ;
FD.georef = hydroFD.georef ;
FD.cellsize = cs;


%%  step 6) construct flow accumulation and streamobj
% first, we convert cell area from [degree2] to [m2] by constructing a
% weight array
W0 = hydroFD;
W0.Z = ones(hydroFD.size);

res = hydroFD.refmat(2);
rows = 1:(W0.size(1,1)); % same amount of rows with W0
rows = rows'; % transpose of rows
cols = 1:(W0.size(1,2));
[R,Q] = ndgrid(rows,cols); % two matrices, same size with different values
[x,y] = pix2map(hydroFD.refmat,R,Q); %calculates map coordinates X, Y from pixel coordinates rows, cols and refmat

% use areaquad.m
e = referenceEllipsoid('WGS84','m'); % return the ellipsoid object with the SemimajorAxis and SemiminorAxis properties in meters
area = areaquad(y-(res/2),x-(res/2),y+(res/2),x+(res/2),e); % convert cell area from [degree2] to [m2]
W0area = W0.*area; %weight of each cell .* its area in m^2;

% calculate flow accumulation with weight
A  = flowacc(FD,W0area); % flow accumulation in [m2]

























