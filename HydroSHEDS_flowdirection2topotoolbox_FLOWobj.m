
%% examine the flow direction data to see if any close point where water not
% flowing to anywhere else.

fd_WGS = GRIDobj('/Users/yanywang/Downloads/hyd_as_dir_30s/hyd_as_dir_30s.tif');

% crop DEM to Mekong watershed
mkshed = shaperead('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/mk_watershed_latlon.shp');
P = polygon2GRIDobj(fd_WGS,mkshed);
hydroFD = fd_WGS;
hydroFD.Z(~P.Z) = nan;
hydroFD = crop(hydroFD);


%% contruct flow diretion receiver and donor arrays 
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

%% construct the flow direction M matrix from the donor and reveiver array, calculate FLOWobj
nrc = numel(hydroFD.Z);
M = sparse(icc,icd,1,nrc,nrc);
cs = hydroFD.cellsize;
siz = hydroFD.size;
FD = FLOWobj(M,'cellsize',cs,'size',siz,'refmat',hydroFD.refmat);
FD.refmat = hydroFD.refmat ;
FD.georef = hydroFD.georef ;

%%
A = flowacc(FD);

%% DEM
DEMwgs = GRIDobj('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/mkhydrodem.tif');

% crop DEM to Mekong watershed
mkshed = shaperead('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/mk_watershed_latlon.shp');
P = polygon2GRIDobj(DEMwgs,mkshed);
DEM = DEMwgs;
DEM.Z(~P.Z) = nan;
DEM = crop(DEM);
DEM = resample(DEM,A);

%% save data into .mat file for feed to catchment_limited_relief function

save('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/Mekong_hydrosheds_flowdirection.mat','FD')
save('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/Mekong_hydrosheds_voidfilledDEM.mat','DEM')
save('/Users/yanywang/ESD Group Dropbox/yanyan wang/MATLAB/HydroSHEDS_flow_direction2Topotoolbox_FLOWobj/Mekong_hydrosheds_flowaccumulation.mat','A')






