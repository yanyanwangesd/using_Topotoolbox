function [divideIX, dividesxyz] = pick_divide_segment(DEM, D)

% This function allows users to hand pick the two ends of an interest
% divide segment by mouse-click on a figure which shows the DEM and the
% DIVIDEobj. The path along the interest divide is found. Function returns the 
% linear index of the points on path and the coordinates.

% This function uses the MATLAB 2021b built-in path finder of shortestpath.
% Please use the right version of MATLAB.

% Input: 1) DEM: GRIDobj of DEM of which the DIVIDEobj is calculated from
%        2) D: DIVIDEobj of the divide object, preferably in type of strahler order
%        3) Interactive user input of the two ends of an interest divide
%        segment. 

% Please read the instructions printed on Command Window to select two
% points by click on figure. 

% Output: 1)divideIX, the linear index of interest divide points that
% refers back to the DIVIDEobj
%         2)dividesxyz, the x, y, and z coordinates of divide points

% Example: DEM = GRIDobj('example.tif');
%          FD  = FLOWobj(DEM,'preprocess','c');
%          ST = STREAMobj(FD,flowacc(FD)*DEM.cellsize^2>10e6);
%          D = DIVIDEobj(FD,ST,'type','strahler');
%          [divideIX, dividesxyz] = pick_divide_segment(DEM, D)
%          figure
%          plot(D)
%          hold on
%          plot(dividesxy(:,1), dividesxy(:,2),'k-','LineWidth',3)


% Author: Yanyan Wang (wangyanyan0607@hotmail.com)
% Date: Feb. 10, 2022



% plot figure for manual pick two ends of interest divide
figure
imageschs(DEM)
hold on
plot(D,'LineWidth',2,'Color','w')
fprintf('\nClick on the divides to choose divide segment head and end point,press RETURN when you finish\n');
[x, y] = ginputc('ShowPoints', true);
close all

% snap picked ends to DIVIDEobj points, and convert to linear index
pickends = nan(2,1); % the linear index of the picked ends in DIVIDEobj
[xd, yd ] = ind2coord(D, D.IX);
for i = 1:2
    dist = (x(i)-xd).^2+(y(i)-yd).^2;
    [~, I] = min(dist);
    pickends(i) = coord2ind(D,xd(I),yd(I));
end

% construct a directed graph with the junction and endpoints in DIVIDEobj
jct = D.jct;
ep = D.ep;

% construct table of connections of divide segments
ordList = D.IX;
strmBreaks = find(isnan(ordList));
edge = nan(length(strmBreaks),5);
jctep = vertcat(jct,ep); %concatenate jct and ep

id1 = 0;
for i = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    edge(i,1) = strmInds(1);
    edge(i,2) = strmInds(end);
    % find the edges that the picked end points are in
    if sum(ismember(pickends, strmInds))
        edge(i,3) = 1;
    end
    [lia1, locb1] = ismember(edge(i,1),jctep);
    [lia2, locb2] = ismember(edge(i,2),jctep);
    if lia1
        edge(i,4) = locb1;
    end
    if lia2
        edge(i,5) = locb2;
    end
    id1 = strmBreaks(i);
end

% construct graph from junctions and endpoints
s0 = edge(:,4);
t0 = edge(:,5);
s = [s0; t0];
t = [t0;s0];
G = graph(s,t);

% find shortest path between the hand-picked two divide ends
sends= edge(~isnan(edge(:,3)),4:5);
sd = sends(1);% interest divide end 1
td = sends(2);% interest divide end 2
P = shortestpath(G, sd, td);

% find the divide pieces stored in DIVIDEobj
divideIX = []; % linear index of divide of interest
edgepair = vertcat(P(1:end-1), P(2:end));
edgepair = jctep(edgepair);

id1=0;
icount=0;
for i = 1:length(strmBreaks)
    strmInds = ordList(id1+1:strmBreaks(i)-1);
    pairends = [strmInds(1), strmInds(end)];
    
    for j = 1:length(edgepair)
        c = union(pairends, edgepair(:,j));
        if numel(c)==2 
            if pairends(1)==edgepair(1,j)
                divideIX = [divideIX; strmInds; NaN];
            else
                divideIX = [divideIX; flip(strmInds); NaN];
            end
            icount = icount+numel(strmInds);
        end
    end    
    id1 = strmBreaks(i);
    
end

% convert linear index into coordinates
[dividesx, dividesy] = ind2coord(D, divideIX); % the x,y coordinates of points on the picked divide segment. 
divideIND = coord2ind(DEM,dividesx, dividesy);

%  elevation of points on the divide segment 
dividesz = DEM.Z(divideIND); 

% summarize x, y, z into one variable
dividesxyz = [dividesx, dividesy,dividesz];
%dividesxyz = [dividesx; dividesy; dividesz];



