% YANYAN WANG ON FEB.6 2023

%% 1)
sx = S.x;
sy = S.y;
dist = S.distance;
IXgrid = S.IXgrid;
chip = CHI.Z(IXgrid);
zp = double(DEM.Z(IXgrid));
zp(isnan(zp)) = 0; % this is for basins covering coastal plain, 
ordList = S.orderednanlist; 
strmBreaks = find(isnan(ordList));

%% 2) % find channel head GRID index
headsGRID = nan(length(strmBreaks),1);
id=0;
for i = 1:length(strmBreaks)    
    headsGRID(i) = ordList(id+1);
    id = strmBreaks(i);    
end
headsGRID = IXgrid(headsGRID);

%% 3) plot profiles 

for i = 1:length(strmBreaks)
    
    % extract the entire channel with the channel head
    Splot = STREAMobj(FD,'channelheads',headsGRID(i));           
    dist2 = Splot.distance;
    IXgrid2 = Splot.IXgrid;
    chip2 = CHI.Z(IXgrid2);
    zp2 =double(DEM.Z(IXgrid2));    
    zp2(isnan(zp2)) = 0; % this is for basins covering coastal plain, 
    ordList2 = Splot.orderednanlist; 
    strmBreaks2 = find(isnan(ordList2));
    id2=0;
    ii=1;
    strmInds2 = ordList2(id2+1:strmBreaks2(ii)-1);   
    % smooth elevation with window size 2000 m
    zp2(strmInds2) = smoothChannelZ(zp2(strmInds2),2000,DEM.cellsize);
    displot = dist2(strmInds2);
    zplot = zp2(strmInds2);
    chiplot = chip2(strmInds2);
    
   
    figure(i)
    % plot river profile
    subplot(2,1,1)
    hold on
    plot(displot/1000,zplot,'-','lineWidth',4,'color','k'); % here change the color 'k' into the color you want 
    xlabel('Distance (km)','FontSize',24); ylabel('Elevation (m)'); 
    set(gca,'FontSize',16)
    axis tight
    box on
    grid on
    % plot Chi profile
    subplot(2,1,2)
    plot( chiplot,zplot,'-','lineWidth',4,'color','r'); hold on % here change the color 'k' into the color you want 
    xlabel('\chi ','FontSize',24); ylabel('Elevation (m)');
    set(gca,'FontSize',16)
    box on
    grid on
    axis tight
    
    
end

