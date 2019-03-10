clc; clear all; close all;
nameall=dir('K:\allsat_adt\7_2010\');
nn=size(nameall,1);
lat_all=ncread(['K:\allsat_adt\7_2010\',nameall(3).name],'NbLatitudes');
lon_all=ncread(['K:\allsat_adt\7_2010\',nameall(3).name],'NbLongitudes');
lat=lat_all(199:364);
lon=lon_all;

h=zeros(length(lat),length(lon),nn-2);
for i=3:nn
    h(:,:,i-2)=ncread(['K:\allsat_adt\7_2010\',nameall(i).name],'Grid_0001',[199,1],[166,1080]);
end
h=permute(h,[2,1,3]);
    [k,klog]=k_across(lat,lon,h);
    nn=size(h,3);
    h_ave=sum(h,3)./nn;



figure(1)
set(gcf,'color','w')
m_proj('miller','lat',[-65 -30],'long',[0 360]);
m_grid('linewi',2,'tickdir','in');
hold on
m_coast('patch',[0.72 0.72 0.72],'edgecolor','none');
m_contourf(lon,lat,klog','linestyle','none');
caxis([2,11])
colormap('jet')
colorbar('southoutside')
hold on
c=m_contour(lon,lat,h_ave',8,'k','linewidth',0.5);
title('ACC SSH K\_across-acc','fontsize',9,'fontname','Arial')
print('-dtiff','-r800','ACC_SSH_K_across-acc');