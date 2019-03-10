function [k,klog]=k_across(lat,lon,h);

nn=size(h,3);
h_ave=sum(h,3)./nn;

[hax,hay]=dist_gradient(h_ave,lat,lon);
abs_delta_h=hax.^2+hay.^2;

hf=zeros(length(lon),length(lat),nn);
for i=1:nn
    hf(:,:,i)=h(:,:,i)-h_ave;
    [hfx(:,:,i),hfy(:,:,i)]=dist_gradient(squeeze(hf(:,:,i)),lat,lon);
    abs_delta_hf(:,:,i)=hfx(:,:,i).^2+hfy(:,:,i).^2;
end
abs_delta_hf_ave=sum(abs_delta_hf,3)./nn;

h_son=sqrt(sum(hf.^2,3)./nn);

for i=1:length(lon)
    for j=1:length(lat)
        k(i,j)=0.32*sw_g(lat(j),0)*h_son(i,j)/abs(sw_f(lat(j)))/(1+8*abs_delta_h(i,j)/abs_delta_hf_ave(i,j));
        
        klog(i,j)=log(k(i,j));
    end
end
