function [img,vx,vy,vz]=dilate_fiber(N,img,vx,vy,vz,thick)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cross?sectional dilation of fiber trajectory
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngh=[];
for i=-1:1
    for j=-1:1
        for k=-1:1
        ngh=[ngh; i,j,k];
        end
    end
end

for th=1:thick
ind=find(img); 
for i=1:length(ind)
    [ii,jj,kk]=ind2sub([N,N,N],ind(i));
    ngh_new=ngh+repmat([ii,jj,kk],size(ngh,1),1); 
    if(ngh_new>0 & ngh_new<=N)
    ind_new=sub2ind([N,N,N],ngh_new(:,1),ngh_new(:,2),ngh_new(:,3));
    ind_tmp=find(~img(ind_new));
    for j=1:length(ind_tmp)
       

    vx(ind_new(ind_tmp(j)))=vx(ii,jj,kk);
    vy(ind_new(ind_tmp(j)))=vy(ii,jj,kk);
    vz(ind_new(ind_tmp(j)))=vz(ii,jj,kk);
    img(ind_new(ind_tmp(j)))=1; 
    end
    end
end
end