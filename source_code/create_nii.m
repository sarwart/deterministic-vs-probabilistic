function [mu msk]=create_nii(img,vx,vy,vz,b,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate dMRI signal using compartment model
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for modeling fiber bundles
mu0=50;
constant=5;
fu=0.2;

lambda1=1.5*10^-3; 
lambda2=0.2*10^-3; 

diff=0.9*10^-3; 


%compartment model
beta=lambda2;
alpha=lambda1-beta;
[r c z f]=size(img);
n=size(g,2);
mu=ones(r,c,z,n)*mu0*constant*exp(-b*diff);
msk=zeros(r,c,z);

 for i=1:f
 ind=find(img(:,:,:,i));
 
 for j=1:length(ind)
 
     [ii,jj, kk]=ind2sub([r,c,z],ind(j));
     v=[vx(ii,jj,kk,i);vy(ii,jj,kk,i);vz(ii,jj,kk,i)];
     D=(alpha*v*v')+(beta*[1 0 0; 0 1 0; 0 0 1]);
     mu(ii,jj,kk,:)=squeeze(mu(ii,jj,kk,:))+diag((exp(-b*g'*D*g)));
     
 end
 
 msk=msk+img(:,:,:,i);
 
 end
 
mul=msk;
check=find(msk);
for i=1:length(check)
[x,y,z] = ind2sub(size(msk),check(i));
mu(x,y,z,:)=mu(x,y,z,:)-mu0*constant*exp(-b*diff); 
mu(x,y,z,:)=mu0*constant*(((1-fu)*mu(x,y,z,:))+(fu*exp(-b*diff)));

end

mu(:,:,:,1)=mu0*constant;
msk(find(msk>0))=1;
 
 
