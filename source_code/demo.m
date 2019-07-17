clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this script implemenats the geneartion of spherical phantom discuused in
%following paper:
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%Parameters Initialization%%%%%%%%%%%%

%%% signal parameters
N_min=1; N_max=100;
N=100; %size of the image NxNxN
grad;%diffusion-weighted gradients
b=3000;%b-value
snr=10;%SNR for the dMRI
voxel=2;%voxel size of dMRI data

%60 nodes on the surface of a sphere
load('atlas.mat');
load('dist.mat'); %distance matrix for given atlas
nodes=max(max(max(atlas)));
conn=zeros(nodes,nodes);

grid_space = linspace(N_min,N_max, N);
[Y,X,Z] = meshgrid(grid_space);
%random seed matrix (should be changed for generating different connectome)
load('seed.mat');

%%%%%%%%%%%Tubular Fiber%%%%%%%%%%%%%%%%
tic
disp('Define fibers');
%for generating a phantom with Cv=50%
fibers=330;

%using generative model for structural connectiivty generation
%Generative models described in the study by Betzel et al (2016) in Neuroimage

% set model type
modeltype = 'matching';%'neighbors';
% set whether the model is based on powerlaw or exponentials
modelvar = [{'powerlaw'},{'powerlaw'}]; 
% choose some model parameters       
params=[0.12,-0.98];
% generate synthetic networks
B = generative_model(seed,dist,fibers,modeltype,modelvar,params);
 Asynth = zeros(nodes,nodes);
 a = zeros(nodes); a(B(:,1)) = 1; a = a + a'; 
 Asynth(:,:) = a; 
index = find(triu(Asynth));
[rr1, rr2]=ind2sub([nodes,nodes], index);
    
        
%parameter for defining tubes
ncurv=[4,5,6,7,8,9];
fcurv=[4,5,6,7,8,9];
complexity=zeros(fibers);
thick=zeros(fibers);
curve=zeros(fibers);
for fiber=1:fibers

   r1=rr1(fiber);
   r2=rr2(fiber);
   thick(fiber)=1;
    if abs(r1-r2)<=8
        complexity(fiber)=1;
        choose=randsample(4,1);
        curve(fiber)=ncurv(choose);
    else
        complexity(fiber)=2;
        choose=randsample(4,1);
        curve(fiber)=fcurv(choose); 
    end
    
[x(:,fiber),y(:,fiber),z(:,fiber),conn] = create_track(X,Y,Z,conn,atlas,r1,r2,thick(fiber),curve(fiber),complexity(fiber));
end
toc


%%%%%%%%%%%%%%Vector estimation for each fiber%%%%%%%%%%%%%%
disp('Create Vectors and create fiber bundles');
tic
img=zeros(N,N,N,size(x,2));
vx=zeros(N,N,N,size(x,2));
vy=zeros(N,N,N,size(x,2));
vz=zeros(N,N,N,size(x,2));

for i=1:size(x,2)
[img(:,:,:,i),vx(:,:,:,i),vy(:,:,:,i),vz(:,:,:,i)]=create_vectors(img(:,:,:,i),X,Y,Z,x(:,i),y(:,i),z(:,i));
[img(:,:,:,i),vx(:,:,:,i),vy(:,:,:,i),vz(:,:,:,i)]=dilate_fiber(N,img(:,:,:,i),vx(:,:,:,i),vy(:,:,:,i),vz(:,:,:,i),thick(i));
end
toc

%%%%%%%%%%%%%%%%%%dMRI Signal Generation%%%%%%%%%%%%%%%%%
tic
disp('Generating dMRI signal');
[mu, mask]=create_nii(img,vx,vy,vz,b,g);
index=find(mask(:,:,:,1));
S0=mu(:,:,:,1);
n=mean(S0(index))/snr;
%Add Rice noise
%https://au.mathworks.com/matlabcentral/fileexchange/14237-rice-rician-distribution
mu=ricernd(mu,n);
toc

%convert to nifti format and save
%https://au.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
tic
disp('Converting to Nifti and Saving');
mask_nii=make_nii(mask,voxel);
save_nii(mask_nii, 'mask.nii.gz');
nii=make_nii(mu,voxel);
save_nii(nii, 'image.nii.gz');
save('conn.mat','conn') 
toc

