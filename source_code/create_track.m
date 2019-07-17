function [x,y,z,conn] = create_track(X,Y,Z,conn,atlas,R1,R2,thick, dist,complexity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defines the tracjectory of fibers
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


index_R1=find(atlas==R1);
index_R2=find(atlas==R2);

ngh=[];

for i=-thick-2:1:thick+2
    for j=-thick-2:1:thick+2
        for k=-thick-2:1:thick+2
        
        temp=[i j k];
        ngh=[ngh; temp];
    
        end
    end
end
 

startp=[];
endp=[];

while isempty(startp)

    i = randsample(length(index_R1),1);
    [ii,jj, kk]=ind2sub(size(atlas),index_R1(i));
    ngh_new=ngh+repmat([ii,jj,kk],size(ngh,1),1); 
    
       if(  ngh_new>0 & ngh_new<=size(atlas,1))
           ind_new=sub2ind(size(atlas),ngh_new(:,1),ngh_new(:,2),ngh_new(:,3));
           if (atlas(ind_new)==0 | atlas(ind_new)==R1)
            startp=[ii jj kk];
   
    break;
           end
       end
    
end

while isempty(endp)

       
i = randsample(length(index_R2),1);
    [ii,jj,kk]=ind2sub(size(atlas),index_R2(i));
       
   
    ngh_new=ngh+repmat([ii,jj,kk],size(ngh,1),1); 
     
       if(ngh_new>0 & ngh_new<=size(atlas,1) )
           ind_new=sub2ind(size(atlas),ngh_new(:,1),ngh_new(:,2),ngh_new(:,3));
           
           if (atlas(ind_new)==0 | atlas(ind_new)==R2)
     
           endp=[ii jj kk];
    
            break;
           end
       end
    
end

midx = floor((startp(1)+endp(1))/2);
midy =floor((startp(2)+endp(2))/2);
midz =floor((startp(3)+endp(3))/2);

fmidx=floor((startp(1)+midx)/2);
fmidy=floor((startp(2)+midy)/2);
fmidz=floor((startp(3)+midz)/2);

smidx=floor((endp(1)+midx)/2);
smidy=floor((endp(2)+midy)/2);
smidz=floor((endp(3)+midz)/2);


neigh=[];
for i=-1:1:1
    for j=-1:1:1
        for k=-1:1:1
        
        if ~(i==0 && j==0 && k==0)
                neigh=[neigh;[i,j,k]]; 
         end
    
        end
    end
end

   temp_mid=dist.*neigh+[midx midy midz];
   
   distances=ones(length(neigh),1)*1000;
    
    for i=1:length(temp_mid)
    
    distances(i)=pdist([50 50 50;temp_mid(i,:)]);
    
    end
    [val ind]=min(distances);
    midx=temp_mid(ind,1);
    midy=temp_mid(ind,2);
    midz=temp_mid(ind,3);
    
if complexity==1 
    
Xa=midx+dist;
Ya=midy+dist;
Za=midz+dist;

elseif complexity==2
    
pmidx=floor((fmidx+midx)/2);
pmidy=floor((fmidy+midy)/2);
pmidz=floor((fmidz+midz)/2);

pmidx2=floor((smidx+midx)/2);
pmidy2=floor((smidy+midy)/2);
pmidz2=floor((smidz+midz)/2);

temp_mid=dist.*neigh+[pmidx pmidy pmidz];
temp_mid1=dist.*neigh-[pmidx pmidy pmidz];
   
   distances=ones(length(neigh),1)*1000;
    
    for i=1:length(temp_mid)
    
    distances(i)=pdist([50 50 50;temp_mid(i,:)]);
    
    end
    [val1 ind1]=min(distances);
    d1x=temp_mid(ind,1);
    d1y=temp_mid(ind,2);
    d1z=temp_mid(ind,3);
    
    
   distances=ones(length(neigh),1)*1000;
    
    for i=1:length(temp_mid1)
    
    distances(i)=pdist([50 50 50;temp_mid1(i,:)]);
    
    end
    [val2 ind2]=min(distances);
    d2x=temp_mid(ind,1);
    d2y=temp_mid(ind,2);
    d2z=temp_mid(ind,3);    

if val2<val1
Xa(1)=d2x;
Ya(1)=d2y;
Za(1)=d2z;
else
Xa(1)=d1x;
Ya(1)=d1y;
Za(1)=d1z;
end

temp_mid=dist.*neigh+[pmidx2 pmidy2 pmidz2];
temp_mid1=dist.*neigh-[pmidx2 pmidy2 pmidz2];
   
   distances=ones(length(neigh),1)*1000;
    
    for i=1:length(temp_mid)
    
    distances(i)=pdist([50 50 50;temp_mid(i,:)]);
    
    end
    [val1 ind1]=min(distances);
    d1x=temp_mid(ind,1);
    d1y=temp_mid(ind,2);
    d1z=temp_mid(ind,3);
    
    
   distances=ones(length(neigh),1)*1000;
    
    for i=1:length(temp_mid1)
    
    distances(i)=pdist([50 50 50;temp_mid1(i,:)]);
    
    end
    [val2 ind2]=min(distances);
    d2x=temp_mid(ind,1);
    d2y=temp_mid(ind,2);
    d2z=temp_mid(ind,3);


if val2<val1
Xa(2)=d2x;
Ya(2)=d2y;
Za(2)=d2z;
else
Xa(2)=d1x;
Ya(2)=d1y;
Za(2)=d1z;
end

elseif complexity==3
    

Xa(1)=fmidx-dist;
Ya(1)=fmidy-dist;
Za(1)=fmidz-dist;

Xa(2)=midx+dist;
Ya(2)=midy+dist;
Za(2)=midz+dist;

Xa(3)=smidx-dist;
Ya(3)=smidy-dist;
Za(3)=smidz-dist;
    
end


if exist('Xa','var') && exist('Ya', 'var') && exist('Za', 'var') 
points_x=[startp(1) Xa endp(1)];
points_y=[startp(2) Ya endp(2)];
points_z=[startp(3) Za endp(3)];
else
points_x=[startp(1)  endp(1)];
points_y=[startp(2)  endp(2)];
points_z=[startp(3)  endp(3)];
end


[x,y,z] = create_curved_fiber(X,Y,Z,points_x,points_y,points_z);

conn(R1,R2)=conn(R1,R2)+1;
conn(R2,R1)=conn(R2,R1)+1;