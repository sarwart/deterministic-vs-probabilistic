function [img,vx,vy,vz]=create_vectors(img,X,Y,Z,x,y,z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimates local tangent vector for each point parameterizing the fiber trajectory
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r ,c, l]=size(img);

vx=zeros(r,c,l);
vy=zeros(r,c,l);
vz=zeros(r,c,l);


ind1=floor(x);
ind2=floor(y);
ind3=floor(z);
try
ind=sub2ind([r,c,l],ind1,ind2,ind3);
catch
end
ind=unique(ind);

for i=1:length(ind)
    
    [ii,jj,kk]=ind2sub([r,c,l],ind(i));
    
        if ii+1<=r && jj+1<=c && kk+1<=l
     
         
        check=find((x>=X(ii,1,1)).*(x<=X(ii+1,1,1)).*(y>=Y(1,jj,1)).*(y<=Y(1,jj+1,1)).*(z>=Z(1,1,kk)).*(z<=Z(1,1,kk+1)));
       
         if (length(check)>1)
         
           p0=[x(check(1));y(check(1));z(check(1))]; 
           p1=[x(check(end));y(check(end));z(check(end))];
        
        
           v=p0-p1;
           v=v/norm(v); 
       
           vx(ii,jj,kk)=v(1);
           vy(ii,jj,kk)=v(2);
           vz(ii,jj,kk)=v(3);

          img(ii,jj,kk)=1;
   
         end
      
         end

end




