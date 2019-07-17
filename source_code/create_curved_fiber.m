function [x,y,z] = create_curved_fiber(X,Y,Z,pt1,pt2,pt3)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generating intermediate points for complete fiber trajectory specification
%Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. 
%Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. 
%Magnetic resonance in medicine, 81(2), pp.1368-1384.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t  = 1:numel(pt1);
Xa=zeros(length(t),1);
Ya=zeros(length(t),1);
Za=zeros(length(t),1);
for i=1:numel(pt1)
    Xa(i)=X(pt1(i),1,1);
    Ya(i)=Y(1,pt2(i),1);
    Za(i)=Z(1,1,pt3(i));
end
ts = linspace(min(t),max(t),10000); 

x = (spline(t,Xa,ts))';
y = (spline(t,Ya,ts))';
z = (spline(t,Za,ts))';






