
# **Generating Spherical Phantom**

## *"Sarwar, T., Ramamohanarao, K. and Zalesky, A., 2019. Mapping connectomes with diffusion mri: deterministic or probabilistic tractography?. Magnetic resonance in medicine, 81(2), pp.1368-1384."*

![alt text](https://github.com/sarwart/Phantoms/blob/master/Image.png)

The study evaluated the state-of-the-art tractography algorithms using simulated diffusion-weighted magnetic resonance imaging (dMRI) dataset. 
The details of simulation can be in paper "https://onlinelibrary.wiley.com/doi/10.1002/mrm.27471"
The provided scripts (main script: "demo.m") can be used to generate spherical phantom with 60 nodes. Node parcellation scheme (atlas.nii.gz) for mapping connetcome is also provided. The dataset provided is simulated with b-value=2000 s/mm2 and provided gradients. Rest of the parameters value can be found in the paper.


A dataset of 5 spherical phantoms can be found in folder "Data". The parcellation and diffusion-weighted gradients along with the masks are also included. 
One test phantom, with only one fiber connetcing a pair of region is also available. This could be used as a reference to deal with flipping (in case) caused by tractography algorithms. 




