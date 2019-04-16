# Application-of-Fractional-Calculus-in-Multispectral-Image-Fusion

[Application of fractional-order differentiation in multispectral image fusion
](https://www.tandfonline.com/doi/abs/10.1080/2150704X.2017.1395963) is the paper for this MATLAB code.


Description
----------
This code provides the fusion of PANchromatic (PAN) and MultiSpectral (MS) images using fractional-order differentation. The steps for fusion is as follows: 

    1) Loading the dataset from its path,
    
    2) Pre-processing steps including downsampling and normalization,
    
    3) Approxiamtion of fractional-order differentiation and creating the superimposed mask,
    
    4) Obtaining the primitive detail map for each spectral band,
        
    5) Applying the created mask into the primitive detail map to obtain the refined map.


Usage
------------
First you need to specify the path of your dataset.
For example:

    addpath QuickBird_Data
    
The main code for this method is the FDIF.m which consists of pre-processing steps as well as obtaining the fusion results. The superimposed mask is created in another M-file called FractionalDiff. 
To _run_ the code, in the _command window_ use this: 

    FDIF.m

Objective Evaluation
----------
For objective assessment of the fusion result, first add the path of objective metric. For example: 

    addpath Objective_Evaluation

Sample Output
-----------
    
    The MS, PAN and pansharpened result of the QuickBird dataset from Sandarbans, Bangladesh. 

![LRMS](https://user-images.githubusercontent.com/48659018/56171542-5284ec00-5fab-11e9-93fb-a973ba1e8014.png)
![PAN](https://user-images.githubusercontent.com/48659018/56171559-603a7180-5fab-11e9-8626-c1103ca22e6d.png)
![Pansharpened](https://user-images.githubusercontent.com/48659018/56179713-5b39ea00-5fcc-11e9-96ca-56f1a3eb16f5.png)


Reference
--------
If you find this code helpful, please cite this paper: 

    A. Azarang and H. Ghassemian, "Application of fractional-order differentiation in multispectral image fusion," Remote Sens. Lett., vol. 9, no. 1, pp. 91-100, Jan. 2018.

Contact
--------
If you have any question regarding the paper, codes and so on, do not hesitate to contact me. 

Arian Azarang  azarang@utdallas.edu
