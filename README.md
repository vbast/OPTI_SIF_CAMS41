# OPTI_SIF_CAMS41

The tool OPTI_SIF_CAMS41 is developed to optimize the relationship between the satellite-derived estimates of sun-induced chlorophyll fluorescence (SIF) and the model evaluations of carbon gross primary production (GPP). The SIF data is assumed to be linearly proportional to GPP:

<img src="https://latex.codecogs.com/svg.latex?SIF=a\cdot{GPP}+b"/>.

Parameters <img src="https://latex.codecogs.com/svg.latex?a"/> and <img src="https://latex.codecogs.com/svg.latex?b"/> are considered to be unqie for each plant functional type (PFT) of the model, they are optimized individually per PFT using the pixels with the dominant PFT content. Different options are available for the optimization:

* two optimization algorithms (iterative gradient-based L-BFGS-B method or one-step inverse matrix calculation);
* data smoothing with the running mean window;
* data split into latitude bands (90N-30N, 30N-10S, 10S-90S);
* independent parameter optimization for different seasons (JFD, MAM, JJA, SON).

For more information see the wiki-page : [Wiki](https://github.com/vbast/OPTI_SIF_CAMS41/wiki/)

Authors: Vladislav Bastrikov, Philippe Peylin<br>Laboratoire des Sciences du Climat et de l’Environnement, Gif-sur-Yvette, France
