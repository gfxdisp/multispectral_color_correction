# Color-Accurate Camera Capture with Multispectral Illumination and Multiple Exposures
[Web page](https://www.cl.cam.ac.uk/research/rainbow/projects/color_accurate/) | [paper](https://www.cl.cam.ac.uk/~rkm38/pdfs/gao2024color_acurate_multispectral_capture.pdf)

<img src="teaser.png" width="100%">

We explore an approach in which we can capture accurate colors with a regular camera by optimizing the spectral composition of the illuminant and capturing one or more exposures. We jointly optimize for the signal-to-noise ratio and for the color accuracy irrespective of the spectral composition of the scene. One or more images captured under controlled multispectral illuminants are then converted into a color-accurate image as seen under the standard illuminant of D65.


## Data
1. `SFU_1993.mat` is the SFU surface reflectance data from the paper [A Data Set for Colour Research](https://www2.cs.sfu.ca/~colour/data/colour_constancy_synthetic_test_data/index.html). By default the dataset consists of 1995 spectra from 380 to 780 nm with a 4 nm sampling interval. We removed 2 repetitive spectra and divide all spectra into different categories, including "additional", "dupont", "krinov", "macbeth", "munsell", "objects" and "all".

2. `illum_D65.mat` is the normalized SPD of D65 illuminant extracted from 390 to 780 nm.

3. `pmcc.mat` is our measured spectral reflectance of the [Preferred Memory Colour Chart](https://thouslite.com/product_detail/1827.html) from 380 to 780 nm.

4. `xyz1931.mat` is the 1931 color matching functions from 360 to 830 nm.

5. `A7R3_390_780_a0.001_r40_iter1.mat` and `IDS_390_780_a0.001_r40_iter5.mat` under `SSF_240506` directory are our estimated camera spectral sensitivity functions for our Sony A7R3 and IDS U3-3800CP-C-HQ.


## Calculate optimized illuminants
`run_optimization.m` is used for calculating the optimized illuminants for our two cameras, i.e. Sony A7R3 and IDS U3-3800CP-C-HQ, for which we measured their spectral sensitivity functions (SSFs) in our lab. If you have measured the SSFs of your own camera, you can also get the optimized illuminants for it. The required parameter is the number of exposures `k` where we normally set to 1, 2 or 3.  


## Test real images captured under our optimized illuminants
`run_testing.m` is used for testing the images captured under our optimized illuminants. We captured the RAW camera values of the color checker under our D65 and optimized illuminants under 1, 2 and 3 exposures. We can obtain the quantitative results based on the provided data.