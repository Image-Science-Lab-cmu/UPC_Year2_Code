This codebase contains two parts --- (1) simulations for under-panel cameras
(UPC) with microlens arrays; (2) simulations for flare removal in UPC
using high-dynamic range (HDR) imaging.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Required MATLAB toolbox:
Parallel Computing Toolbox

Note: 
BM3D are compiled on macOS system. You might need to recompile BM3D
packages based on your system requirements.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Part I: UPC + MLA
This part of code is inside 'simulate_UPC_MLA' folder.
    compute_SV_PSFs.m:       This script computes the spatially-varying PSFs
                             under UPC with MLA setup. System parameters
                             are saved in 'UPC_MLA_functions/MLA_Params.mat'.
    quantitatie_evaluation.m:This script simulates images captured under
                             UPC+MLA camera, and deblur using Wiener deconv.
                             The evaluation is conduct on 30 images in the
                             'test_data' folder and various SNRs.

Example:
    cd simulate_UPC_MLA;
    compute_SV_PSFs([1]);
    quantitative_evaluation([1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Part II: Flare removal using HDR imaging
This part of code is inside 'simulate_UPC_HDR' folder.
    mainCaptureConfig.m:    This script specifies the pixel pattern of the
                            under-panel camera, and HDR scene used in
                            the simulations.
    mainCapture.m:          This script simulates an exposure stack 
                            captured under the specified UPC. And it 
                            composites an HDR image using this exp stack.
    mainDeblur.m:           This script recovers sharp image from HDR capture
                            by iteratively solving least-square with 
                            estimated spatially-varying noise.

Example:
    cd simulate_UPC_HDR;
    mainCapture;
    mainDeblur;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%