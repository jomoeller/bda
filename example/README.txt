###################################################################################
# Example for the BCC Defect Analysis                                             #
# Penny-shaped crack on (010) plane in Fe modelled with the Mendelev-II potential #
# Downloaded from http://jomoeller.github.io/bda                                  #
#                                                                                 #
# Jo Moeller, Thu Apr  7 14:11:36 CEST 2016                                       # 
# Please cite the BDA:                                                            # 
#     J.J. Moeller and E. Bitzek                                                  #       
#     BDA: A novel method for identifying defects in body-centered cubic crystals #       
#     MethodsX 3 (2016), 279-288                                                  #       
#     http://dx.doi.org/10.1016/j.mex.2016.03.013                                 #  
###################################################################################

1. Unzip the splitted .gz file:
===============================

$> cat example.chkpt.gz.part-* > example.chkpt.gz
$> gunzip example.chkpt.gz

2. Run the BDA with ovitos:
===========================

$> ovitos /path/to/bda/ovitos_bcc-defect-analysis.py -c example.chkpt -b 0 0 0 -p Men-II

Note: the command "ovitos" might be different on your system, e.g., /opt/ovito-2.6.1-37-x86_64/bin/ovitos

3. Open the analyzed file with ovito:
=====================================

$> ovito example.chkpt.bda &

4. Load the BDA color scheme
============================

In ovito do the following:
* enable the 'color coding' modifier in the modification pipeline on the right;
* change the color scheme by selecting 'Load custom color map' in 'Color gradient:';
* navigate to the image 'colorbar_cubehelix_0-6.png' in the root directory of the BDA repository;
* adjust the range of the color map such that it contains the numbers from 0 to 6. 

Et voila!

Appendix: Commands to create the splitted .gz file:
===================================================

$> gzip nov_cut_nve_defrate1e8_psc_r100_Fe_Men-II_x100y010z001_800x1000x800_e0465.00038.chkpt
$> split -b 25m nov_cut_nve_defrate1e8_psc_r100_Fe_Men-II_x100y010z001_800x1000x800_e0465.00038.chkpt.gz example.chkpt.gz.part-
