# BDA
****************************************************************************************
**BDA - A novel method for identifying defects in body-centered cubic crystals**
****************************************************************************************

The accurate and fast identification of crystallographic defects plays a key role for the analysis of atomistic simulation output data. For face-centered cubic (fcc) metals, most existing structure analysis tools allow for the direct distinction of common defects, such as stacking faults or certain low-index surfaces. For body-centered cubic (bcc) metals, on the other hand, a robust way to identify such defects is currently not easily available.

We therefore introduce a new method for analyzing atomistic configurations of bcc metals, the BCC Defect Analysis (BDA). It uses existing structure analysis algorithms and combines their results to uniquely distinguish between typical defects in bcc metals.

In essence, the BDA method offers the following features:
* Identification of typical defect structures in bcc metals.
* Reduction of erroneously identified defects by iterative comparison to the defects in the atom's neighborhood.
* Availability as ready-to-use Python script for the widespread visualization tool OVITO (http://ovito.org).

The BDA method was developed and written by:   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Johannes J. Möller (johannes.moeller@fau.de),  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Materials Science and Engineering, Institute I,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU), Germany.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; http://www.gmp.ww.uni-erlangen.de

The Python script can be downloaded from GitHub under the link given above. If you used the BDA method to analyze your simulation results, please cite the BDA in your publications as follows:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; J.J. Möller and E. Bitzek  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; BDA: A novel method for identifying defects in body-centered cubic crystals  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; MethodsX 3 (2016), 279-288  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; http://dx.doi.org/10.1016/j.mex.2016.03.013

If possible, please also include a link to this website:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; http://jomoeller.github.io/bda/

The author greatfully acknowledges the supervision and support of   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Prof. Dr. Erik Bitzek (erik.bitzek@fau.de),  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Materials Science and Engineering, Institute I, FAU Erlangen-Nürnberg, Germany.

****************************************************************************************
