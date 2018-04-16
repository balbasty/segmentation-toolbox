# Segmentation and template construction of neuroimaging data

This code performs segmentation and template construction of neuroimaging data (CT and MRI). It is an extension of the segmentation framework in the SPM software.

## Folder structure
This code requires a well-defined folder structure of the input data. What follows is an example of this structure. 

Let's say we have eight subjects who each have T1-, T2- and PD-weighted MR data; and an image with categorical labels of tissue types. If the root directory is called *data1* then the folder structure should be:

|-- *data1*  
|-- |-- *S1*  
|-- |-- |-- *scans*  
|-- |-- |-- | -- *T1*  
|-- |-- |-- | -- *T2*  
|-- |-- |-- | -- *PD*  
|-- |-- |-- *labels*    
.  
.  
.  
|-- |-- *S8*  
|-- |-- |-- *scans*  
|-- |-- |-- | -- *T1*  
|-- |-- |-- | -- *T2*  
|-- |-- |-- | -- *PD*  
|-- |-- |-- *labels*   

The folder that contains the imaging data for each subjects must be named *scans* and the folder that contains the labels *labels*. The categorical labels are assumed to be stored in one image.

## Basic use cases

A basic use case is in the main function *build_template.m*. It uses data available in the shared folder *Ashburner_group*.

## Dependencies

This project has strong dependencies to SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [`auxiliary-functions` toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).

## TODO

* Finish label implementation (rem ~= lab)
* Missing data
* Improve registration (w. Yael)
* Make results same on cluster and locally (when so, test summing to gr and H from push_resp)
