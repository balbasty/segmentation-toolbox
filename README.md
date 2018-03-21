# Segmentation and template construction of neuroimaging data

This code performs segmentation and template construction of neuroimaging data (CT and MRI). It is an extension of the segmentation framework in the SPM software.

## Dependencies

This project has strong dependencies to SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [`auxiliary-functions` toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).

## TODO

* Try k-means on histogram accumulated from all CT images
* Finish label implementation (rem ~= lab)
* Missing data
* Improve registration (w. Yael)
* init_obj should directly read pars.segment
* Introduce use of more functions from auxiliary-functions 
* Make results same on cluster and locally (when so, test summing to gr and H from push_resp)
