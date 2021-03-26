
<!-- [![DOI](https://zenodo.org/badge/..../.svg)](https://doi.org/###/zenodo.###) -->

# RAEBURN_Nuclear-versus-cytoplasmic-proteostasis

This repository contains the analysis code associated with the **Proteostasis in Cellular Compartments** project, led by Candice Raeburn. This manuscript has been submitted for publication under the title *The cytoplasmic milieu is more able to mitigate protein homeostasis imbalance than the nucleus: a prime role of the heat shock 40 and 70 chaperone system*.

This manuscript has been submitted as a preprint via BioRxiv [here](biorxiv/link). A link to the final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python 3.7. For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v **8.0**](https://www.graphpad.com/scientific-software/prism/).

## Raw data

For convenience, example datasets are provided here under the ```data``` folder. These data may be used to explore the workflows presented in the ```src``` folder as described below.

## Workflow

Initial preprocessing of the raw microscopy files to extract TIFF images from proprietary formats was completed in [Fiji][1]<sup>[1]</sup> or [ImageJ][2] <sup>[2]</sup> using the [BioFormats importer][3] <sup>[3]</sup> (available via the drag-and-drop interface). Stacked or individual TIFF files were then exported without modification for further processing where necessary, examples of which are provided within the ```data``` folder and referred to from hereon as raw data.

Individual analyses are presented within the ```src``` folder. Where processing order is important for individual analyses, scripts have been numbered and should be run in order before any unnumbered counterparts.

In general, scripts performing similar functions have been labelled consistently across different analyses. Briefly, the main function of each script type is as follows:

| Script name         | Description                                                                                                                                  |
|---------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| initial_cleanup     | Clean filenames of raw data exported from ImageJ to remove extraneous information and provide machine-readable labels                        |
| cellpose            | Define regions of interest for cells, nuclei (and where appropriate aggregates) using the CellPose package                                   |
| define_masks        | Provide defined ROIs as masks overlayed on the original image, with optional user engagement to manually fine-tune segmentation/ROIs         |
| pixel_collection    | Collect intensity value in each channel for individual pixels associated with each ROI type                                                  |
| summary_calculation | Generate calculations for measure of interest e.g. background corrected FRET values or relative enrichment in the nucleus versus the cytosol |

## References

[1]: https://imagej.net/ImageJ2

1. Schindelin, Johannes, Ignacio Arganda-Carreras, Erwin Frise, Verena Kaynig, Mark Longair, Tobias Pietzsch, Stephan Preibisch, et al. “Fiji: An Open-Source Platform for Biological-Image Analysis.” Nature Methods 9, no. 7 (July 2012): 676–82. https://doi.org/10.1038/nmeth.2019.

[2]: https://imagej.net/Fiji

2. Rueden, Curtis T., Johannes Schindelin, Mark C. Hiner, Barry E. DeZonia, Alison E. Walter, Ellen T. Arena, and Kevin W. Eliceiri. “ImageJ2: ImageJ for the next Generation of Scientific Image Data.” BMC Bioinformatics 18, no. 1 (November 29, 2017): 529. https://doi.org/10.1186/s12859-017-1934-z.

[3]: https://docs.openmicroscopy.org/bio-formats/5.8.2/users/imagej/installing.html

3. Linkert, Melissa, Curtis T. Rueden, Chris Allan, Jean-Marie Burel, Will Moore, Andrew Patterson, Brian Loranger, et al. “Metadata Matters: Access to Image Data in the Real World.” Journal of Cell Biology 189, no. 5 (May 31, 2010): 777–82. https://doi.org/10.1083/jcb.201004104.
