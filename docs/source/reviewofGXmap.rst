
Review of Gas Exchange Mapping
====================================
A quick review of how gas exchange map is produced.

The process is detailed in
`Quantitative analysis of hyperpolarized 129Xe gas transfer MRI <https://aapm.onlinelibrary.wiley.com/doi/full/10.1002/mp.12264/>`_

Briefly, there are 4 steps. First, the raw FIDs are read out from the Siemens .twix file, which is the raw imaging file. Second, image reconstruction was conducted for non-Cartesian 3D radial images. The data was gridded and interpolated back to Cartesian coordinate and then a FFT was used. Third, the reconstructed gas and dissolved phased images were used to make separate gas, barrier and RBC images. With them created, a binning step was used with thresholds derived from healthy populations. The last step was producing clinical report. The report was first rendered into a html file. The html file was then turned into a pdf and lastly a pptx file. Currently the pptx file is sent out to our clients.

.. _fig2:

    .. image::  _static/Slide3.PNG
       :width: 100%

    Figure 3. A diagram of Gas exchange mapping that is conducted by Zmap main program.
