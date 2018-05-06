System Structure
====================================
A quick introduction on how the system is structured, where the data is

Where is the data and what does it include
--------------------------------------------------
The data is in either **Patients** or **Patients_old** folder. You can access them by:

::

    cd ~/Patients
    cd ~/Patients_old

The subject data is all grouped in a folder. **The matlab file (SUBJECT_ID.mat) in each folder has everything you need**.

How is the program structured
--------------------------------------------------
The figure explains how the programs are structured.

.. _fig4:

    .. image::  _static/Slide5.PNG
       :width: 100%

    Figure 4. The program is constructed in 3 layers. The bottom utility layer contains all independent functions that provide fundamental functionalities. The middle layer uses the utility functions to construct class for sophisticated task, such as our clinical subject class and rat class. Each class encapsulated their own method functions. The top layer has the process script for your to run different tasks with selected subjects.
