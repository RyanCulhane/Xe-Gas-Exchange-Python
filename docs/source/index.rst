.. Zmap documentation master file, created by
   sphinx-quickstart on Fri May  4 16:56:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Zmap's documentation!
====================================

Zmap is a Python based program developed to facilitate \ :sup:`129`\ Xe MRI scan and analyze images. The goal of Zmap is to build a **fully automated** analysis pipeline from the image data created from MR scanner to a **clinical report** that contains all meaningful information from this subject.

Currently Zmap provides functions in calibration and image process.

1. The calibration requires only the calibration twix file from Siemens, and the result will be a text file created at the folder of the calibration file. A  SMS text will also be sent to the pre-registered clients.

2. The final clinical report, generated after the image processing, will be sent all through email to the pre-registered clients.


Quick Start
=============
A quick introduction on how to use the fundamental functionalities on Ziyi-PC.



* :ref:`getting_started`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Welcome to Zmap:

   gettingstarted
