#!/bin/bash

# Mount into the share drive

sudo mount -t cifs -o domain=DHE,username=zw73,vers=1.0 //duhsnas-pri.dhe.duke.edu/xenon_MRI_raw /media/rawdata

sudo mount -t cifs -o domain=DHE,username=zw73,vers=1.0 //duhsnas-pri.dhe.duke.edu/duhs_radiology/Private/TeamXenon /media/sharedrive
