
Getting Started
====================================
A quick introduction on how to use the fundamental functionalities on Ziyi-PC.

Boot Ziyi-PC and Connect to clinical sharefolder
--------------------------------------------------
You are usually lucky enough to skip this step as the PC is always on. But just in case it is not.

Press the power button on the PC machine. If leave it there, the computer shall boot into Ubuntu system by default (There is also a windows co-existing in the PC, which you may not want to use). When the screen lights up, select user **Ziyi**, and type in password which I believe you know.

The system kicks off on the desktop with a menu bar on the left. To connect to the share folder, open up a terminal. You can find the terminal as a black rectange icon from the menu bar on the left. Type in the following commands:

::

    cd ~/Gas_Exchange
    ./MountDrive.sh

The **MountDrive.sh** will connect the sharefolder through ssh. It will typically ask you for 3 password. The first is the root user password which is the same as the lock screen password. The second and third will be the password to connect to the sharefolder (the .sh script uses my account, you can change it to yours so you may connect with your credential).

Process one or multiple subjects
----------------------------------
You can manually run one or multiple subjects with script **GX_Process_batch.py**. The following commands open up the script.

::

    cd ~/Gas_Exchange
    atom GX_Process_batch.py

In atom, you can modify the script **GX_Process_batch.py** to run one or multiple subjects. Specify the data directory and a list of subject IDs (need to exist in the directory) that you want to process.

::

    data_dir = '/home/ziyiw/Patients/'
    Subject_IDs = ['002104','002105']

When finished, save the file. Go back to the terminal and type the following command to run.

::

    python GX_Process_batch.py

Change mapping parameters such as reference thresholds
-------------------------------------------------------
You can manually most of the mapping parameters (including the reference thresholds) in the script **GX_defineColormaps.py**. The following commands open up the script.

::

    cd ~/Gas_Exchange
    atom GX_defineColormaps.py

**The reference thresholds**, **the montage and histogram plotting parameters** and **the healthy reference values** shown in the report table are recorded in this script. They are imported to the main program in runtime. Specifically, you can change the reference thresholds with the following lines:

.. code-block:: python

    mean_vent = 0.5098
    std_vent = 0.1899
    thre_vent = [mean_vent-2*std_vent, mean_vent-std_vent, mean_vent, mean_vent+std_vent, mean_vent+2*std_vent]

    mean_bar = 0.4859
    std_bar = 0.1484
    thre_bar = [mean_bar-2*std_bar, mean_bar-std_bar, mean_bar, mean_bar+std_bar, mean_bar+2*std_bar, mean_bar+3*std_bar, mean_bar+4*std_bar]

    mean_rbc = 0.2592
    std_rbc = 0.1038
    thre_rbc = [mean_rbc-2*std_rbc, mean_rbc-std_rbc, mean_rbc, mean_rbc+std_rbc, mean_rbc+2*std_rbc]

The **thre_vent, thre_bar, thre_rbc** are the thresholds that are finally used.
