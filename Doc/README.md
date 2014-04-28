Installing R for use with ArcGIS Toolboxes and Scripts
======================================================


1. Make sure that you have installed ArcGIS 10.0 or later and Python
    1. Note that Python should be installed by default when you install
ArcGIS.

    2.  Check for the location of several important directories, which
are installed to `C:\\Program Files` or `C:\\Program Files (x86)` depending
on your operating system.  The directories you should be able to locate
include:
        - `...\ArcGIS\Desktop10.x\bin`
        - `...\ArcGIS\Desktop10.x\arcpy`
        - `...\ArcGIS\Desktop\ArcToolbox\Scripts`
        - `..\Python26\ArcGIS10.0` or `..\Python27\ArcGIS10.1` or `..\Python\ArcGIS10.2`

    3. When you install ArcGIS and Python the above directories should be added automatically to your Python 
       Path (an environment variable, which you will learn how to set in the next set of steps)
        - To make sure that these directories were added to your Path environment variable automatically, 
        open ArcMap, open the Python Window (from the Geoprocessing Dropdown), and type the following:

        ```
                >>> import sys
                >>> import sys.path
        ```

        - The output should include the following (the exact output will be different, but should include the pathnames listed above)
    ![arcgis-python-path](./InstallDoc_files/image001.png)
    - If you do not see the directories listed (which is unlikely), add them in Step XX below
    - Close ArcMap 

2. Install R

    1. Go to the R homepage ([http://www.r-project.org](http://www.r-project.org)) and download the most recent version 
    by clicking the "Download R" link in the Getting Started section, and then choosing the location that is closest to you
        - Choose to Download and Install R (NOT the Source Code for all Platforms), choose the appropriate operating system, choose to download the Base install, and then choose to Download R (this will be the most recent version of R)
    2. Follow the standard steps for installation
        - IMPORTANT: Take note of the directory where you choose to install R and write it here (or on a separate piece of paper so that you don't forget... you will need it to add the directory to the Path environment variable) <br> R installation path: <input type="text" size="60" /> 

        - Accept all of the defaults, and run the installation 

3. Install the necessary R Libraries (maptools, clustTool, SM, and Design/rms)

    1. For R versions < 2.14, you will need the Design module.  For R &gt;= 2.14 this package was renamed rms.  The R script should be robust to which ever you have on your system.
    2. Open R, either from the desktop shortcut that was created (if you chose to create one), or from the Start menu
    3. Go to the Packages menu, then choose to Install Packages, then choose the location closest to you for download
    4. Find clustTool, maptools, SM, and Design in the list of packages and select them both (you can select multiple packages using the Ctrl key), and click OK
    5. Close R, without saving the workspace image
    6. You are now ready to run the ArcGIS Sample R Scripts in ArcMap (which can be found in the Geoprocessing Resource Center in the Model and Script tool gallery here)

4. Add the R Tools toolbox and start using it

    1. Open ArcMap
    2. Open ArcToolbox: from the Geoprocessing Window select ArcToolbox
    3. Right click inside of ArcToolbox, and choose Add Toolbox
    4. Navigate to the location where you saved the Using R in ArcGIS 10 zip file and unzipped everything, and choose R Tools.tbx
    5. Double-click on the Point Clustering tool, and fill out of the necessary parameters
        - This tool can only take Shapefiles as inputs, and can only output Shapefiles
    6. Run the tool
        - Ignore the R windows that may pop open, that is just part of the way that the tool functions

Modifying the Path Environment Variable
---------------------------------------

Current versions of this script automatically detect and set the correct path for your R installation. However, if the scripts are unable to find your installation of R, these steps allow you to directly set the location of your R installation. You can follow these same steps to add ArcGIS directories to your path if they did not appear in step 1.3.

1. Navigate to the location of the R install on your computer (that you wrote on the line above) and open the folder labeled Bin
2. Copy the full path location of the bin folder (by default it would install here: `C:\Program Files (x86)\R\R-2.11.1\bin`) 
3.  Go to your Start Menu, find the My Computer icon (or just Computer icon, depending on what operating system you are using), right click on it, and choose the Properties option

    ![properties-start-menu](./InstallDoc_files/image003.png)

4. For both Windows Vista and Windows 7, this will open a dialog where you will then click the Advanced System Settings link on the left side (highlighted in red).  For Windows XP, this will open the System Properties dialog box, from which you will choose the Advanced tab (highlighted in red in)
    ![properties-menus](./InstallDoc_files/image005.png)

5. The next window will have a button on the bottom that says Environment Variables, click it and the following dialog will appear

    ![environment-variables](./InstallDoc_files/image007.png)

6. Find the environment variable called Path (or PATH).  This may be in the User Variables, and it may be in the System Variables.  If you are the only person that uses your machine, this distinction should not matter.  If your machine is shared, and you use the Path variable that is a user variable, then you will be the only person impacted by the change.  If you want your changes to be shared by all users, make sure to use the System Variable.  

7. Select the Path Variable, and click Edit

    ![environment-variables](./InstallDoc_files/image009.png)

7. At the end of the list of directories, add a semi-colon (;) with NO space after it, and then paste in the R Bin directory that you copied earlier (should look something like this: `C:\Program Files (x86)\R\R-2.11.1\bin`) 

8. Click OK

9.  Your R Bin directory should be successfully added to your Path environment variable
