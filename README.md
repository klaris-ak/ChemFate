# Welcome to ChemFate!
* ChemFate has been published in *Chemosphere*, please go [here](https://dx.doi.org/10.1016/j.chemosphere.2020.126897) for detailed model documentation.
## ChemFate (Chemical Fate and Transport models)
A level IV dynamic chemical fate and transport tool that can handle four classes of chemicals:
  1) **organoFate** - non-ionizable organic chemicals
  2) **ionOFate** - ionizable organic chemicals
  3) **metalFate** - metals
  4) **nanoFate** - nanomaterials

## 1. Run Model from Python Code
If you would like to run the model using Python, you can download the whole folder on this page by clicking on 'Clone or download' button.
* **run_ChemFate.py** is the python script to call to run the code.

## 2. Run Model from Executable File (.exe)
To run ChemFate through a simple GUI (ChemFate.exe), 
please go [here](https://drive.google.com/file/d/12zNlE2hnfWgw7UkB04ENNXbtYqvNP7t5/view?usp=sharing) to download the tool.
* Please make sure to download the whole zipped folder.

### Below are the steps to run ChemFate tool:
1. Unzip the folder. Please keep the folder structure and keep the files in their specified folder.
	* **ChemFate.exe**: ChemFate executable file to run the tool.
	* **Input** folder: contain the input files.
	* **Output** folder: store the output files.
	* **Model_Data** folder: the data files that are needed during the tool execution.

2. Double click on '**ChemFate.exe**'. It may take a few seconds or longer to open the ChemFate GUI window below.

![ChemFate.exe](https://github.com/klaris-ak/ChemFate/blob/master/Images/chemfate.png "ChemFate.exe")

### Here are a few tips:
1. Please keep the input files in the 'Input' folder.
	- ChemParam_nonionizableOrganic.xlsx
	- ChemParam_ionizableOrganic.xlsx
	- ChemParam_metal.xlsx
	- ChemParam_nanomaterial.xlsx
	- ChemRelease.xlsx
	- Region.xlsx

2. Please don't change the names of files in the 'Input' Folder.
	- if you want to change to another chemical or another region,
	you can just replace the content inside those files.

3. Your model results will be in the 'Output' Folder.

4. If you are going to run ChemFate for 10 yrs of data (2005 1 1 to 2014 12 31),
	it may take around 3-5 mins.

If you have any questions, please contact Dr. Arturo A. Keller (arturokeller@ucsb.edu), Keller Lab, Bren School, UC Santa Barbara, USA
