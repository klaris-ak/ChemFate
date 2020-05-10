# Welcome to ChemFate!
* ChemFate has been published in *Chemosphere*, please go [here](https://dx.doi.org/10.1016/j.chemosphere.2020.126897) for detailed model documentation.
## ChemFate (Chemical Fate and Transport models)
ChemFate predicts daily chemical environmental concentrations for four classes of chemicals in nine environmental compartments (e.g., air, freshwater, coastal seawater, natural soil, urban soil, agricultural soil with and without applied biosolids). ChemFate comprises four different models:

  1) **organoFate** - a model for non-ionizable organic chemicals
  2) **ionOFate** - a model for ionizable organic chemicals
  3) **metalFate** - a model for metal ions
  4) **nanoFate** - a model for nanomaterials

## 1. Run Model from Python Code
If you would like to run the model using Python, you can download the whole folder on this page by clicking on 'Clone or download' button.
* **run_ChemFate.py** is the python script to call to run the code.

## 2. Run Model from Executable File (.exe)
To run ChemFate through a simple GUI (ChemFate.exe), 
please go [here](https://drive.google.com/file/d/12zNlE2hnfWgw7UkB04ENNXbtYqvNP7t5/view?usp=sharing) to download the tool.
* Please make sure to download all four folders/files within the ChemFate folder. You can download all in one zipped file.

### 1) Below are the steps to run ChemFate tool:
1. Unzip the folder. Please keep the folder structure.
	* **ChemFate.exe**: ChemFate executable file to run the tool.
	* **Input** folder: contain the input files.
	* **Output** folder: store the output files.
	* **Model_Data** folder: the data files that are needed during the tool execution.

2. Double click on '**ChemFate.exe**'. It may take a few seconds or longer to open the ChemFate GUI window below.

![ChemFate.exe](https://github.com/klaris-ak/ChemFate/blob/master/Images/chemfate.png "ChemFate.exe")

### 2) Here are a few tips:
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

## 3. Model Output Results Example
ChemFate generates the following outputs:
| Content       | File Tyle     |
| ------------- |:-------------:|
| Daily concentration of each chemical’s various forms (if applicable) in each compartment      | XLSX |
| Daily mass of each chemical’s various forms (if applicable)  in each compartment     | XLSX       |
| Heatmap showing comparison of each chemical’s overall average concentration for all compartments| PNG     |
| Line graphs showing each chemical’s concentration variability over time| PNG     |

### 1) Heatmap example

![heatmap](https://github.com/klaris-ak/ChemFate/blob/master/Images/metal_heatmap_example.png "heatmap")

### 2) Line graph example

![linegraph](https://github.com/klaris-ak/ChemFate/blob/master/Images/neutral_bulk_example.png "linegraph")

*If you have any questions, please contact Dr. Arturo A. Keller (arturokeller@ucsb.edu),Bren School, UC Santa Barbara, USA*
