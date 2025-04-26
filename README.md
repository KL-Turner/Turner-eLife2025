# Type-I nNOS neurons orchestrate cortical neural activity and vasomotion

This document outlines the steps necessary to generate the figures for the manuscript **Type-I nNOS neurons orchestrate cortical neural activity and vasomotion** by Kevin L Turner, Dakota F Brockway, Md Shakhawat Hossain, Keith R Griffithm Denver I Greenawalt, Qingguang Zhang, Kyle W Gheres, Nicole A Crowley, and Patrick J Drew.

---
## Generating the figures
---
The code in this repository generates all main figures, supplemental information, and statistics. It also contains the analysis code for all pre-processing steps.

Begin by downloading the code repository by clicking the green **Code** button at the top of the page and selecting **Download ZIP**. 
* Code repository location: https://github.com/KL-Turner/Turner-eLife2025

![](https://private-user-images.githubusercontent.com/30758521/437243955-d1bc8eb4-60c5-4db9-a08e-3e473504e80e.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDU1NDAyNDksIm5iZiI6MTc0NTUzOTk0OSwicGF0aCI6Ii8zMDc1ODUyMS80MzcyNDM5NTUtZDFiYzhlYjQtNjBjNS00ZGI5LWEwOGUtM2U0NzM1MDRlODBlLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA0MjUlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNDI1VDAwMTIyOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWQ4OGEzODFiNWMyMzJhMTE5MmY2Y2IxYjY2YmZjYmZhMzYyMWU3YmU5ZTYzY2U4NDQ2ZjVjMzdhODIxNTllOGYmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.XHN6xsouhaH43EB6i-fwLzzWghQTZ5J9tl4E6PDh6Mw)

Next, the data (~8 GB) can be downloaded from the following location: https://doi.org/10.5061/dryad.tb2rbp0bq

![](https://private-user-images.githubusercontent.com/30758521/437244656-6544bbcf-6631-42a0-9380-fc69af9e0e74.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDU1NDA1MzgsIm5iZiI6MTc0NTU0MDIzOCwicGF0aCI6Ii8zMDc1ODUyMS80MzcyNDQ2NTYtNjU0NGJiY2YtNjYzMS00MmEwLTkzODAtZmM2OWFmOWUwZTc0LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA0MjUlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNDI1VDAwMTcxOFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTE1YjkzMTA5MGY4MTQxMmY2MDdlZDBjNWRjNzIyMjM4ZjVmZTc5ZThjZmQyMTY4MTU4ODBiMzMzN2RiNzliNDAmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.n78OWnmj9iSagTyX9h07qZYpONHMHyTW2vM7M4K5gr8)

After downloading both the code and data, unzip both folders by right clicking and selecting *Extract all* or double clicking each folder. The unzipped folders should look something like this.

![](https://private-user-images.githubusercontent.com/30758521/437248046-0e733d7a-10d7-4372-b861-b4ddca99c5ad.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDU1NDE3NDgsIm5iZiI6MTc0NTU0MTQ0OCwicGF0aCI6Ii8zMDc1ODUyMS80MzcyNDgwNDYtMGU3MzNkN2EtMTBkNy00MzcyLWI4NjEtYjRkZGNhOTljNWFkLnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA0MjUlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNDI1VDAwMzcyOFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWVlZGQ0M2I5OGJjZjY1YTIxYjFlNzVhMDM5MDBkM2Y1NTI1OTk2Mjg5ZjNhY2IxOWQ1YTkzZTc0ZGU5YjA2MTYmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.9ZYoweCJ6KdfMv3eUteT4XpR6gtNRnUvZmHZ44R57yE)

The Dryad link contains several pre-analyzed structures that can be used to immediately generate the figures without re-analyzing any data. To generate the figures, begin by moving all of the files from the unzipped Dryad folder (the **Results_** folders) into the code repository. It should look something like this:

![](https://private-user-images.githubusercontent.com/30758521/437248578-c5cdcec0-7358-425c-ba28-79595cdac095.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDU1NDE5NjAsIm5iZiI6MTc0NTU0MTY2MCwicGF0aCI6Ii8zMDc1ODUyMS80MzcyNDg1NzgtYzVjZGNlYzAtNzM1OC00MjVjLWJhMjgtNzk1OTVjZGFjMDk1LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA0MjUlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNDI1VDAwNDEwMFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPWMwMzViODhjMGRjZGRlZGRiOWIxY2M5YmUzOTk4ZGIxMWYzMDY0ZDI5NDkzNWZiNGJhZmFkYzE5MjI3Y2ZmY2MmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.opxtFi7wCxzoyYcbdtdQZZEFzyzjNvXAWsheH0Q6hak)

## Results_Turner (.mat - MATLAB2024b; .xlsx - MS Excel)
The repository includes pre-analyzed structures for each analysis performed in the paper, each of these has the prefix "Results_" followed by the analysis performed.
* Results_ArousalStateProb_Ephys.mat
* Results_ArousalStateProb_GCaMP.mat
* Results_Baseline_2P.mat
* Results_BaselineShift_2P.mat
* Results_BilatCoher_Ephys.mat
* Results_BilatCoher_GCaMP.mat
* Results_CrossCorr_EGFP.mat
* Results_CrossCorr_Ephys.mat
* Results_CrossCorr_GCaMP.mat
* Results_Derivative_Ephys.mat
* Results_Diameter_2P.mat
* Results_Evoked_2P.mat
* Results_Evoked_Ephys.mat
* Results_Evoked_GCaMP.mat
* Results_Evoked_Pulse.mat
* Results_GFP.mat
* Results_HbTCrossCorr_Ephys.mat
* Results_HRF_Ephys.mat
* Results_IntSig_Ephys.mat
* Results_IntSig_GCaMP.mat
* Results_IntSig_Pulse.mat
* Results_ModelAccuracy_Ephys.mat
* Results_ModelAccuracy_GCaMP.mat
* Results_NeuralHemoCoher_Ephys.mat
* Results_NeuralHemoCoher_GCaMP.mat
* Results_PearsonCorr_Ephys.mat
* Results_PearsonCorr_GCaMP.mat
* Results_PowerSpec_2P.mat
* Results_PowerSpec_Ephys.mat
* Results_PowerSpec_GCaMP.mat
* Results_PowerSpec_LFP.mat
* Results_PreWhitenedPowerSpec_2P.mat
* Results_PreWhitenedPowerSpec_Ephys.mat
* Results_PreWhitenedPowerSpec_GCaMP.mat
* Results_PupilArea_Ephys.mat
* Results_PupilBlinkInterval_Ephys.mat
* Results_PupilCoher_Ephys.mat
* Results_PupilCrossCorr_Ephys.mat
* Results_PupilEvoked_Ephys.mat
* Results_Transitions_Ephys.mat
* Results_Transitions_GCaMP.mat
* Results_WhiskBehav_Ephys.mat
* Results_WhiskBehav_GCaMP.mat

Data sets used for representative examples (Fig. 2)
* T165_210222_12_45_38_ProcData.mat
* T165_210222_12_45_38_SpecDataA.mat
* T165_210223_12_37_03_ProcData.mat
* T165_210223_12_37_03_SpecDataA.mat
* T192_RH_210518_11_41_23_005_A08_MergedData.mat
* T192_RH_210518_11_41_23_005_A08_SpecData.mat
* T233_220206_13_07_08_ProcData.mat
* T233_220206_13_07_08_SpecDataA.mat
* T233_220206_14_10_44_ProcData.mat
* T233_220206_14_10_44_SpecDataA.mat

Immunohistochemistry and Histological characterization data
* DAPI_Counts.xlsx
* DAPI_Fluorescence.xlsx
* Diaphorase_Counts.xlsx
* GFAP_Fluorescence.xlsx
* IBA1_Counts.xlsx
* MasterMouseList.xlsx
* NeuN_Fluorescence.xlsx
* nNOS_Counts.xlsx

## Results_Zhang (.mat - MATLAB 2024b)
14 Folders corresponding to experimental IDs: T192, T200, T205, T206, T208, T209, T211, T212, T213, T216, T217, T218, T219, T225. Each folder contains the analyzed data from three runs (all .mat structures with name, condition, file date) for spectroscopy experiments during locomotion.

## Results Hossain (.mat - MATLAB 2024b; .png - ImageJ/FIJI; .xlsx - MS Excel;)
* 55 data sets (AnimalID_Run1.mat or AnimalID_Run1.png) of open-field, exploratory * behavior.
* ExperimentalRecords_Keys.xlsx - summary of animal information
* OFTTable.mat - Analysis of open field exploratory behavior

From here, open MATLAB and nativate to the code's folder. Open the function **MainScript_eLife2025.m**. Click the green **Run** button and the figures will then take a few minutes to generate.

**Software/System Requirements:** Code was written and tested with MATLAB 2024b. Running **MainScript_eLife2025.m** took ~5 minutes to run on a 2021 Macbook Pro with M1 Pro chipset and 16 Gb RAM.

The Mainscript will automatically save all of the MATLAB figures and statistical read-outs, which will be saved in a folder called MATLAB Figs & Stats. We are happy to provide any raw data and unprocessed histological images upon request.

LabVIEW code used to acquire the data can be found at: https://github.com/DrewLab/LabVIEW-DAQ

---
## Acknowledgements
---
* Chronux Toolbox: http://chronux.org/
* ConvertTDMS.m by Brad Humphreys: https://github.com/humphreysb/ConvertTDMS
* colors.m by John Kitchin: http://matlab.cheme.cmu.edu/cmu-matlab-package.html
* BrainMesh.m by Yaoyao Hao: https://github.com/Yaoyao-Hao/BrainMesh
* multiWaitbar.m by Ben Tordoff: https://www.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin

#### Contact Patrick Drew or Kevin Turner (kevin_turner@brown.edu) with any questions/comments.

