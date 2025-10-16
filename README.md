# The cost-effectiveness of an AI stroke imaging tool: The Model

> This is the R project repository for the paper: The health and economic impact of an Artificial Intelligence Software for Stroke Imaging and treatment decisions along the patient pathway in England: A cost-effectiveness analysis

## Description

This repo is a cost-effectiveness model, including all the parameter inputs, R code and outputs, and data files used to estimate the cost-effectiveness of an Artificial Intelligence Software for Stroke Imaging and treatment decisions along the patient pathway in England.

Note the accompanying paper is going through review processes and therefore this repository is subject to change until publication.

## How to Use
Please see the "User_Guide_V1.doc" in this repository.

## Contents

Folder | Description
-----|------------
[inputs](inputs)| parameter input data
[outputs](outputs)| output data from running the model
[model_code](model_code)| the files to run the model. Start --> finish running order should be in prefix numerical order. These run the models at the national level. Note you have to run this order (1_ --> 4_) for the parameter csv table to be formatted properly. For local or regional level run "local_level.R" and "regional_level.R" respectively. Model functions are in the model_functions.R script, with an additional variation used in scenario analyses. 
[figure_creation](model_code/figure_creation)| A sub-folder within model_code which creates figures


> The main project folder also has the following:

-  **Report_Dynamic_Choices.Rmd**: this file is used to create reports based on user selections.

-  **Report_Dynamic_Choices.pdf**: this file is updated when "Report_Dynamic_Choices.Rmd" is run.

- **LICENSE** : this outlines the license applicable for this code (MIT)

- **.gitignore**: this outlines which file types to ignore in GitHub commits

- **README.md**: this file 

Note if you run the "Report_Dynamic_Choices.Rmd" file, a file called "Report_Dynamic_Choices.log" is also created,
this file is generated during the report compilation process due to LaTex and can be ignored.

## How to Cite this Work

### This is work in progress and has not been finalised/signed-off - please do not use yet
Manuscript citation: TBC


## 👂 Further information & Feedback

Please feel free to raise an issue on this GitHub, using the Issues functionality on the GitHub Repository, though the funding this project has been finished so responses and capacity for changes may be limited. 
