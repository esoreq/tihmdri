![TihmDRI_tools logo](img/TihmDRI_tools-logo-left-01.png)


`TihmDRI_tools` are currently the easiest way to gain accsess to the data collected by the TIHM projects.

The tools are an ongoing effort by  Dr Eyal Soreq to standardize data-driven analysis for different types of data domains (IoT, behavioural, physiological, etc.) being collected as part of the ongoing TIHM study to advance our understanding of the different manifestations of dementia.
The package was designed to either work on the Imperial cluster or on your personal computer using Python >3.6. 

 
**Installation Options**
---
1. Install with [`pip`](https://pypi.org/project/stronghold/)
    + `$ pip install TihmDRI_tools`

**Usage**
---
See example Jupyter notebooks 


**Inputs**
---

**TIHM data sets**
---

3 zip files
 
Filename | Description | Duration
--- | --- | ---
tihm10.zip | Data collected by TIHM IoT system during TIHM project | (2016-2018)
tihm15.zip | Data collected by TIHM IoT system during TIHM 1.5 project (extension) | (2018-2019)
tihmdri.zip | Data collected by TIHM IoT system during DRI project | (2019-present) - exported once a week and copied to RDS

The TIHM system has undergone numerous iterations during this time. These exported datasets are extracted from historic backups taken at the end of each project. The original databases are MongoDB and do not have a consistent schema. The CSV exports harmonise these into a simplified, consistent tabular format. 


**How to Contribute**
---

1. Clone repo and create a new branch: `$ git checkout https://github.com/esoreq/TihmDRI_tools -b name_for_new_branch`.
2. Make changes and test
3. Submit Pull Request with comprehensive description of changes
