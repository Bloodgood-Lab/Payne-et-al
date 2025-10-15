# Payne et al. — Analysis Code Repository

This repository contains the analysis code used for the analyses presented in *Payne et al.*, unpublished manuscript.

---

## Overview

The repository includes one main MATLAB workflow:

**`inVivoProject_masterWorkflow_v1_20240718.m`**

This file serves as the backbone of the analysis pipeline.  
It calls a variety of custom functions, which can be found in the **`functions`** folder.

---

## Workflow Structure

Within the main workflow, the major analyses from the paper are organized into sections.  
Each section includes:

1. **Inputs** – The data or files used to run that section.  
2. **Settings** – Parameters and configurations used to generate the results.  
   - Alternative settings are included in comments where applicable.  
3. **Outputs** – The processed data and any corresponding plots.  

As part of the output, each section generates a `.mat` file that contains:
- A **data structure** with all processed variables accumulated to that point.  
- A **settings structure** with the parameters used up to that section.

Subsequent sections load the `.mat` file from the previous section, so the workflow should be run **in series**.  
Once a section has been completed, its results are stored in the `.mat` file and do **not** need to be re-run.

---

## Environment

This code was originally written and executed in **MATLAB R2015b**.

---

## Data Availability

Data are available upon request.

---

## Contact

For questions or requests, please contact:  
**Anja Payne** — [anja.e.payne@gmail.com](mailto:anja.e.payne@gmail.com)
