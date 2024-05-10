# The Role of Policy and Module Manufacturing Learning in Industrial Decarbonization by Small Modular Reactors
Max Vanatta<sup>1</sup>, Michael T. Craig<sup>1,2</sup>, William R. Stewart<sup>3</sup>  
<sup>1</sup> University of Michigan, School for Environment and Sustainability, [ASSET Lab](https://www.assetlab.org)  
<sup>2</sup> University of Michigan, Department of Industrial and Operations Engineering  
<sup>3</sup> [Boston Atomics, LLC.](https://www.bostonatomics.com)  

## General Description
This is the code and support files for NEUP research for industrial heat decarbonization pathways using SMRs with learning and policy contexts. This work builds upon previous work by Vanatta et al. [^1]. This Repo is a work in progress. For more inforamtion please contact me, Max, directly.

The purpose of this code is to develop deployment pathways for SMRs at industiral facitlies from NREL's analysis of GHGRP facilities while considering the impact of policy and learning rates (factory and onsite). 
Functionally, this code will 1) accept a set of conditions in which it 2) pulls industrial facility data, 3) performs a rolling horizon, per-facility mixed integer quadradic program (miqp) to find the most profitable SMR deployment, and 4) outputs the resulting system decisions. 

1. The input conditions for the analysis are the context of the model run such as the number of facilties to be tested, the number of total years for each optimization, number of iterations to perform the rolling horizon over, etc... These inputs are the foundation of the scenario analysis performed in the publication of this work.
2. With the selected number of facilities and a specified facitity file (multiple included here) the model will generate the landscape of potential investment sites. The facilites provided are derived from the work by Colin McMillan's team at NREL[^2] which took the EPA's Greenhouse Gas Reporting Program data[^3] and converted emissions to idnutrial process ehat demands by fuel type. These facilties are those which emit more than 25,000 mt CO<sub>2</sub>e per year. To substitute your own facilties, it is important to match the requisite column names. Template provided, but this is not recommended as it would require NAICS codes which correspond to the NREL work to be processed.
3. Each year the model runs every facility, and sub temperature range, and finds the most profitable investment in competition with natural gas in isolation. The deployments are stacked in order of profitability and deployed until the annual manufacturing limit is reached, remaining modules are left unbuilt. The model adjusts the capital costs due to the built modules by SMR design for the following year/iteration. The process repeats yearly until the iteration count (RollCount) is complete.
4. The output of the model can be interpreted by the included jupyter notebooks which take the build decisions of the model and output some useful figures and data files. More data analysis is available upon request.

## Basic Operation
To use the code, locate all included files within the same directory and ensure that you have installed all dependencies:
Pandas, numpy, matplotlib, seaborn, plotly, pyomo, branca, geopandas
  
If the use is for unsubsidized deployment, V10_GL_0403c.py is the better code to run as it is more stable. If subsidies are the aim, V12_GL.py is required.
ALL OTHER VERSIONS ARE TROBLESHOOTING (I apologize for the mess there)
### Inputs
Run the V10, V12 files with the following required command line arguments:  
  **--Rand_n** _int_, Number of random facilities to sample (randomization is consistent between runs, ie, facility 1 will always have the same id)  
  **--YearStart** _int_, Starting Year of analysis  
  **--YearCount** _int_, Number of years to include in the discrete optimization   
  **--HrCount** _int_, Number of hours per year  
  **--NGP** _str_, Natural Gas Price ($/MMBtu)  
  **--RollStep** _int_, Rolling horizon step size in years  
  **--RollCount** _int_, Rolling horizon count  
  **--Max_Builds** _str_, Minimum build requirements formatted as a string  
  **--ITC** _float_, Investment tax credit percentage (0-1)  
  **--FolderName** _str_, Makes the filename and directory  
  **--IR** _float_, Interest Rate  
  **--Lifetime** _int_, SMR book lifetime in years  
  **--ITC_Duration** _int_, Years of ITC operation  
  **--Min_Builds** _str_, Minimum build requirements (just use 0,0,0,0) 

And the optional arguements: 
  **--NGP_Variable** _str_, Is natural gas price changing  
  **--NGP_V** _str_, Natural gas prices over the years  
  **--Max_Increase** _float_, Percent of module maximum growth per year  
  **--LearningRate** _str_, Learning Rate Scenario  
  **--Gen_File** _str_, Name of generator parameter file, default = 'LC_base_4pack.csv'  
  **--Init_N** _str_, Number of SMR modules previously built  
  **--CarbonTax** _float_, Carbon Tax in dollars per ton CO2  
  **--FacilityFile** _str_, File for the facilities', default = '2015_NG_National.csv'  
  (V12 only) **--Subsidy** _float_, Per SMR type subsidy quantity which will be distributed out to each until consumed, default = 0  
  **--CC** _float_, Capital Cost Modifier for sensitivites, default = 1  

Example:  
  `python ./V10_GL_0403c.py --Rand_n 925 --YearStart 2030 --YearCount 40 --HrCount 1 --NGP 10 --RollStep 20 --RollCount 20 --Min_Builds 0,0,0,0 --ITC 0.0 --ITC_Duration 3 --FolderName Git_example_10 --IR 0.07 --Lifetime 30 --Max_Increase 0.20 --Max_Builds 2,8,2,1`

### Outputs
The output of this code will be:  
  (folder) 'FolderName'largeoutputs: This is the per-heat range, per-facility resutls of natural gas supply and SMR output.  
  'FolderName'_ binaries_loggin_All.csv: This is the large collection of post-rolling horizon limitation Build, Construction, and Operation binaries.  
  'FolderName'_ binaries_loggin.csv: This is the reduced set of post-rolling horizon limitation Build binaries for selsected facilities.  
  profit.csv: This file tracks the best objective function for each facility for each rolling horizon iteration.  
  demand_selection.csv: This file tracks the selected temperature band for each facility for each rolling horizon iteration.  

## Acknowledgements
For funding, we thank Idaho National Laboratory’s Emerging Energy Markets Analysis initiative, and the U.S. Department of Energy Office of Nuclear Energy’s Nuclear Energy University Program under contract number DE-NE0008976.

## References
[^1]:Vanatta, M., Patel, D., Allen, T., Cooper, D. & Craig, M. T. Technoeconomic analysis of small modular reactors decarbonizing industrial process heat. Joule 7, 713–737 (2023). https://www.sciencedirect.com/science/article/pii/S2542435123001216?via%3Dihub  
[^2]:McMillan, C. & Ruth, M. Industrial Process Heat Demand Characterization. https://data.nrel.gov/submissions/91 (2018) doi:10.7799/1461488.  
[^3]: US EPA. Inventory of U.S. Greenhouse Gas Emissions and Sinks | US EPA. https://www.epa.gov/ghgemissions/inventory-us-greenhouse-gas-emissions-and-sinks (2022).  
