# SMR_Learning
This is the code and support files for NEUP research for industrial heat decarbonization pathways using SMRs with learning and policy contexts. This Repo is a work in progress. For more inforamtion please contact me, Max, directly.

The prupose of this code is to develop deployment pathways for SMRs at industiral facitlies from NREL's analysis of GHGRP facilities while considering the impact of policy and learning rates (factory and onsite). 

To use the code, locate all included files within the same directory and ensure that you have installed all dependencies:
Pandas, numpy, matplotlib, seaborn, plotly, pyomo

Run the V10, V12 files with the following required command line arguments:

'--Rand_n', type=int, Number of random facilities to sample (randomization is consistent between runs, ie, facility 1 will always have the same id)
'--YearStart', type=int, Starting Year of analysis
'--YearCount', type=int, Number of years to include in the discrete optimization 
'--HrCount', type=int, Number of hours per year
'--NGP', type = str, Natural Gas Price ($/MMBtu)
'--RollStep', type=int, Rolling horizon step size in years
'--RollCount', type=int, Rolling horizon count
'--Max_Builds', type=str, Minimum build requirements formatted as a string
'--ITC', type=float, Investment tax credit percentage (0-1)
'--FolderName', type=str, Makes the filename and directory
'--IR', type=float, Interest Rate
'--Lifetime', type=int, SMR book lifetime in years

And the optional arguements:

'--Min_Builds', type=str, Minimum build requirements (just use 0,0,0,0)
'--ITC_Duration', type=int, Years of ITC operation
'--NGP_Variable', type=str, Is natural gas price changing
'--NGP_V', type=str, Natural gas prices over the years
'--Max_Increase', type=float, Percent of module maximum growth per year
'--LearningRate', type=str, Learning Rate Scenario
'--Gen_File', type=str, Name of generator parameter file, default = 'LC_base_4pack.csv'
'--Init_N', type=str, Number of SMR modules previously built
'--CarbonTax', type=float, Carbon Tax in dollars per ton CO2
'--FacilityFile', type=str, File for the facilities', default = '2015_NG_National.csv'
'--Subsidy', type=float, Per SMR type subsidy quantity which will be distributed out to each until consumed, default = 0
'--CC', type=float, Capital Cost Modifier for sensitivites, default = 1

Example:
python ./V10_GL.py --Rand_n 925 --YearStart 2030 --YearCount 40 --HrCount 1 --NGP 10 --RollStep 20 --RollCount 20 --Min_Builds 0,0,0,0 --ITC 0.0 --ITC_Duration 3 --FolderName Git_example_10 --IR 0.07 --Lifetime 30 --Max_Increase 0.20 --Max_Builds 2,8,2,1

The output of this code will be: 
(folder) 'FolderName'largeoutputs: This is the per-heat range, per-facility resutls of natural gas supply and SMR output.
'FolderName'_ binaries_loggin_All.csv: This is the large collection of post-rolling horizon limitation Build, Construction, and Operation binaries.
'FolderName'_ binaries_loggin.csv: This is the reduced set of post-rolling horizon limitation Build binaries for selsected facilities.
profit.csv: This file tracks the best objective function for each facility for each rolling horizon iteration.
demand_selection.csv: This file tracks the selected temperature band for each facility for each rolling horizon iteration.

