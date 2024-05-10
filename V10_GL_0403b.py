## V5 This version is designed to use the per facility opt within the framework of a rolling horizon 

# Additionally, this will include the three temperature batching method per facility. 
# Therefore will add in a filter method in the results evaluation to select which of the three temperature ranges has the best profit. 
# could also do it by including a set in the optimizzation to select which one, but thi is additioanl opt wieght.

# V5b
# Also dropping the retirement binary as the window is sufficientlt small (5+5+30) For this go to V4
# Need to then filter out the already maxed out facilities in the filter process. 

## V6
# This update inputs the maximum per year installations. 

## V7
# Reincorporating the retirement binary. Testing if this is per facility or module.
from sys import base_exec_prefix
import pandas as pd
import numpy as np
#import cplex
import os

import matplotlib.pyplot as plt
import seaborn as sns
import math
import random

import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

import pyomo.environ as pe
import pyomo.opt
from pyomo.opt import SolverStatus

from datetime import datetime
s = datetime.now()
print(s,'start')

import argparse


parser = argparse.ArgumentParser(description='Import the input variables for the model run')
parser.add_argument('--Rand_n', type=int, help='Number of random facilities to sample',required=False)
parser.add_argument('--YearStart', type=int, help='Starting Year of analysis',required=True)
parser.add_argument('--YearCount', type=int, help='number of years to include',required=True)
parser.add_argument('--HrCount', type=int, help='number of hours per year',required=True)
parser.add_argument('--NGP', type = str, help='Natural Gas Price',required=True)
parser.add_argument('--RollStep', type=int, help='Rolling horizon step size',required=True)
parser.add_argument('--RollCount', type=int, help='Rolling horizon length',required=True)
parser.add_argument('--Min_Builds', type=str, help='Minimum build requirements',required=False)
parser.add_argument('--Max_Builds', type=str, help='Minimum build requirements',required=True)
parser.add_argument('--ITC', type=float, help='Investment tax credit',required=True)
parser.add_argument('--FolderName', type=str, help='Makes the filename and directory',required=True)
parser.add_argument('--IR', type=float, help='Interest Rate',required=True)
parser.add_argument('--Lifetime', type=int, help='SMR lifetime',required=True)
parser.add_argument('--ITC_Duration', type=int, help='Years of ITC operation',required=False)
parser.add_argument('--NGP_Variable', type=str, help='Is natural gas price changing',required=False)
parser.add_argument('--NGP_V', type=str, help='Natural gas prices over the years',required=False)
parser.add_argument('--Max_Increase', type=float, help='Percent of module maximum growth per year',required=False)
parser.add_argument('--LearningRate', type=str, help='Learning Rate Scenario',required=False)
parser.add_argument('--Gen_File', type=str, help='Name of generator parameter file',required=False)
parser.add_argument('--Init_N', type=str, help='Number of SMR modules previously built',required=False)
parser.add_argument('--CarbonTax', type=float, help='Carbon Tax in dollars per ton CO2',required=False)
parser.add_argument('--FacilityFile', type=str, help='File for the facilities',required=False)
parser.add_argument('--CC', type=float, help='Capital Cost Modifier',required=False, default = 1)

# Sample
#python ./V10_GL.py --Rand_n 1 --YearStart 2030 --YearCount 20 --HrCount 1 --NGP 50 --RollStep 1 --RollCount 3 --Min_Builds 0,0,0,0 --Max_Builds 1,1,1,1 --ITC 0.30 --FolderName LC_testing --IR 0.07 --Lifetime 15 --ITC_Duration 10 --Max_Increase 0.0

args = parser.parse_args()


try: 
        os.mkdir(args.FolderName) 
except OSError as error: 
    print('Directory Already Exists, Overwriting')

try: 
        os.mkdir(args.FolderName+'/'+args.FolderName+'_largeoutputs') 
except OSError as error: 
    print('Directory Already Exists, Overwriting')

if args.Gen_File is None:
    Gen_File = 'LC_base_4pack.csv'
else:
    Gen_File = args.Gen_File
print(Gen_File)

Gen_Count = pd.read_csv(Gen_File).shape[0]

print(pd.read_csv(Gen_File).shape[0])

if args.Init_N is None:
    Init_N = [0]*Gen_Count
    args.Init_N = 'NONE'
else:
    Init_N = [int(n) for n in args.Init_N.split(',')]
print('Starting model with '+ args.Init_N +'prebuilt modules (already learned)')

if args.FacilityFile is None:
    args.FacilityFile = '32_823_2015_NG_National.csv'
else:
    args.FacilityFile = args.FacilityFile
print(args.FacilityFile)


import Facility_Processing_2015 as FP

class EconomicAssessment_OneGen:
    def __init__(self):
        
        #### sets ####
        self.F = None                         # Facility index [0,1,2....]
        self.G = None                         # SMR type index [0,1,2....]
        self.Y = None                         # Year indicies [0,1,2....]
        self.T = None                         # Timestep indexes [0,1,2....8759]
        self.M = None                         # Maximum number of modules for any generator per facility [0,1,2... N]
        
        #### paramteres ####
        
        # Context paramters/System
        self.pNG_Cost = None                  # Cost of natural gas used throughout [f,y,t]
        self.pTLMP = None                     # NG price of heat $/MWht [f,y,t]
        self.pTempdemanded = None             # Minimum temperature which the reactor must supply for the given TLMP single value [f]
        self.pTEnergyDemand = None            # Thermal energy required per hour in MWt [f,t]
        self.ITC = args.ITC                   # Ivestment Taac Credit
        # Engineering parameters for generators
        self.pOnoroffinitial = 0              # initial state, default, 1 = On [f,g,m]
        self.pStartupfixedcost = None         # $/MW [g]
        self.pTVOC = None                     # Variable O&M ($/MWht) [g]
        self.pCap = None                      # Maximum Thermal capacity MWt [g]
        self.pThermaTransfereff = None        # V11 Thermal transsfer efficicnecy between coolant and workign fluid. [g]
        self.pMSL = None                      # Minimum stable thermal load (MWt) [g]
        self.pMDT = None                      # Minimum down time (hours) [g]
        self.pGenabovemininitial = 0          # Initial generation above the min (MWt) [g]
        self.pRGenInitial = 0                 # Initial generation of reactor (MWt) [g]
        self.pRamprate = None                 # Maximum ramp (up and down) (MWt/hr) [g]
        self.pAHF = None                      # Available heat fraction for NG input (temp dependant) (MWt/MWt)[f]
        self.pTDLoss = None                   # Transportation and distribution losses [f]

        self.pRollStep = None                 # number of year s for horizon to keep 

        
        self.pCAPEX_OS = None                 # Capital costs of the generator onsite $/MWe [g]
        self.pCAPEX_Mod = None                # Capital costs of the generator module $/MWe [g]
        self.pMod_Percent = None              # Degree of modularity (% factory module cost verusus total) [g]

        self.pFOPEX = None                    # Fixed O&M costs of the generator $/MWe [g]
        self.pModulesize = None               # Module size in MWe [g]
        self.pOutlettemp = None               # SMR outlet temperature [g]
        self.pWorkingTemp = None              # SMR working fluid temp [g]
        self.pMaxMod = None                   # maximum number of modules [g]
        
        self.pIR = args.IR                    # Discount/ Interest Rate 7% (OMB Circular) 
        self.pLifetime = args.Lifetime        # Lifetime, years
        self.pDR = self.pIR                   # Annual discount Rate value [y]
        self.pLR_OS = None                    # onsite learning rate [g]
        self.pLR_Mod = None                   # factory module learnign rate [g]
        self.pConstTime = None                # construction time (onsite and off) [g]
        self.pMod_N= None                     # Learned discount per module installed [N]


        #### variables ####
        
        # Float variables

        self.vGen = None                      # Heat generation (MWht) [f,g,y,m,t]
        self.vRGen = None                     # Heat generation (MWht) [f,g,y,m,t]
        self.vRGenabovemin = None             # Thermal Generation above MSL (MWht) [f,g,y,m,t]
        
        # Binary variables
        self.vTurnon = None                   # Does the generator turn on (1) or not (0) [f,g,y,m,t]
        self.vTurnoff = None                  # Gnerator turns off(1) or not (0) at time t [f,g,y,m,t]
        self.vOnoroff = None                  # Generator state at time t [f,g,y,m,t]
        self.vTOnoroff = None                 # Heat generation on (1) or off (0) [f,g,y,m,t]


        
        # results/objective variables  --These are not explicit model components and ar calculated from other variable values
        self.rProfit = None                   # Hourly profits ($) [m,t]
        
        self.vTotalProfits = None             # Total Profits ($) [-] single value
        self.rFeasible = None                 # is this configurations feasible     

        
        #V18 
        self.rACC = None                      # Annualized capital cost
        self.rFOMC = None                     # Fixed cost
        self.rVC = None                       # Variable cost
        self.rFC = None                       # Fuel cost
        self.rEP = None                       # Electricity profit
        self.rAVC = None                      # Avoided cost of natural gas
        self.rSUC = None                      # Startup costs
        self.RampCost = 1.0
        # LC_V1
        
        
    def BuildModel(self, linear= False):
        print('Building Model')

        model = pe.ConcreteModel()
        
        #### SETS ####     
        model.F = pe.Set(initialize = self.F, doc = "Facilities")
        model.G = pe.Set(initialize = self.G, doc = "SMR type")
        model.Y = pe.Set(initialize = self.Y, doc = "Years")
        model.T = pe.Set(initialize = self.T, doc = "Timestep")
        model.M = pe.Set(initialize = self.M, doc = "Modules")
        model.N = pe.Set(initialize = self.M, doc = "N for the Learning")
        
        
        #### PARAMETERS ####
        
        ## Market factors
        model.pTLMP = pe.Param(model.F,model.Y, model.T,initialize = {(self.F[z],self.Y[i],self.T[j]):self.pTLMP[z,i,j] for z in range(len(self.F)) for i in range(len(self.Y)) for j in range(len(self.T))} , doc = "NG Cost $/MWht")
        model.pCarbonTax = pe.Param(model.F, initialize = {self.F[i]:self.CarbonTax[i] for i in range(len(self.F))}, doc = "Carnon tax in dollars per MWh")

        ## Facility Details
        model.pTEnergyDemand = pe.Param(model.F, model.Y, model.T,initialize = {(self.F[f],self.Y[i],self.T[j]):self.pTEnergyDemand[f,i,j] for f in range(len(self.F)) for i in range(len(self.Y)) for j in range(len(self.T))}, doc = "Energy demand MWht")
        model.pAHF = pe.Param(model.F, initialize = {self.F[i]:self.pAHF[i] for i in range(len(self.F))}, doc = "Annual discount rate value")
        model.pTDLoss = pe.Param(model.F, initialize = {self.F[i]:self.pTDLoss[i] for i in range(len(self.F))}, doc = "Annual discount rate value")
        #model.pMinBuilds = pe.Param(model.Y, initialize = {self.Y[i]:self.pMinBuilds[i] for i in range(len(self.Y))}, doc = "Annual Minimum Builds")
        model.pPeak = pe.Param(model.F, initialize = {self.F[i]:self.pPeak[i] for i in range(len(self.F))}, doc = "Peak demand in test period")
        
        ## Cost Structure
        model.pDR = pe.Param(model.Y, initialize = {self.Y[i]:self.pDR[i] for i in range(len(self.Y))}, doc = "Annual discount rate value")
        model.pMod_pen = pe.Param(model.G,model.M, model.M,initialize = {(self.G[g],self.M[i],self.M[j]):self.pMod_pen[g,i,j] for g in range(len(self.G)) for i in range(len(self.M)) for j in range(len(self.M))}, doc = "Duplicate year build penalty")
        model.pCAPEX_OS = pe.Param(model.G,model.M, initialize = {(self.G[j],self.M[i]):self.pCAPEX_OS[j,i] for j in range(len(self.G)) for i in range(len(self.M))}, doc = "Onsite CApex")
        model.pCAPEX_Mod = pe.Param(model.G, initialize = {(self.G[j]):self.pCAPEX_Mod[j] for j in range(len(self.G))}, doc = "Onsite CApex")
        model.pStartupfixedcost = pe.Param(model.G, initialize = {self.G[i]:self.pStartupfixedcost[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pTVOC = pe.Param(model.G, initialize = {self.G[i]:self.pTVOC[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pFOPEX = pe.Param(model.G, initialize = {self.G[i]:self.pFOPEX[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        
        ## Enginerring Params
        model.pCap = pe.Param(model.G, initialize = {self.G[i]:self.pCap[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pThermaTransfereff = pe.Param(model.G, initialize = {self.G[i]:self.pThermaTransfereff[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pMSL = pe.Param(model.G, initialize = {self.G[i]:self.pMSL[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pMDT = pe.Param(model.G, initialize = {self.G[i]:self.pMDT[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pRGenInitial = pe.Param(model.G, initialize = {self.G[i]:self.pRGenInitial[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pRamprate = pe.Param(model.G, initialize = {self.G[i]:self.pRamprate[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pOutlettemp = pe.Param(model.G, initialize = {self.G[i]:self.pOutlettemp[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pWorkingTemp = pe.Param(model.G, initialize = {self.G[i]:self.pWorkingTemp[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pConstTime = pe.Param(model.G, initialize = {self.G[i]:self.pConstTime[i] for i in range(len(self.G))}, doc = "Annual discount rate value")
        model.pThermalSat = pe.Param(model.F,model.G, initialize = {(self.F[j],self.G[i]):self.pThermalSat[j,i] for j in self.F for i in self.G}, doc = "Thermal Quality sat.")
        
        ## Inherited Params from Prvious opt windows 
        model.pInitialB = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pInitialB[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Initial B value")
        model.pInitialC = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pInitialC[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Initial C value")
        model.pInitialOp = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pInitialOp[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Initial Op value")
        model.pInitialRet = pe.Param(model.F, model.M,initialize = {(self.F[f],self.M[m]):self.pInitialRet[f,m] for f in range(len(self.F)) for m in range(len(self.M))}, doc = "Initial Ret value")
        model.pInitialUndone = pe.Param(model.F, model.M,initialize = {(self.F[f],self.M[m]):self.pInitialUndone[f,m] for f in range(len(self.F)) for m in range(len(self.M))}, doc = "Initial Ret value")

        model.pCountB = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pCountB[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Count of B values")
        model.pCountC = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pCountC[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Count of C values")
        model.pCountOp = pe.Param(model.F,model.G, model.M,initialize = {(self.F[f],self.G[g],self.M[m]):self.pCountOp[f,g,m] for f in range(len(self.F)) for g in range(len(self.G)) for m in range(len(self.M))}, doc = "Count of Op values")
        
        #### VARIABLES ####
        
        ## SMR Binaries (Capacity)
        model.vB = pe.Var(model.F,model.G,model.Y,model.M,within = pe.Binary, doc = "Start building binary indicator")
        model.vB_done = pe.Var(model.F,model.G,model.Y,model.M,within = pe.Binary, doc = "End building binary indicator")
        model.vC = pe.Var(model.F,model.G,model.Y,model.M,within = pe.Binary, doc = "ongoing construction binary indicator")
        model.vOp = pe.Var(model.F,model.G,model.Y,model.M,within = pe.Binary, doc = "Operational year binary indicator")
        model.vGeneratorChoice = pe.Var(model.F,model.G, within = pe.Binary, doc = "which SMR type is being chosen")
        
        ## SMR operations
        model.vGenabovemininitial = pe.Var(model.F,model.G,model.Y,model.M,within = pe.NonNegativeReals, doc = "Gen above min of last hour of previous year")
        model.vRGenInitial = pe.Var(model.F,model.G,model.Y,model.M,within = pe.NonNegativeReals, doc = "Gen of last hour of previous year")
        model.vGen = pe.Var(model.F,model.G,model.Y,model.M,model.T,within=pe.NonNegativeReals, doc = " Heat generation (MWht) ")
        model.vRGen = pe.Var(model.F,model.G,model.Y,model.M,model.T,within=pe.NonNegativeReals, doc = " Total Reactor generation (MWht) ")
        model.vRGenabovemin = pe.Var(model.F,model.G,model.Y,model.M,model.T,within=pe.NonNegativeReals, doc = "Thermal Generation above MSL (MWht)")
        model.vRetire = pe.Var(model.F,model.Y,model.M,within = pe.Binary, initialize =0, doc = "0 if the unit is retired, 1 if operational or unbuilt.")
        model.vUndone = pe.Var(model.F,model.Y,model.M,within = pe.Binary, initialize =1, doc = "0 if the unit is retired, 1 if operational or unbuilt.")
        
        ## Natural gas
        model.vNG = pe.Var(model.F,model.Y,model.T, within = pe.NonNegativeReals, doc = "Hourly natural gas usage (MWht)")
        
        ## Cost variables - Learning
        model.vCAPEX_OS = pe.Var(model.F,model.G,model.M,within=pe.NonNegativeReals, doc = "Module onsite capital costs")
        model.vCAPEX_Mod = pe.Var(model.F,model.G,model.M,within=pe.NonNegativeReals, doc = "Module onsite capital costs")
        model.vCAPEX_OS_noFut = pe.Var(model.F,model.G,model.M,within=pe.NonNegativeReals, doc = "Module onsite capital costs")
        
        
        print('Params & Vars established')
        
        
        #### OBJECTIVES ####
        def CRF(model,g):
            return (self.pIR/(1-((1+self.pIR)**(-1*(self.pLifetime+model.pConstTime[g])))))
        
        def D(model):
            return (len(self.T)/8760)
        
        def TRev(model):
            return sum(sum(sum(sum(sum(model.vGen[f,g,y,m,t]*model.pTLMP[f,y,t]*model.pDR[y] + model.vGen[f,g,y,m,t]*model.pCarbonTax[f]*model.pDR[y]
                                       for t in model.T) for y in model.Y) for m in model.M) for g in model.G) for f in model.F)

        def SMR_Op_Costs(model):
            return sum(sum(sum(
                sum(sum(model.vRGen[f,g,y,m,t]*model.pTVOC[g] for t in model.T)*model.pDR[y] +
                     model.pFOPEX[g]*model.pCap[g]*model.vOp[f,g,y,m]*D(model)*model.pDR[y]  
                    for y in model.Y) for m in model.M) for g in model.G) for f in model.F)
        
        def SMR_Capex_Costs_OS(model):
            return (sum(sum(sum(sum((model.vC[f,g,y,m]+model.vOp[f,g,y,m])*model.pDR[y] for y in model.Y) * 
                        model.vCAPEX_OS[f,g,m]*CRF(model,g)*D(model)*(1-args.ITC)
                            for m in model.M) for g in model.G) for f in model.F))
        
        def SMR_Capex_Costs_Mod(model):
            return (sum(sum(sum(sum(
                    (model.vC[f,g,y,m]+model.vOp[f,g,y,m])*model.pCAPEX_Mod[g]*model.pDR[y]*D(model)*CRF(model,g)*(1-args.ITC)
                    for y in model.Y) for m in model.M) for g in model.G) for f in model.F))

        def Obj_Profit(model):
            return (TRev(model)-(SMR_Op_Costs(model)+ SMR_Capex_Costs_Mod(model)+SMR_Capex_Costs_OS(model)))
        
        
        model.Obj_Profit = pe.Objective(rule = Obj_Profit, sense=pe.maximize, doc = "Maximize the profits by balancing thermal and electric generation")
        

        #### CONSTRAINTS ####
        '''
        def MinBuild(model,y):
            return sum(sum(sum(model.vB[f,g,y,m] for f in model.F) for g in model.G) for m in model.M) >= model.pMinBuilds[y]
        model.MinBuild = pe.Constraint(model.Y, rule = MinBuild, doc = "Minimum Builds must occur per year")
        '''
        '''
        def MinBuild(model):
            return sum(sum(sum(sum(model.vB[f,g,y,m] for f in model.F) for g in model.G) for m in model.M) for y in model.Y)  >= 1
        model.MinBuild = pe.Constraint(rule = MinBuild, doc = "Minimum Builds must occur per year")
        '''
        def MaxBuildCap(model,f,g):
            return sum(sum(model.vB[f,g,y,m]for m in model.M) for y in model.Y)*model.pCap[g] <= model.pPeak[f]+model.pCap[g]
        model.MaxBuildCap = pe.Constraint(model.F,model.G, rule = MaxBuildCap, doc = "Limiting the number of modules onsite")
        
        def GenChoice(model,f,g,y,m):
            return model.vB[f,g,y,m]+model.pCountB[f,g,m] <= model.vGeneratorChoice[f,g]
        model.GenChoice = pe.Constraint(model.F, model.G,model.Y,model.M, rule = GenChoice, doc = "The status of the plant must be the previous state plus the turnon/off")
        
        def GenLimit(model,f):
            return sum(model.vGeneratorChoice[f,g] for g in model.G) == 1
        model.GenLimit = pe.Constraint(model.F, rule = GenLimit, doc = "The status of the plant must be the previous state plus the turnon/off")
        '''
        def CStatus(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vC[f,g,y,m] == (model.pInitialC[f,g,m] + model.pInitialB[f,g,m]- model.vB_done[f,g,y,m])
            return model.vC[f,g,y,m] == (model.vC[f,g,model.Y.prev(y),m] + model.vB[f,g,model.Y.prev(y),m]- model.vB_done[f,g,y,m])
        model.CStatus = pe.Constraint(model.F, model.G,model.Y,model.M, rule = CStatus, doc = "The status of the plant must be the previous state plus the turnon/off")
        '''
        def OS_LearningCosts(model,f,g,m):
            if args.LearningRate == 'NPP':
                return  model.vCAPEX_OS[f,g,m] == (((sum(sum(model.pMod_pen[g,m,n]*model.vB[f,g,y,n] for n in model.M)* model.vB[f,g,y,m] for y in model.Y)))* -1 * (self.pCAPEX_OS[g,m])) + self.pCAPEX_OS[g,m]
            return  model.vCAPEX_OS[f,g,m] ==  (((sum(sum(model.pMod_pen[g,m,n]*model.vB[f,g,y,n] for n in model.M)* model.vB[f,g,y,m] for y in model.Y))) * (self.pCAPEX_OS[g,m])) + self.pCAPEX_OS[g,m]
        model.OS_LearningCosts = pe.Constraint(model.F,model.G,model.M, rule = OS_LearningCosts, doc = "finds the learned onsite CAPEX for each module")
        
        def ThermalEquality(model,f,y,t):
            return sum(sum(model.vGen[f,g,y,m,t] for m in model.M)for g in model.G) + model.vNG[f,y,t] == model.pTEnergyDemand[f,y,t] ######## Need to evaluate whether the equaltiy is per genreator as well. ie all gens must equal the thermal req with gen choice binary
        model.ThermalEquality = pe.Constraint(model.F, model.Y,model.T, rule = ThermalEquality, doc = "Limits the thermal output (and therefore the 'revenue from heat')")

        def RGenEquality(model,f,g,y,m,t):
            return (model.vGen[f,g,y,m,t]) <= model.vRGen[f,g,y,m,t]*model.pThermaTransfereff[g]
        model.RGenEquality = pe.Constraint(model.F, model.G,model.Y,model.M,model.T,rule = RGenEquality, doc = "reactor output must be less than the nameplate capacity")
        
        def RStatus(model,f,g,y,m,t):
            return model.vRGen[f,g,y,m,t]<= model.pCap[g]*model.vOp[f,g,y,m]
        model.RStatus = pe.Constraint(model.F, model.G,model.Y,model.M,model.T,rule = RStatus, doc = "The sum of the electric and thermal generation must not surpass maximum generation")
        
        def TAllow(model,f,g,y,m,t):
            return model.vOp[f,g,y,m] <= model.pThermalSat[f,g]
        model.TAllow = pe.Constraint(model.F, model.G,model.Y,model.M,model.T, rule = TAllow, doc = "Must satisfty heat constraint")
        
        def RGenAbove(model,f,g,y,m,t):
            return model.vRGenabovemin[f,g,y,m,t]  == (model.vRGen[f,g,y,m,t])-(model.pMSL[g])*model.vOp[f,g,y,m]
        model.RGenAbove = pe.Constraint(model.F,model.G,model.Y,model.M,model.T, rule = RGenAbove, doc = "This defines the generation above the minimum stable load")
       
        
        def ConstLength(model,f,g,m):
            return sum(model.vC[f,g,y,m] for y in model.Y)+model.pCountC[f,g,m] == self.pConstTime[g]*(sum(model.vB[f,g,y,m] for y in model.Y)+model.pCountB[f,g,m])
        model.ConstLength = pe.Constraint(model.F,model.G,model.M, rule = ConstLength, doc = "Consturction time must eqaul the const binary sum")
        ''''''
        def MultiState(model,f,y,m):
            return model.vUndone[f,y,m] + sum(model.vOp[f,g,y,m] + model.vC[f,g,y,m]+ model.vB[f,g,y,m] for g in model.G) + model.vRetire[f,y,m] == 1
        model.MultiState = pe.Constraint(model.F, model.Y,model.M, rule = MultiState, doc = "Build, construction, and operation are mutually exclusive")
                
        def RetirementCont(model,f,y,m):
            if y == model.Y[1]:
                return model.vRetire[f,y,m] >= model.pInitialRet[f,m]
            return model.vRetire[f,y,m] >= model.vRetire[f,model.Y.prev(y),m]
        model.RetirementCont = pe.Constraint(model.F, model.Y,model.M, rule = RetirementCont, doc = "Cannot return/reretire")
        
        def PreCont(model,f,y,m):
            if y == model.Y[1]:
                return model.vUndone[f,y,m] <= model.pInitialUndone[f,m]
            return model.vUndone[f,y,m] <= model.vUndone[f,model.Y.prev(y),m]
        model.PreCont = pe.Constraint(model.F, model.Y,model.M, rule = PreCont, doc = "Cannot return/reretire")
        
        def LifetimeLimt(model,f,g,m):
            return sum(model.vOp[f,g,y,m] for y in model.Y)+model.pCountOp[f,g,m] <= self.pLifetime
        model.LifetimeLimt = pe.Constraint(model.F, model.G,model.M, rule = LifetimeLimt, doc = "Operation limited to lifetime")
        
        ''''''
        def BuildOrder(model,f,g,y,m):
            if m == model.M[1]:
                return pe.Constraint.Feasible
            return model.vB[f,g,y,m]+model.pCountB[f,g,m]  <= model.vB[f,g,y,model.M.prev(m)]+model.pCountB[f,g,model.M.prev(m)]
        model.BuildOrder = pe.Constraint(model.F,model.G,model.Y,model.M, rule = BuildOrder, doc = "Build from module 0 to 12 in order for clarity")
        
        def BuildLoophole(model,f,g,y,m):
            return model.vOp[f,g,y,m] <= sum(model.vB[f,g,y,m] for y in model.Y)+model.pCountB[f,g,m]
        model.BuildLoophole = pe.Constraint(model.F,model.G,model.Y,model.M, rule = BuildLoophole, doc = "No operations while building")
        ##########################
        '''
        def OperationStatus(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vOp[f,g,y,m] == (model.pInitialOp[f,g,m]+model.pInitialC[f,g,m])
            return model.vOp[f,g,y,m] <= model.vOp[f,g,model.Y.prev(y),m]+ model.vC[f,g,model.Y.prev(y),m] # -model.vRetire[f,y,m]
        model.OperationStatus = pe.Constraint(model.F,model.G,model.Y,model.M, rule = OperationStatus, doc = "Operation does not end")
        '''
        def ContinuousUse(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vC[f,g,y,m] <= model.pInitialC[f,g,m] + model.pInitialB[f,g,m] # model.vOp[f,g,y,m] + 
            return  model.vOp[f,g,y,m] + model.vC[f,g,y,m] >= model.vC[f,g,model.Y.prev(y),m]
        model.ContinuousUse = pe.Constraint(model.F, model.G,model.Y,model.M, rule = ContinuousUse, doc = "Cannot operate during construction")
        
        def ContinuousUse6(model,f,g,y,m):
            if y == model.Y[1]:
                return pe.Constraint.Feasible
            return  sum(model.vB[f,g,y,m]+model.vC[f,g,y,m] for g in model.G) + model.vUndone[f,y,m] <= sum(model.vB[f,g,model.Y.prev(y),m] + model.vC[f,g,model.Y.prev(y),m] for g in model.G) + model.vUndone[f,model.Y.prev(y),m]
        model.ContinuousUse6 = pe.Constraint(model.F, model.G,model.Y,model.M, rule = ContinuousUse6, doc = "Cannot operate during construction")
        '''
        def ContinuousUse5(model,f,g,y,m):
            if y == model.Y[1]:
                return pe.Constraint.Feasible
            return  sum(model.vB[f,g,y,m] for g in model.G) + model.vUndone[f,y,m] <= sum(model.vB[f,g,model.Y.prev(y),m] for g in model.G) + model.vUndone[f,model.Y.prev(y),m]
        model.ContinuousUse5 = pe.Constraint(model.F, model.G,model.Y,model.M, rule = ContinuousUse5, doc = "Cannot operate during construction")
        '''

        #####
        
        
        '''
        def ContinuousUse3(model,f,y,m):
            if y == model.Y[1]:
                return sum(model.vB[f,g,y,m] for g in model.G) + model.vUndone[f,y,m] >= model.pInitialUndone[f,m]
            return pe.Constraint.Feasible #sum(model.vB[f,g,y,m] for g in model.G) + model.vUndone[f,y,m]  >= sum(model.vB[f,g,model.Y.prev(y),m] for g in model.G)+ model.vUndone[f,model.Y.prev(y),m]
        model.ContinuousUse3 = pe.Constraint(model.F, model.Y,model.M, rule = ContinuousUse3, doc = "Cannot operate during construction")
        

        def OperationStatus(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vOp[f,g,y,m] == (model.pInitialOp[f,g,m]+model.pInitialC[f,g,m])
            return model.vOp[f,g,y,m] <= model.vOp[f,g,model.Y.prev(y),m]+ model.vC[f,g,model.Y.prev(y),m] # -model.vRetire[f,y,m]
        model.OperationStatus = pe.Constraint(model.F,model.G,model.Y,model.M, rule = OperationStatus, doc = "Operation does not end")
        
        
        def BuildLoophole2(model,f,g,m):
            return  sum(model.vB[f,g,y,m] for y in model.Y)+model.pCountB[f,g,m] <= 1
        model.BuildLoophole2 = pe.Constraint(model.F,model.G,model.M, rule = BuildLoophole2, doc = "Only one build")
        
        
        def ContinuousUse(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vC[f,g,y,m] <= model.pInitialC[f,g,m]+ model.pInitialB[f,g,m] # model.vOp[f,g,y,m] + 
            return  model.vOp[f,g,y,m] + model.vC[f,g,y,m] >= model.vC[f,g,model.Y.prev(y),m]
        model.ContinuousUse = pe.Constraint(model.F, model.G,model.Y,model.M, rule = ContinuousUse, doc = "Cannot operate during construction")
        
        
        def ContinuousUse4(model,f,y,m):
            if y == model.Y[1]:
                return sum(model.vC[f,g,y,m] for g in model.G) == sum(model.pInitialB[f,g,m] for g in model.G)
            return sum(model.vC[f,g,y,m] + model.vOp[f,g,y,m] for g in model.G) + model.vRetire[f,y,m] >= sum(model.vC[f,g,model.Y.prev(y),m] + model.vOp[f,g,model.Y.prev(y),m] for g in model.G) + model.vRetire[f,model.Y.prev(y),m]
        model.ContinuousUse4 = pe.Constraint(model.F, model.Y,model.M, rule = ContinuousUse4, doc = "Cannot operate during construction")
        '''
        ###################


        '''
        def BuildOrder2(model,f,g,m):
            return sum(model.vOp[f,g,y,m]for y in model.Y)+model.pCountC[f,g,m] >= sum(model.vB[f,g,y,m] for y in model.Y)+model.pCountB[f,g,m]
        model.BuildOrder2 = pe.Constraint(model.F,model.G,model.M, rule = BuildOrder2, doc = "Build from module 0 to 12 in order for learning")
        '''



        '''
        def ContinuousUse(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vOp[f,g,y,m] + model.vC[f,g,y,m]+ (model.vRetire[f,y,m]) >= model.pInitialOp[f,g,m]+model.pInitialC[f,g,m]+ (model.pInitialRet[f,m])
            return model.vOp[f,g,y,m] + model.vC[f,g,y,m]+ (model.vRetire[f,y,m]) >= model.vOp[f,g,model.Y.prev(y),m]+model.vC[f,g,model.Y.prev(y),m]+ (model.vRetire[f,model.Y.prev(y),m])
        model.ContinuousUse = pe.Constraint(model.F, model.G,model.Y,model.M, rule = ContinuousUse, doc = "Cannot operate during construction")
        '''
        '''
        def ContinuousUse2(model,f,y,m):
            if y == model.Y[1]:
                return sum(model.vOp[f,g,y,m] for g in model.G) + model.vRetire[f,y,m] == sum(model.pInitialOp[f,g,m] for g in model.G)
            return sum(model.vOp[f,g,y,m] for g in model.G) + model.vRetire[f,y,m] >= sum(model.vOp[f,g,model.Y.prev(y),m] for g in model.G)
        model.ContinuousUse2 = pe.Constraint(model.F, model.Y,model.M, rule = ContinuousUse2, doc = "Cannot operate during construction")
        '''
        '''
        def RetireLoop(model,f,m):
            return  sum(sum(model.vB[f,g,y,m] for y in model.Y)+model.pCountB[f,g,m] for g in model.G) >= sum(model.vRetire[f,y,m] for y in model.Y) + model.pInitialRet[f,m]
        model.RetireLoop = pe.Constraint(model.F,model.M, rule = RetireLoop, doc = "Only one build")
        
        
        def ContinuousOpHelper(model,f,g,y,m):
            if y == model.Y[1]:
                return pe.Constraint.Feasible #model.vOp[f,g,y,m]-model.pInitialOp[f,g,m]+model.pInitialRet[f,m] == 0 
            return model.vOp[f,g,y,m]-model.vOp[f,g,model.Y.prev(y),m]+model.vRetire[f,y,m] == 0
        model.ContinuousOpHelper = pe.Constraint(model.F, model.G, model.Y,model.M, rule = ContinuousOpHelper, doc = "Operation limited to lifetime")
        
        
        def RetirementCont(model,f,g,y,m):
            if y == model.Y[1]:
                return model.vRetire[f,y,m] <= model.pInitialRet[f,m]
            return model.vRetire[f,y,m] <= model.vRetire[f,model.Y.prev(y),m]
        model.RetirementCont = pe.Constraint(model.F, model.G,model.Y,model.M, rule = RetirementCont, doc = "Cannot return/reretire")
        
        def BuildLoophole3(model,f,g,y,m):
            return model.vB[f,g,y,m]+model.vOp[f,g,y,m] + model.vC[f,g,y,m] <= model.vRetire[f,y,m]
        model.BuildLoophole3 = pe.Constraint(model.F, model.G,model.Y,model.M, rule = BuildLoophole3, doc = "closing build loopholes")
        
        model.pInitialB.pprint()
        model.pInitialC.pprint()
        model.pInitialOp.pprint()
        model.pInitialRet.pprint()
        model.pCountB.pprint()
        model.pCountC.pprint()
        model.pCountOp.pprint()
        '''
        self.model = model
        print('Model Built')
        
        
    def SolveModel(self, solver='gurobi', linear = False):
        self.BuildModel(linear = linear)
        print('Solving...')
        #print(sorted(pyomo.IOptSolver._factory_cls.keys()))
        opt = pyomo.opt.SolverFactory(solver,tee = True)
        opt.options.mipgap = 0.025
        opt.options['timelimit'] = 2000

        #results = opt.solve(self.model, tee = False, logfile='CPLEX_LC_V1_test.log')
        results = opt.solve(self.model, tee = True, logfile=args.FolderName+'/log_LC_C_V5b.log')
        print('>>Solver status is {} and solver termination condition is {}'.format(results.solver.status,results.solver.termination_condition))
        
        if results.solver.termination_condition =='optimal':
            print(pe.value(self.model.Obj_Profit))
            self.vTotalProfits = pe.value(self.model.Obj_Profit)
            self.rFeasible = True
            self.ResultData()
            #self.ResultOutput()
        elif results.solver.termination_condition =='maxTimeLimit':
            print(pe.value(self.model.Obj_Profit))
            self.rFeasible = True
            self.ResultData()
            #self.ResultOutput()
            
        else:
            self.vRGen_pd = pd.DataFrame(np.zeros((2,2)))
            self.vNG_pd = pd.DataFrame(np.zeros((2,2)))
            self.rFeasible = False
            self.vTotalProfits= None
            self.rACC = 0
            self.rFOMC = 0
            self.rVC = 0
            self.rFC = 0
            self.rEP = 0
            self.rAVC = 0
            self.rSUC = 0
            self.ModCount = 0
            self.ResultData_inf()


            
    def ResultData_inf(self):
        
        PH = np.zeros((1,4,40,12,1))
        self.vGen = PH
        self.vRGen = PH
        self.vNG = np.zeros((1,40,1)) 
        self.vCAPEX_OS = np.zeros((12,4,1)) 
        
        self.vB = np.zeros((4,40,12,1)) 
        self.vC = np.zeros((4,40,12,1))
        self.vOp = np.zeros((4,40,12,1))
        self.vGeneratorChoice = np.zeros((4,1)) 
        self.vRet = np.zeros((40,12,1)) 
        self.vUnd = np.ones((40,12,1))

    def ResultData(self):
        

        self.vGen = np.array([[[[[pe.value(self.model.vGen[f,g,y,m,t]) for t in self.T] for m in self.M] for y in self.Y] for g in self.G] for f in self.F])
        self.vRGen = np.array([[[[[pe.value(self.model.vRGen[f,g,y,m,t]) for t in self.T] for m in self.M] for y in self.Y] for g in self.G] for f in self.F])
        self.vNG = np.array([[[pe.value(self.model.vNG[f,y,t]) for t in self.T] for y in self.Y] for f in self.F])
        self.vCAPEX_OS = np.array([[[pe.value(self.model.vCAPEX_OS[f,g,m])for m in self.M] for g in self.G] for f in self.F])
        self.vNG_pd = pd.DataFrame(np.sum(self.vNG,axis = 2))
        self.vRGen_pd = pd.DataFrame(np.sum(self.vRGen ,axis = (1,3,4)))
        print('OS_AFTER',self.vCAPEX_OS)
        
        self.vB = np.array([[[[round(pe.value(self.model.vB[f,g,y,m])) for y in self.Y] for m in self.M]  for g in self.G] for f in self.F])
        self.vC = np.array([[[[round(pe.value(self.model.vC[f,g,y,m])) for y in self.Y] for m in self.M]  for g in self.G] for f in self.F])
        self.vOp = np.array([[[[round(pe.value(self.model.vOp[f,g,y,m])) for y in self.Y] for m in self.M]  for g in self.G] for f in self.F])
        self.vGeneratorChoice = np.array([[round(pe.value(self.model.vGeneratorChoice[f,g]))for g in self.G] for f in self.F])
        self.rGenerator = [self.pGenerators[list(x).index(1)] for x in list(self.vGeneratorChoice)]
        self.vRet = np.array([[[round(pe.value(self.model.vRetire[f,y,m])) for y in self.Y] for m in self.M] for f in self.F])
        self.vUnd = np.array([[[round(pe.value(self.model.vUndone[f,y,m])) for y in self.Y] for m in self.M] for f in self.F])


        self.rB_Gen = [self.vB[x,list(self.vGeneratorChoice[x]).index(1),:,:] for x in range(len(self.vGeneratorChoice))]
        self.rC_Gen = [self.vC[x,list(self.vGeneratorChoice[x]).index(1),:,:] for x in range(len(self.vGeneratorChoice))]
        self.rOp_Gen = [self.vOp[x,list(self.vGeneratorChoice[x]).index(1),:,:] for x in range(len(self.vGeneratorChoice))]
        self.vGen_Gen = [self.vGen[x,list(self.vGeneratorChoice[x]).index(1),:,:,:] for x in range(len(self.vGeneratorChoice))]
        self.vRGen_Gen = [self.vRGen[x,list(self.vGeneratorChoice[x]).index(1),:,:,:] for x in range(len(self.vGeneratorChoice))]
        
        
        self.B_df = pd.DataFrame()
        self.C_df = pd.DataFrame()
        self.Op_df = pd.DataFrame()
        fi = 0
        while fi < len(self.rB_Gen):
            Gen = self.rGenerator[fi]
            F_num = self.fID[fi]
            Temp = self.pTempdemanded[fi]
            for i in range(len(self.rB_Gen[fi])):
                if self.rB_Gen[fi][i].sum()==0:
                    pass
                else:
                    self.B_df[('_').join((str(F_num),str(Temp),Gen,str(i)))] = self.rB_Gen[fi][i]
                    self.C_df[('_').join((str(F_num),str(Temp),Gen,str(i)))] = self.rC_Gen[fi][i]
                    self.Op_df[('_').join((str(F_num),str(Temp),Gen,str(i)))] = self.rOp_Gen[fi][i]
            fi+=1
        self.State_df = self.B_df+(self.C_df*2)+(self.Op_df*3)
        
        #print('B',self.vB)
        #print('NG', self.vNG)
        #print('C',self.vC)
        #print('Op',self.vOp)
        print('vGeneratorChoice', self.vGeneratorChoice)
        #print('CAPEX_OS',self.vCAPEX_OS)
        
        #self.vTotalProfits = pe.value(self.model.Obj_Profit)
        
        #self.rACC_OS = sum((self.vOp+self.vC)*self.vCAPEX_OS)*(len(self.T)/8760)*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))
        #print('ACC_OS', self.rACC_OS)8760
        #self.rACC_Mod = sum((self.vOp+self.vC)*self.pCAPEX_Mod*self.pCap)*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))
        #print('ACC_Mod', self.rACC_Mod)
        #self.rSUC = sum(self.vTurnon)*self.pStartupfixedcost
        #self.rFOMC = sum((self.vOp)*self.pFOPEX*self.pCap)*(len(self.T)/8760)*(self.pIR/(1-((1+self.pIR)**(-1*self.pLifetime))))
        #self.rVC =  sum(self.vRGen)*self.pVOC 
        #self.rFC =  sum(self.vRGen)*self.FuelCost

        
    def ResultOutput(self):
        self.vRGen_PD = pd.DataFrame()
        self.vGen_PD = pd.DataFrame()
        self.vTurnon_PD = pd.DataFrame()
        self.vOnoroff_PD = pd.DataFrame()
        
        self.Output = pd.DataFrame()

        m = list(self.M)
        for x in m:
            self.Output['Reactor_Gen [MWht]_'+str(x)] = self.vRGen_PD[x]
            self.Output['Thermal_Gen [MWht]_'+str(x)] = self.vGen_PD[x]
        self.Output['SOC'] = self.vTESSOC[:,self.vGenerator.index(1)]
        self.Output

def ParamsVarsPD(base,sites,F_pd,NG_Lmps,startHr, 
    hrcount, years, startYr, Conditions = None, FacN = None, Step = 1, n2 = None, iterationN = 0):
    
    sites.reset_index(drop = True, inplace = True)
    F_pd.reset_index(drop = True, inplace = True)
    
    #### Sets ####
    
    base.F = F_pd.index.tolist()
    base.G = sites.index.tolist()
    base.Y = list(range(years))
    base.M = list(range(12))
    base.T = list(range(hrcount))

    base.pTlimit = F_pd['Temp_Req'].tolist()   

    ## SMR facility pairing [aindexed on f and g]
    base.pThermalSat = []
    base.pWorkingTemp = sites['Working Fluid Temp (C)'].tolist()

    TS = True
    for f in base.F:
        for g in base.pWorkingTemp:
            if g<=base.pTlimit[f]:
                TS = False
    #if not TS:
    #base.M = list(range(1))
    base.TS = TS

    for f in base.F:
        temp_l = []
        for g in base.G:
            temp_l.append(int(base.pWorkingTemp[g]>=base.pTlimit[f]))
            
        base.pThermalSat.append(temp_l)
    base.pThermalSat = np.array(base.pThermalSat)
    
    ## Facility Params [f indexed only]
    base.pTempdemanded = F_pd['Temp_degC'].tolist()    
    print(int(base.pTempdemanded[0]), ' Deg C')
    base.pAHF = [(-0.00038*(f) + 0.90556) for f in base.pTempdemanded]
    base.pTDLoss = [0.95]* len(base.F)
    
    base.fID = F_pd['FACILITY_ID'].tolist()
    base.fNAICS = F_pd['MECS_NAICS'].tolist()
    base.fState = F_pd['STATE'].tolist()
    base.fIndustry = F_pd['Industry'].tolist()
    base.fAvgDem = F_pd['Thermal MWh/hr'].tolist()
    
    base.pTEnergyDemand = []
    base.fhrProfiles = []
    base.pPeak = []
    temp_l = []
    
    f = 0
    for n in base.fNAICS:
        temp_l = []
        hrProfiles = []
        for yi in range(years):
            hrP = FP.QuickProfile(startYr+yi, n)
            
            hrProfiles.append(hrP)
            TotalHeat = base.fAvgDem[f]*hrcount
            TotalHours = sum(hrP[startHr:startHr+hrcount])
            CorrectiveFactor = TotalHeat/TotalHours
            NormalizedHour = [hrP[startHr+i] * CorrectiveFactor for i in range(hrcount)]
            temp_l.append(NormalizedHour)
        base.pPeak.append(max(NormalizedHour))
        base.fhrProfiles.append(hrProfiles)
        base.pTEnergyDemand.append(temp_l)
        f+=1
    print(base.pPeak[0],'Peak MW')
    #base.fhrProfiles = np.array(base.fhrProfiles)
    base.pTEnergyDemand = np.array(base.pTEnergyDemand)
    
    ## Market conditions [indexed on f, y, t]
    if args.NGP_Variable == 'True':
        base.pTLMP = []
        NGP_V_list = [float(ngp) for ngp in args.NGP_V.split(',')]
        if len(NGP_V_list) == args.RollStep*args.RollCount+args.YearStart:
            for f in base.F:
                temp_l = []
                for y in base.Y:
                    temp_l.append([NGTempCostCurve(base.pTempdemanded[f],NG_Cost = NGP_V_list[y+(iterationN*args.RollStep)])]*len(base.T))
                base.pTLMP.append(temp_l)
            base.pTLMP = np.array(base.pTLMP)
        else:
            for f in base.F:
                temp_l = []
                for y in base.Y:
                    if y >= len(NGP_V_list):
                        temp_l.append([NGTempCostCurve(base.pTempdemanded[f],NG_Cost = float(NG_Lmps))]*len(base.T))
                    else:
                        temp_l.append([NGTempCostCurve(base.pTempdemanded[f],NG_Cost = NGP_V_list[y])]*len(base.T))
                base.pTLMP.append(temp_l)
            base.pTLMP = np.array(base.pTLMP)
        print('NGPs',base.pTLMP)

    else:
        base.pTLMP = []
        for f in base.F:
            temp_l = []
            for y in base.Y:
                temp_l.append([NGTempCostCurve(base.pTempdemanded[f],NG_Cost = float(NG_Lmps))]*len(base.T))
            base.pTLMP.append(temp_l)
        base.pTLMP = np.array(base.pTLMP)
    
    if args.CarbonTax is None:
        base.CarbonTax = [0]*len(base.F)
    else:
        base.CarbonTax = []
        for f in base.pAHF:
            CT = (3.412 * (1/f) * 52.91 / 1000) * args.CarbonTax  # 3.412 - MMBtu to MWh, 52.91 kg CO2/mmbtu EPA
            base.CarbonTax.append(CT)
    print('Carbon Tax $/MWh',base.CarbonTax)
    ## SMR market conditions
    MB = [int(mbs) for  mbs in args.Max_Builds.split(',')]
    base.pMaxMod_annual = MB
    if args.LearningRate == 'NPP':
        base.pLR_Mod = (1-sites['Mod Learning_Low']).tolist() 
        base.pLR_OS = (1-sites['OS Learning_Low']).tolist()
    elif args.LearningRate == 'opt':
        base.pLR_Mod = (1-sites['Mod Learning_High']).tolist()
        base.pLR_OS = (1-sites['OS Learning_High']).tolist()
    elif args.LearningRate == 'Zero':
        base.pLR_Mod = [1]*len(base.G)
        base.pLR_OS = [1]*len(base.G)
    elif args.LearningRate == 'NewMid':
        base.pLR_Mod = [0.84]*len(base.G)
        base.pLR_OS = [1]*len(base.G)
    elif args.LearningRate == 'OnsiteOnly':
        base.pLR_Mod = [1]*len(base.G)
        base.pLR_OS = [0.95]*len(base.G)
    else:
        base.pLR_Mod = (1-sites['Mod Learning_Med']).tolist()
        base.pLR_OS = (1-sites['OS Learning_Med']).tolist()
    
    #print(base.pLR_Mod)
    ## SMR engineering parameters [indexed on g]
    base.pGenerators = sites['Sites'].tolist()
    base.pCap = sites['Power in MWt'].tolist()
    base.pMSL = sites['MSL in MWt'].tolist()
    base.pThermaleff = sites['Thermal Efficiency'].tolist()
    base.pThermaTransfereff = sites['Thermal Transfer Efficiency'].tolist()
    base.pRGenInitial = [0]*len(base.G)
    base.pRamprate = sites['Ramp Rate (MW/hr)'].tolist()
    base.pOutlettemp = sites['Outlet Temp (C)'].tolist()
    base.pWorkingTemp = sites['Working Fluid Temp (C)'].tolist()
    base.pConstTime = [3]*len(base.G)
    base.pMDT = sites['MDT in hours'].tolist()
    base.pOnoroffinitial = [0]*len(base.G)
    base.pGenabovemininitial = [0]*len(base.G)
    
    ## SMR facility pairing [aindexed on f and g]
    base.pThermalSat = []
    for f in base.F:
        temp_l = []
        for g in base.G:
            temp_l.append(int(base.pWorkingTemp[g]>=base.pTlimit[f]))
        base.pThermalSat.append(temp_l)
    base.pThermalSat = np.array(base.pThermalSat)
    
    
    #### initial Conditions Based on the repeating
    
    ## Build Binary [indexed on f, g, y, m]
    base.pInitialB = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))
    base.pInitialC = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))
    base.pInitialOp = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))
    base.pInitialRet = np.array([[0] * len(base.M)]*len(base.F))
    base.pInitialUndone = np.array([[1] * len(base.M)]*len(base.F))
    base.pCountB = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))
    base.pCountC = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))
    base.pCountOp = np.array([[[0]*len(base.M)]*len(base.G)]*len(base.F))

    
    base.pModConstructed  = [0]*len(base.G)

    if Conditions != None:
        VB_l, VC_l, VOP_l, VB_c, VC_c,VOP_c, VRET_l, VUND_l = Conditions

        
        base.pInitialB = np.array([VB_l[n2][:,:,Step-1]])
        
        base.pInitialC = np.array([VC_l[n2][:,:,Step-1]])
        base.pInitialOp = np.array([VOP_l[n2][:,:,Step-1]])
        base.pInitialRet = np.array([VRET_l[n2][:,Step-1]]) ########################################
        base.pInitialUndone = np.array([VUND_l[n2][:,Step-1]]) ########################################
        base.pCountB = np.array([VB_c[n2]])
        base.pCountC = np.array([VC_c[n2]])*np.array([VB_c[n2]])
        base.pCountOp = np.array([VOP_c[n2]])*np.array([VB_c[n2]])

        print('VUND',base.pInitialUndone)
        print('IB',base.pInitialB)
        
        print('IC',base.pInitialC )
        print('IOP',base.pInitialOp )
        print('IRET',base.pInitialRet)
        print('IU',base.pInitialUndone)
        print('CB',base.pCountB )
        print('CC',base.pCountC)
        print('COP',base.pCountOp )

        base.pModConstructed  = list(np.sum(VB_c,axis = (0,2)))
        print(base.pModConstructed)

    
    ## SMR Costs [indexed on g]
    base.pVOC = sites['VOM in $/MWh-e']*sites['Thermal Efficiency']
    base.FuelCost = sites['FC in $/MWh-e']*sites['Thermal Efficiency']
    base.pTVOC = (base.pVOC+base.FuelCost).tolist()
    base.pFOPEX = (sites['FOPEX $/kWt']*1000).tolist()
    base.pStartupfixedcost = sites['Startupfixedcost in $'].tolist()
    base.pMod_Percent = sites['Modularity']
    base.pCAPEX_OS_base = (sites['CAPEX $/kWt']*1000*args.CC)*(1-sites['Modularity']).tolist()
    base.pCAPEX_Mod_base = ((sites['CAPEX $/kWt']*1000*args.CC)*sites['Modularity']).tolist()
    
    print(base.pCAPEX_OS_base)
    
    ## SMR costs by module at each facility (onsite learning rates)
    base.pCAPEX_OS = np.array([[base.pCAPEX_OS_base[g]*base.pCap[g]*((m+1)**(math.log(base.pLR_OS[g])/math.log(2))) for m in base.M] for g in base.G])
    totalN = len(base.M)*len(base.F)
    base.pCAPEX_Mod = np.array([base.pCAPEX_Mod_base[g]*base.pCap[g]*((base.pModConstructed[g]+Init_N[g]+1)**(math.log(base.pLR_Mod[g])/math.log(2))) for g in base.G])
    
    print('OS', base.pCAPEX_OS)
    print('Mod', base.pCAPEX_Mod)

    ## Learning rate Oniste 
    temp_l = []
    for g in base.G:
        penalties = [(x+1)**(math.log(base.pLR_OS[g])/math.log(2)) for x in range(len(base.M))]
        penalties_d = [penalties[x]-penalties[x+1] for x in range(len(penalties)-1)]+[0]
        P = [[0]*len(penalties)]
        for i in range(len(base.M)-1):
            mask = [0]*len(penalties)
            n = 0
            for v in penalties_d:
                if n <= i:
                    mask[n] = v
                n+=1
            P.append(mask)
        temp_l.append(P)
    base.pMod_pen = np.array(temp_l)
    
    ## General cost [indexed on y]
    base.pDR = [1/((1+base.pIR)**(y)) for y in range(years)] #[1]* years #
    '''
    mins = [0]*len(base.Y)

    for l in base.Y:
        if l<len((args.Min_Builds).split(',')):
            mins[l]+= int((args.Min_Builds).split(',')[l])


    base.pMinBuilds = mins
    '''


def ImportData(sites, Facilities):
    Sites_pd = pd.read_csv(sites) 
    F_pd = pd.read_csv(Facilities) 
    return Sites_pd, F_pd

def NGTempCostCurve(Temp,NG_Cost = 3.0, AHF_Coeffs = [0,-0.00038,0.90556]):
    HHV = 54 # MJ/kg
    Density = 0.68 # kg/m3
    cfTom3 = 35.31 # Unit conversion
    AHF = AHF_Coeffs[0]*(int(Temp)^2) + AHF_Coeffs[1]*int(Temp) + AHF_Coeffs[2] # avaialble Heat fraction - Deep Patel Equation
    
    HHV = HHV*Density*(1/cfTom3)*(1/1000000)*(1000) # returns  TJ/thousand cf 
    Cost = NG_Cost*(1/HHV)*(1/AHF)*(1/277.778) # returns the Cost in $/MWh
    return Cost

def mainFullOpt(NGCost = 4, hrcount = 8760, dem = 100,sites = Gen_File, Facilities = None, startyear = 2020,
                startHr = 0, years = 5, linear = False, FacN = None, Initial = True,Step = 1, Conditions = None, n2 = None, iterationN = None):
    if Initial:
        s = datetime.now()
        #S, F = ImportData(sites, Facilities)
        S = pd.read_csv(sites) 
        F = Facilities

        if FacN != None:
            F = pd.DataFrame(F.iloc[FacN]).T.reset_index(drop= True)
        Assessment = EconomicAssessment_OneGen()

        ParamsVarsPD(Assessment,S,F,NGCost,startHr, hrcount, years,startyear, iterationN = iterationN)
        Assessment.SolveModel(linear = linear)
    else:
        s = datetime.now()
        S = pd.read_csv(sites) 
        F = Facilities

        if FacN != None:
            F = pd.DataFrame(F.iloc[FacN]).T.reset_index(drop= True)
            
        Assessment = EconomicAssessment_OneGen()

        ParamsVarsPD(Assessment,S,F,NGCost,startHr, hrcount, years,startyear, Conditions = Conditions, FacN = FacN, Step = Step, n2 = n2, iterationN = iterationN)
        Assessment.SolveModel(linear = linear)
    
    return Assessment

def FacilityImportProcess():
    print(os.getcwd())
    FacilitiesFile = args.FacilityFile # '2015_NG_National.csv' #'32_823_2015_NG_National.csv' # '2015_NG_National_trunc.csv' #
    Facs = pd.read_csv(FacilitiesFile)
    Facs

    Facs_base = Facs[['CITY', 'COUNTY', 'COUNTY_FIPS', 'END_USE', 'FACILITY_ID',
           'FINAL_NAICS_CODE', 'MECS_NAICS', 'REPORTING_YEAR', 'STATE', 
           'UNIT_NAME', 'UNIT_TYPE', 'Lat', 'Lon', 'LLType',
           'geometry', 'FacilityName', 'Industry', 'Company']]

    #  'Temp_degC', 'Thermal MWh/hr'

    Facs['3-5'] = True
    Facs['5-8'] = True
    Facs.loc[Facs['Temp 300']-Facs['Temp 550'] == 0, '3-5'] = False
    Facs.loc[Facs['Temp 550']-Facs['Temp 850'] == 0, '5-8'] = False
    Facs

    Facs_300 = Facs_base.copy()
    Facs_550 = Facs_base.copy()
    Facs_850 = Facs_base.copy()

    Facs_300['Thermal MWh/hr'] = Facs['Demand 300']*1.1*277.778/8760
    Facs_300['Temp_degC'] = Facs['Temp 300']
    Facs_300['Temp_Req'] = 300
    Facs_300['Indexer'] = range(len(Facs_300.index))

    Facs_550['Thermal MWh/hr'] = Facs['Demand 550']*1.1*277.778/8760
    Facs_550['Temp_degC'] = Facs['Temp 550']
    Facs_550['Temp_Req'] = 550
    Facs_550['Indexer'] = range(len(Facs_550.index))
    Facs_550 = Facs_550.loc[Facs['3-5']]

    Facs_850['Thermal MWh/hr'] = Facs['Demand 850']*1.1*277.778/8760
    Facs_850['Temp_degC'] = Facs['Temp 850']
    Facs_850['Temp_Req'] = 850
    Facs_850['Indexer'] = range(len(Facs_850.index))
    Facs_850 = Facs_850.loc[Facs['5-8']]

    Facs_total = pd.concat((Facs_300,Facs_550,Facs_850)).reset_index(drop = True)

    Facs_total['FullIndex'] = range(len(Facs_total.index))

    F_IDS = Facs.FACILITY_ID.unique().tolist()
    Inds = []
    for f in F_IDS:
        Inds.append(Facs_total.loc[Facs_total['FACILITY_ID']==f,'FullIndex'].tolist())
    Fac_dict = {F_IDS[x]:Inds[x] for x in range(len(F_IDS))}

    return Facs_total, Fac_dict

def graphAll(VB_graph,VC_graph,VOP_graph, VB_graph_sum,VC_graph_sum,VOP_graph_sum):

    VB_graph.reset_index(drop = True,inplace = True)
    VC_graph.reset_index(drop = True,inplace = True)
    VOP_graph.reset_index(drop = True,inplace = True)
    VB_graph_sum.reset_index(drop = True,inplace = True)
    VC_graph_sum.reset_index(drop = True,inplace = True)
    VOP_graph_sum.reset_index(drop = True,inplace = True)

    ### Output the generator sums

    sns.set(rc = {'figure.figsize':(8,4)})
    sns.set_style("whitegrid")
    g7 = sns.lineplot(data = VB_graph_sum.reset_index(drop = True),x = 'Year', y = 'Value', hue = 'Generator', style = 'Iteration', palette = 'tab20')

    g7.set(ylabel='Modules Built')
    g7.set(xlabel='Year')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.show()
    plt.tight_layout()
    g7.figure.savefig(args.FolderName + '/'+args.NGP+'_'+str(args.YearCount)+'_Learning_B_20_Summer.png', dpi = 700, bbox_inches = "tight")
    plt.cla()
    plt.clf()


    sns.set(rc = {'figure.figsize':(8,4)})
    sns.set_style("whitegrid")
    g8 = sns.lineplot(data = VC_graph_sum.reset_index(drop = True),x = 'Year', y = 'Value', hue = 'Generator', style = 'Iteration', palette = 'tab20')

    g8.set(ylabel='Modules Under Construction')
    g8.set(xlabel='Year')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.show()
    plt.tight_layout()
    g8.figure.savefig(args.FolderName + '/'+args.NGP+'_'+str(args.YearCount)+'_Learning_C_20_Summed.png', dpi = 700, bbox_inches = "tight")
    plt.cla()
    plt.clf()


    sns.set(rc = {'figure.figsize':(8,4)})
    sns.set_style("whitegrid")
    g9 = sns.lineplot(data = VOP_graph_sum.reset_index(drop = True),x = 'Year', y = 'Value', hue = 'Generator', style = 'Iteration', palette = 'tab20')

    g9.set(ylabel='Modules Operating')
    g9.set(xlabel='Year')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.show()
    plt.tight_layout()
    g9.figure.savefig(args.FolderName + '/'+args.NGP+'_'+str(args.YearCount)+'_Learning_Op_20_Summed.png', dpi = 700, bbox_inches = "tight")
    plt.cla()
    plt.clf()

Facs_total, Fac_dict = FacilityImportProcess()


rand_n = args.Rand_n

def Optimizations(Facs_total, Fac_dict, rand_n = args.Rand_n):

    SampleSet = list(Fac_dict.keys())[:rand_n] 

    Profits = []
    VB_l = []
    VC_l = []
    VOP_l = []
    VRET_l = []
    VG_l = []

    Demand_pd = pd.DataFrame()
    Profit_pd = pd.DataFrame()
    B_pd = pd.DataFrame()



    ### Start
    n = 0

    Pre_Profits = []
    Pre_Max_Porfs = []
    Pre_VB_l = []
    Pre_VC_l = []
    Pre_VOP_l = []
    Pre_VRET_l = []
    Pre_VUND_l = []
    Pre_VG_l = []
    Pre_G_ind = []
    temp_DemN = []
    temp_Prof = []

    for s in SampleSet:
        print('#######################', n, s)
        temp_Profits = []
        temp_VB_l = []
        temp_VC_l = []
        temp_VOP_l = []
        temp_VG_l = []
        temp_VRET_l = []
        temp_VUND_l = []
        built = []

        temp_NG = []
        

        for i in Fac_dict[s]:
            FR = mainFullOpt(NGCost = args.NGP, hrcount = args.HrCount, sites = Gen_File, 
                             startyear = args.YearStart,  startHr = 0, years = args.YearCount,
                             Facilities = Facs_total, FacN = i, iterationN = 0)
            try:
                float(FR.vTotalProfits)
            except:
                FR.vTotalProfits = -1000000000000

            temp_Profits.append(FR.vTotalProfits)
            #FR.vNG_pd.to_csv(args.FolderName+'/'+args.FolderName+'_largeoutputs/'+str(s)+'_'+str(i)+'_0_vNG.csv')
            #FR.vRGen_pd.to_csv(args.FolderName+'/'+args.FolderName+'_largeoutputs/'+str(s)+'_'+str(i)+'_0_vRGen.csv')

            if FR.vTotalProfits > 1:
                built.append(1)
            else:
                built.append(0)
            temp_VB_l.append(FR.vB[0])
            temp_VC_l.append(FR.vC[0])
            temp_VOP_l.append(FR.vOp[0])
            temp_VG_l.append(FR.vGeneratorChoice[0])
            temp_VRET_l.append(FR.vRet[0])
            temp_VUND_l.append(FR.vUnd[0])

        print(built)
        if sum(built)>=1:
            MaxProf_Dem = temp_Profits.index(max(temp_Profits))
            temp_DemN.append(MaxProf_Dem)
            temp_Prof.append(temp_Profits[MaxProf_Dem])
            Pre_Profits.append(temp_Profits[MaxProf_Dem])
            Pre_VB_l.append(temp_VB_l[MaxProf_Dem])
            Pre_VC_l.append(temp_VC_l[MaxProf_Dem])
            Pre_VOP_l.append(temp_VOP_l[MaxProf_Dem])
            Pre_VG_l.append(temp_VG_l[MaxProf_Dem])
            Pre_Max_Porfs.append(temp_Profits.index(max(temp_Profits)))
            Pre_G_ind.append(list(temp_VG_l[MaxProf_Dem]).index(1))
            Pre_VRET_l.append(temp_VRET_l[MaxProf_Dem])
            Pre_VUND_l.append(temp_VUND_l[MaxProf_Dem])
            Fac_dict[s] = [Fac_dict[s][MaxProf_Dem]] ###########
        else:
            temp_DemN.append(4)
            temp_Prof.append(0)
            Pre_Profits.append(temp_Profits[0])
            Pre_VB_l.append(temp_VB_l[0])
            Pre_VC_l.append(temp_VC_l[0])
            Pre_VOP_l.append(temp_VOP_l[0])
            Pre_VG_l.append(temp_VG_l[0])
            Pre_VRET_l.append(temp_VRET_l[0])
            Pre_VUND_l.append(temp_VUND_l[0])
            Pre_Max_Porfs.append(0)
            Pre_G_ind.append(8)
        n+=1

    Demand_pd['0'] = temp_DemN
    Profit_pd['0'] = temp_Prof
    


    limit = [int(mbs) for  mbs in args.Max_Builds.split(',')]
    toBuilt_indexes =[]
    for gi in range(Gen_Count):
        mask_gi = np.array([gi == si for si in Pre_G_ind])
        sorted_profs = np.argsort(mask_gi*np.array(Pre_Profits))[::-1]
        built_n = sum(mask_gi)
        if built_n > int(limit[gi]):
            toBuild = int(limit[gi])
        else:
            toBuild = built_n
        toBuilt_indexes+=list(sorted_profs[:toBuild])

    toBuilt_indexes

    built = []
    gen_ind = []
    mod_count = []

    fc = 0
    for fb in Pre_VB_l[:len(SampleSet)]:
        mb = np.sum(fb, axis = (1,2))
        mi = np.sum(fb, axis = (0,2))
        built.append(mi)
        mod_count.append(mb.sum())
        mg = np.clip(mb, 0,1)
        if sum(mg)>0:
            gen_ind.append(list(mg).index(1))
        else:
            gen_ind.append(10)


    toBuilt_indexes_f =[]
    toBuilt_index_m = []

    for gi in range(Gen_Count):
        mask_gi = np.array([gi == gt for gt in gen_ind])

        built_g = sum(mask_gi)
        built_n = sum(mask_gi*np.array(mod_count))
        if built_n > int(limit[gi]):
            toBuild = int(limit[gi])
        else:
            toBuild = built_n
        print(built_n, toBuild)

        sorted_profs = np.argsort(mask_gi*np.array(Pre_Profits))[::-1]
        b_h = 0
        for sp in sorted_profs[:built_g]:

            # filter for the first year builds only
            # If no first year builds, then paste the 'if b_h...' section for no build

            # If first years, then create mask of first year build, multiply by v_3d_2.
            ## Sum up the first year and add that instead of toBuild-b_h in second clause
            years = Pre_VB_l[sp].shape[2]
            mods = Pre_VB_l[sp].shape[1]
            gens = Pre_VB_l[sp].shape[0]

            if (Pre_VB_l[sp][gi][:,0].sum() == 0):
                print('No first year builds',Pre_VB_l[sp][gi])

                Pre_VB_l[sp][gi]*=0
                years = Pre_VB_l[sp].shape[2]
                v2 = np.array(np.array([0]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h))))
                v2 = v2[:mods]
                v_3d_2 = np.array([[m]*years for m in v2])
                Pre_VC_l[sp][gi]*=v_3d_2
                Pre_VOP_l[sp][gi]*=v_3d_2
                temp_und = Pre_VUND_l[sp]
                Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                Pre_VRET_l[sp]*=v_3d_2
                print('no build_FY0', v2.sum())

                b_h+=0
                
            else:
                    

                FY0_build = Pre_VB_l[sp][gi][:,0].sum()
                print('FY0',FY0_build)

                v = np.array(np.array([1]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h))))
                v2 = np.array(np.array([1]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h))))
                    
                v_FY0 = np.array(np.array([1]*(FY0_build)+[0]*(len(built[sp])-(FY0_build))))
                v2_FY0 = np.array(np.array([1]*(FY0_build)+[0]*(len(built[sp])-(FY0_build))))
                    
                print('v',v)

                v_mask = v[:mods]* v_FY0
                v2_mask = v2[:mods]* v2_FY0

                nB = v_mask.sum()


                if b_h >= toBuild:
                    Pre_VB_l[sp][gi]*=0
                    years = Pre_VB_l[sp].shape[2]
                    vn = np.array(np.array([0]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h))))
                    v_3d_2 = np.array([[m]*years for m in vn])
                    Pre_VC_l[sp][gi]*=v_3d_2
                    Pre_VOP_l[sp][gi]*=v_3d_2
                    temp_und = Pre_VUND_l[sp]
                    Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                    Pre_VRET_l[sp]*=v_3d_2
                    print('no build', vn.sum())

                    b_h+=0

                elif sum(built[sp])>toBuild-b_h:

                    try:
                        v_3d = np.array([[m]*years for m in v_mask])
                        Pre_VB_l[sp][gi]*=v_3d
                    except:
                        v = v[:mods]
                        v_3d = np.array([[m]*years for m in v_mask])
                        Pre_VB_l[sp][gi]*=v_3d
                        print('EXCEPTION IN SIZES HIT, ln 1555')
                    print(v_3d)
                    
                    v2_mask = v2_mask[:mods]
                    v_3d_2 = np.array([[m]*years for m in v2_mask])
                    Pre_VC_l[sp][gi]*=v_3d_2
                    Pre_VOP_l[sp][gi]*=v_3d_2
                    temp_und = Pre_VUND_l[sp]
                    Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                    
                    Pre_VRET_l[sp]*=v_3d_2
                    b_h+=nB ###################################
                    print('partial build',nB,'to', b_h)
                else:
                    b_h+=sum(built[sp])
                    print('full build', sum(built[sp]),'to', b_h)
        print('Built per gen',b_h)

    VB_l = np.stack(Pre_VB_l[:len(SampleSet)])

    VC_l = np.stack(Pre_VC_l[:len(SampleSet)])
    VOP_l = np.stack(Pre_VOP_l[:len(SampleSet)])
    VG_l = np.stack(Pre_VG_l[:len(SampleSet)])
    VRET_l = np.stack(Pre_VRET_l[:len(SampleSet)])
    VUND_l = np.stack(Pre_VUND_l[:len(SampleSet)])
    print(VUND_l)
    Step = args.RollStep
    VB_c = np.sum(VB_l[:,:,:,:Step], axis = 3)

    VC_c = np.sum(VC_l[:,:,:,:Step], axis = 3)
    VOP_c = np.sum(VOP_l[:,:,:,:Step], axis = 3)
    Conditions = VB_l, VC_l, VOP_l,  VB_c, VC_c,VOP_c, VRET_l, VUND_l
    VB_l, VC_l, VOP_l,  VB_c, VC_c,VOP_c, VRET_l, VUND_l = Conditions

    pd.DataFrame(VB_l[:,:,:,0].sum(2)).to_csv(f'{args.FolderName}/{args.FolderName}_largeoutputs/VB_0.csv')
    #B_pd['0'] = Pre_VB_l
    #B_pd.to_csv(args.FolderName + '/'+'VB.csv')

    VB_log = []
    VC_log = []
    VOP_log = []
    VRET_log = []

    VB_log.append(VB_l)
    VC_log.append(VC_l)
    VOP_log.append(VOP_l)
    VRET_log.append(VRET_l)
    
    
    #### Iterating post
    it = 1
    it_max = args.RollCount

    

    while it <= it_max:
        if isinstance(args.Max_Increase, float):
            limit= [(l_i * (1+args.Max_Increase)) for l_i in limit]
            print('New Module Cap',limit)
        if it*Step >= args.ITC_Duration:
            args.ITC = 0.0
        Profits_h = []
        VB_l_h = []
        VC_l_h = []
        VOP_l_h = []
        VG_l_h = []

        Pre_Profits = []
        Pre_Max_Porfs = []
        Pre_VB_l = []
        Pre_VC_l = []
        Pre_VOP_l = []
        Pre_VRET_l = []
        Pre_VUND_l = []
        Pre_VG_l = []
        Pre_G_ind = []

        temp_DemN = []
        temp_Prof = []



        print('#########################', list(np.sum(Conditions[4],axis = (0,2))))
        n_helper = 0

        for s in SampleSet:
            print('####################### Facility ID', s, Fac_dict[s])
            temp_Profits = []
            temp_VB_l = []
            temp_VC_l = []
            temp_VOP_l = []
            temp_VG_l = []
            temp_VRET_l = []
            temp_VUND_l = []
            built = []
            exist = []

            for i in Fac_dict[s]:
                FR = mainFullOpt(NGCost = args.NGP, hrcount = args.HrCount, sites = Gen_File, 
                                 startyear = args.YearStart+(Step*it),  startHr = 0, years = args.YearCount,
                                 Facilities = Facs_total, FacN = i, Conditions = Conditions, Initial = False,
                                 Step = Step, n2 = n_helper, iterationN = it)

                try:
                    float(FR.vTotalProfits)
                except:
                    FR.vTotalProfits = -1000000000000
                temp_Profits.append(FR.vTotalProfits)
                #FR.vNG_pd.to_csv(args.FolderName+'/'+args.FolderName+'_largeoutputs/'+str(s)+'_'+str(i)+'_'+str(it)+'_vNG.csv')
                #FR.vRGen_pd.to_csv(args.FolderName+'/'+args.FolderName+'_largeoutputs/'+str(s)+'_'+str(i)+'_'+str(it)+'_vRGen.csv')


                #print('Und',FR.vUnd)
                print('B', FR.vB)
                #print('op', FR.vOp)
                #print('C', FR.vC)
                #print('Ret', FR.vRet)
                #FR.vNG.tofile(args.FolderName+'/'+args.FolderName+'_largeoutputs/'+str(s)+'_'+str(i)+'_'+str(it)+'_vNG')
                try:
                    print(FR.vTotalProfits)
                    if FR.vTotalProfits > 1:
                        built.append(1)
                    else:
                        built.append(0)


                    temp_VB_l.append(FR.vB[0])
                    temp_VC_l.append(FR.vC[0])
                    temp_VOP_l.append(FR.vOp[0])
                    temp_VG_l.append(FR.vGeneratorChoice[0])
                    
                    temp_VRET_l.append(FR.vRet[0])
                    temp_VUND_l.append(FR.vUnd[0])
                except:
                    print('EXCEPTION ISSUE')
                    FR.vTotalProfits = 0
                    built.append(0)
                    temp_VB_l.append(FR.vB[0])
                    temp_VC_l.append(FR.vC[0])
                    temp_VOP_l.append(FR.vOp[0])
                    temp_VG_l.append(FR.vGeneratorChoice[0])  
                    temp_VRET_l.append(FR.vRet[0])
                    temp_VUND_l.append(FR.vUnd[0])
                    

            if sum(built)>=1:
                MaxProf_Dem = temp_Profits.index(max(temp_Profits))
                Pre_Profits.append(temp_Profits[MaxProf_Dem])
                Pre_VB_l.append(temp_VB_l[MaxProf_Dem])
                Pre_VC_l.append(temp_VC_l[MaxProf_Dem])
                Pre_VOP_l.append(temp_VOP_l[MaxProf_Dem])
                Pre_VG_l.append(temp_VG_l[MaxProf_Dem])
                Pre_Max_Porfs.append(temp_Profits.index(max(temp_Profits)))
                Pre_G_ind.append(list(temp_VG_l[MaxProf_Dem]).index(1))
                Pre_VRET_l.append(temp_VRET_l[MaxProf_Dem])

                Pre_VUND_l.append(temp_VUND_l[MaxProf_Dem])
                temp_DemN.append(MaxProf_Dem)
                temp_Prof.append(temp_Profits[MaxProf_Dem])

                Fac_dict[s] = [Fac_dict[s][MaxProf_Dem]] 
            elif sum(exist)>=1:
                MaxProf_Dem = temp_Profits.index(max(temp_Profits))
                Pre_Profits.append(temp_Profits[MaxProf_Dem])
                Pre_VB_l.append(temp_VB_l[MaxProf_Dem])
                Pre_VC_l.append(temp_VC_l[MaxProf_Dem])
                Pre_VOP_l.append(temp_VOP_l[MaxProf_Dem])
                Pre_VG_l.append(temp_VG_l[MaxProf_Dem])
                Pre_Max_Porfs.append(temp_Profits.index(max(temp_Profits)))
                Pre_G_ind.append(list(temp_VG_l[MaxProf_Dem]).index(1))
                Pre_VRET_l.append(temp_VRET_l[MaxProf_Dem])
                Pre_VUND_l.append(temp_VUND_l[MaxProf_Dem])
                temp_DemN.append(5)
                temp_Prof.append(temp_Profits[MaxProf_Dem])

                Fac_dict[s] = [Fac_dict[s][MaxProf_Dem]] 
            else:
                Pre_Profits.append(-1000000000000000) #temp_Profits[0])
                Pre_VB_l.append(temp_VB_l[0])
                Pre_VC_l.append(temp_VC_l[0])
                Pre_VOP_l.append(temp_VOP_l[0])
                Pre_VG_l.append(temp_VG_l[0])
                Pre_Max_Porfs.append(0)
                Pre_G_ind.append(8)
                Pre_VRET_l.append(temp_VRET_l[0])
                Pre_VUND_l.append(temp_VUND_l[0])
                temp_DemN.append(4)
                temp_Prof.append(0)
            print(built)
            n_helper+=1
        Demand_pd[str(it)] = temp_DemN
        Profit_pd[str(it)] = temp_Prof
        #B_pd[str(it)] = Pre_VB_l

        Demand_pd.to_csv(args.FolderName + '/'+'demand_selection.csv')
        Profit_pd.to_csv(args.FolderName + '/'+'profit.csv')
        #B_pd.to_csv(args.FolderName + '/'+'VB.csv')

        ##########
        toBuilt_indexes =[]
        for gi in range(Gen_Count):
            mask_gi = np.array([gi == si for si in Pre_G_ind])
            sorted_profs = np.argsort(mask_gi*np.array(Pre_Profits))[::-1]
            built_n = sum(mask_gi)
            if built_n > int(limit[gi]):
                toBuild = int(limit[gi])
            else:
                toBuild = built_n
            toBuilt_indexes+=list(sorted_profs[:toBuild])

        toBuilt_indexes

        built = []
        gen_ind = []
        mod_count = []

        fc = 0
        for fb in Pre_VB_l[:len(SampleSet)]:
            mb = np.sum(fb, axis = (1,2))
            mi = np.sum(fb, axis = (0,2))
            built.append(mi)
            mod_count.append(mb.sum())
            mg = np.clip(mb, 0,1)
            if sum(mg)>0:
                gen_ind.append(list(mg).index(1))
            else:
                gen_ind.append(10)

        toBuilt_indexes_f =[]
        toBuilt_index_m = []

        for gi in range(Gen_Count):
            mask_gi = np.array([gi == gt for gt in gen_ind])

            built_g = sum(mask_gi)
            built_n = sum(mask_gi*np.array(mod_count))
            if built_n > int(limit[gi]):
                toBuild = int(limit[gi])
            else:
                toBuild = built_n
            print(built_n, toBuild)

            sorted_profs = np.argsort(mask_gi*np.array(Pre_Profits))[::-1]
            b_h = 0
            for sp in sorted_profs[:built_g]:

                # filter for the first year builds only
                # If no first year builds, then paste the 'if b_h...' section for no build

                # If first years, then create mask of first year build, multiply by v_3d_2.
                ## Sum up the first year and add that instead of toBuild-b_h in second clause
                years = Pre_VB_l[sp].shape[2]
                mods = Pre_VB_l[sp].shape[1]
                gens = Pre_VB_l[sp].shape[0]

                if (Pre_VB_l[sp][gi][:,0].sum() == 0) & (it != it_max):
                    print('No first year builds',Pre_VB_l[sp][gi])

                    Pre_VB_l[sp][gi]*=0
                    years = Pre_VB_l[sp].shape[2]
                    v2 = np.array(np.array([1]*VB_c[sp].sum()+[0]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h)-VB_c[sp].sum())))
                    v2 = v2[:mods]
                    v_3d_2 = np.array([[m]*years for m in v2])
                    Pre_VC_l[sp][gi]*=v_3d_2
                    Pre_VOP_l[sp][gi]*=v_3d_2
                    temp_und = Pre_VUND_l[sp]
                    Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                    Pre_VRET_l[sp]*=v_3d_2
                    print('no build_FY0', v2.sum())

                    b_h+=0
                
                else:
                    

                    FY0_build = Pre_VB_l[sp][gi][:,0].sum()
                    print('FY0',FY0_build)

                    v = np.array(np.array([0]*VB_c[sp].sum()+[1]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h)-VB_c[sp].sum())))
                    v2 = np.array(np.array([1]*VB_c[sp].sum()+[1]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h)-VB_c[sp].sum())))
                    
                    v_FY0 = np.array(np.array([0]*VB_c[sp].sum()+[1]*(FY0_build)+[0]*(len(built[sp])-(FY0_build)-VB_c[sp].sum())))
                    v2_FY0 = np.array(np.array([1]*VB_c[sp].sum()+[1]*(FY0_build)+[0]*(len(built[sp])-(FY0_build)-VB_c[sp].sum())))
                    
                    print('v',v)

                    v_mask = v[:mods]* v_FY0
                    v2_mask = v2[:mods]* v2_FY0

                    nB = v_mask.sum()


                    if b_h >= toBuild:
                        Pre_VB_l[sp][gi]*=0
                        years = Pre_VB_l[sp].shape[2]
                        vn = np.array(np.array([1]*VB_c[sp].sum()+[0]*(toBuild-b_h)+[0]*(len(built[sp])-(toBuild-b_h)-VB_c[sp].sum())))
                        v_3d_2 = np.array([[m]*years for m in vn])
                        Pre_VC_l[sp][gi]*=v_3d_2
                        Pre_VOP_l[sp][gi]*=v_3d_2
                        temp_und = Pre_VUND_l[sp]
                        Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                        Pre_VRET_l[sp]*=v_3d_2
                        print('no build', vn.sum())

                        b_h+=0

                    elif sum(built[sp])>toBuild-b_h:
                        print('VB_l sum', VB_c[sp].sum())
                        
                        
                        try:
                            v_3d = np.array([[m]*years for m in v_mask])
                            Pre_VB_l[sp][gi]*=v_3d
                        except:
                            v = v[:mods]
                            v_3d = np.array([[m]*years for m in v_mask])
                            Pre_VB_l[sp][gi]*=v_3d
                            print('EXCEPTION IN SIZES HIT, ln 1555')
                        print(v_3d)
                    
                        v2_mask = v2_mask[:mods]
                        v_3d_2 = np.array([[m]*years for m in v2_mask])
                        Pre_VC_l[sp][gi]*=v_3d_2
                        Pre_VOP_l[sp][gi]*=v_3d_2
                        temp_und = Pre_VUND_l[sp]
                        Pre_VUND_l[sp] = np.clip(((1-v_3d_2)+temp_und),0,1)
                    
                        Pre_VRET_l[sp]*=v_3d_2
                        b_h+=nB ###################################
                        print('partial build',nB,'to', b_h)
                    else:
                        b_h+=sum(built[sp])
                        print('full build', sum(built[sp]),'to', b_h)
            print('Built per gen',b_h)

        VB_l_h = np.stack(Pre_VB_l[:len(SampleSet)])

        VC_l_h = np.stack(Pre_VC_l[:len(SampleSet)])
        VOP_l_h = np.stack(Pre_VOP_l[:len(SampleSet)])
        VG_l_h = np.stack(Pre_VG_l[:len(SampleSet)])
        VRET_l_h = np.stack(Pre_VRET_l[:len(SampleSet)])
        VUND_l_h = np.stack(Pre_VUND_l[:len(SampleSet)])

        #VB_l_h = np.array(VB_l_h)
        #VC_l_h = np.array(VC_l_h)
        #VOP_l_h = np.array(VOP_l_h)
        #VG_l_h = np.array(VG_l_h)

        VB_c_h = np.sum(VB_l_h[:,:,:,:Step], axis = 3)+VB_c
        VC_c_h = (np.sum(VC_l_h[:,:,:,:Step], axis = 3)+VC_c) *VB_c_h   ###### 7/7/23 troubleshooting exports

        VOP_c_h = (np.sum(VOP_l_h[:,:,:,:Step], axis = 3)+VOP_c ) *VB_c_h

        #print(VC_c_h)
        #print(VC_l_h)
        Conditions = VB_l_h, VC_l_h, VOP_l_h,  VB_c_h, VC_c_h,VOP_c_h, VRET_l_h, VUND_l_h
        #print(Conditions)
        VB_log.append(VB_l_h)
        pd.DataFrame(VB_l_h[:,:,:,0].sum(2)).to_csv(f'{args.FolderName}/{args.FolderName}_largeoutputs/VB_{it}.csv')  # args.FolderName+'/'+args.FolderName+'_largeoutputs'
        VC_log.append(VC_l_h)
        VOP_log.append(VOP_l_h)
        VRET_log.append(VRET_l_h)
        it+=1

        VB_l = VB_l_h

        VB_c = VB_c_h
        VC_c = VC_c_h
        VOP_c = VOP_c_h
    Demand_pd.to_csv(args.FolderName + '/'+'demand_selection.csv')
    Profit_pd.to_csv(args.FolderName + '/'+'profit.csv')
    ##### End optimization, begin processing

    VC_multi = []
    VB_multi = []
    VOP_multi = []


    Step = args.RollStep
    it_max = args.RollCount
    it = 0
    for z in VC_log:
        tempf = []
        for f in z:
            tempg = []
            for g in f:
                tempm = []
                for m in g:
                    new = [0]*(it*Step)+list(m)+[0]*((it_max-it)*Step)
                    tempm.append(new)
                tempg.append(tempm)
            tempf.append(tempg)
        tempf = np.array(tempf)
        print(tempf.shape)
        VC_multi.append(tempf)
        it+=1
    VC_sum = []
    for z in range(len(VC_multi)):
        VC_sum.append(VC_multi[z])


    VC_sum

    VC_df = pd.DataFrame()
    VC_graph = pd.DataFrame()
    it_c = 0

    for it in VC_sum:
        iteration = it_c
        it_c +=1
        f_c = 0
        for f in it:
            facility = f_c
            f_c+=1
            g_c = 0
            for g in f:
                generator = g_c
                g_c+=1
                m_c = 0
                for m in g:
                    module = m_c
                    m_c +=1 
                T = [iteration,facility,generator]
                temp_pdf = pd.DataFrame()
                temp_pdf['Value'] = np.sum(g,axis = 0)
                temp_pdf['Iteration'] = iteration
                temp_pdf['Facility'] = facility
                temp_pdf['Generator'] = generator
                temp_pdf['Year'] = range(g.shape[1])
                VC_graph = pd.concat((VC_graph,temp_pdf))





    #### VB 
    it = 0
    for z in VB_log:
        tempf = []
        for f in z:
            tempg = []
            for g in f:
                tempm = []
                for m in g:
                    new = [0]*(it*Step)+list(m)+[0]*((it_max-it)*Step)
                    tempm.append(new)
                tempg.append(tempm)
            tempf.append(tempg)
        tempf = np.array(tempf)
        print(tempf.shape)
        VB_multi.append(tempf)
        it+=1
    VB_sum = []

    for z in range(len(VB_multi)):
        VB_sum.append(VB_multi[z])

    VB_df = pd.DataFrame()
    VB_graph = pd.DataFrame()
    it_c = 0

    for it in VB_sum:
        iteration = it_c
        it_c +=1
        f_c = 0
        for f in it:
            facility = f_c
            f_c+=1
            g_c = 0
            for g in f:
                generator = g_c
                g_c+=1
                m_c = 0
                for m in g:
                    module = m_c
                    m_c +=1 

                T = [iteration,facility,generator]
                temp_pdf = pd.DataFrame()
                temp_pdf['Value'] = np.sum(g,axis = 0)
                temp_pdf['Iteration'] = iteration
                temp_pdf['Facility'] = facility
                temp_pdf['Generator'] = generator
                temp_pdf['Year'] = range(g.shape[1])
                VB_graph = pd.concat((VB_graph,temp_pdf))


    ###### OP

    it = 0
    for z in VOP_log:
        tempf = []
        for f in z:
            tempg = []
            for g in f:
                tempm = []
                for m in g:
                    new = [0]*(it*Step)+list(m)+[0]*((it_max-it)*Step)
                    tempm.append(new)
                tempg.append(tempm)
            tempf.append(tempg)
        tempf = np.array(tempf)
        print(tempf.shape)
        VOP_multi.append(tempf)
        it+=1
    VOP_sum = []

    for z in range(len(VOP_multi)):
        VOP_sum.append(VOP_multi[z])

    VOP_df = pd.DataFrame()
    VOP_graph = pd.DataFrame()
    it_c = 0

    for it in VOP_sum:
        iteration = it_c
        it_c +=1
        f_c = 0
        for f in it:
            facility = f_c
            f_c+=1
            g_c = 0
            for g in f:
                generator = g_c
                g_c+=1
                m_c = 0
                for m in g:
                    module = m_c
                    m_c +=1 

                T = [iteration,facility,generator]
                temp_pdf = pd.DataFrame()
                temp_pdf['Value'] = np.sum(g,axis = 0)
                temp_pdf['Iteration'] = iteration
                temp_pdf['Facility'] = facility
                temp_pdf['Generator'] = generator
                temp_pdf['Year'] = range(g.shape[1])
                VOP_graph = pd.concat((VOP_graph,temp_pdf))
    VC_df
    VC_graph['Type'] = 'C'
    VB_graph['Type'] = 'B'
    VOP_graph['Type'] = 'Op'
    pd.concat((VB_graph,VC_graph,VOP_graph)).to_csv(args.FolderName + '/'+ args.NGP+'_'+str(args.YearCount)+'_binaries_logging.csv')


    VC_graph_sum = pd.DataFrame()
    it_c = 0

    for it in VC_sum:
        iteration = it_c
        it_c +=1
        g_c = 0
        it_gens = np.sum(it,axis= (0,2))
        for g in it_gens:
            generator = g_c
            g_c+=1
            temp_pdf = pd.DataFrame()
            temp_pdf['Value'] = g
            temp_pdf['Iteration'] = iteration
            temp_pdf['Generator'] = generator
            temp_pdf['Year'] = range(g.shape[0])
            VC_graph_sum = pd.concat((VC_graph_sum,temp_pdf))
    print(VC_graph_sum)


    #### VB 

    VB_graph_sum = pd.DataFrame()
    it_c = 0

    for it in VB_sum:
        iteration = it_c
        it_c +=1
        g_c = 0
        it_gens = np.sum(it,axis= (0,2))
        for g in it_gens:
            generator = g_c
            g_c+=1
            temp_pdf = pd.DataFrame()
            temp_pdf['Value'] = g
            temp_pdf['Iteration'] = iteration
            temp_pdf['Generator'] = generator
            temp_pdf['Year'] = range(g.shape[0])
            VB_graph_sum = pd.concat((VB_graph_sum,temp_pdf))

    VOP_graph_sum = pd.DataFrame()
    it_c = 0

    for it in VOP_sum:
        iteration = it_c
        it_c +=1
        g_c = 0
        it_gens = np.sum(it,axis= (0,2))
        for g in it_gens:
            generator = g_c
            g_c+=1
            temp_pdf = pd.DataFrame()
            temp_pdf['Value'] = g
            temp_pdf['Iteration'] = iteration
            temp_pdf['Generator'] = generator
            temp_pdf['Year'] = range(g.shape[0])
            VOP_graph_sum = pd.concat((VOP_graph_sum,temp_pdf))

    VC_graph_sum['Type'] = 'C'
    VB_graph_sum['Type'] = 'B'
    VOP_graph_sum['Type'] = 'Op'

    pd.concat((VB_graph_sum,VC_graph_sum,VOP_graph_sum)).to_csv(args.FolderName + '/'+ args.NGP+'_'+str(args.YearCount)+'_binaries_logging.csv')


    VB_graph['Type'] = 'B'
    VC_graph['Type'] = 'C'
    VOP_graph['Type'] = 'Op'

    pd.concat((VB_graph,VC_graph,VOP_graph)).to_csv(args.FolderName + '/'+ args.NGP+'_'+str(args.YearCount)+'_binaries_logging_ALL.csv')

    graphAll(VB_graph,VC_graph,VOP_graph, VB_graph_sum,VC_graph_sum,VOP_graph_sum)

Optimizations(Facs_total, Fac_dict, rand_n = args.Rand_n)