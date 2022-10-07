**************************************************************************************************************************
*****************************************Programs for the Study*********************************************************

*******************Part I Simulation Setting I to Setting IV*************************
1. Generate Datasets for setting I to IV. R
    - Used to generate 200 datasets for each setting

2. Macro program setting I to IV. sas
   - Estimate the model parameters under settings I to IV

 
3. Summary the fitting results from setting I to IV.R
   -Summary the estimation of model parameters under settings I to IV
   - Resutls shown in Table S1-S4



*******************Part II Self Defined Function*************************
1. Lambdainv.R
    -Function used to find vector T such that Lambda(T)=s for a vector input s

2. myprod.R
    -Fucntion used to product of two step function

3. mysimRec
   -Function using Cinlar's inversion Method to generate Non-homogeneous Poisson process

4. datagen.R
   - Function used to generate the dataset with known covariates X and it's involved in computing the true NDE/NIE

5. datagen_M.R
   -Function used to generate recurrent events M0 and M1

6. datagen_X.R
   -Function used to generate datasets with simulated X and Z

7. datagenQ13_M.R
   -Function used to generate recurrent events M_Q1 and M_Q3
   -Q1 and Q3 are two comparing categories for NDE/NIE

8. mymed.R
   -Function used to compute the NDE/NIE for simulation settings

9. mymed2.R
   -Function used to compute the true NDE/NIE

10. mymed_Q13.R
   -Function used to compute the NDE/NIE with Z containing two categories Q1 and Q3

11. S.R
    - Used to compute the survival functions

12. S00.R
     - Compute the survival functions S00

13. S01.R
     - Compute the survival functions S01

14. S10.R
     - Compute the survival functions S10
 
15. S11.R
     - Compute the survival functions S11
 
   

*******************Part III Estimate NDE/NIE Under Simulation Setting I to IV*************************

1. Folder "NDE NIE from Simulation I"
    *3-1. SimulationI_parallel.R
              - Function used to estimate NDE/NIE under simulation setting I
    *3-2. SimulationI_parallel_boot.R
             -Function of the bootstrap for NDE/NIE under simulation setting I
     *3-3. SimulationI_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting I


2. Folder "NDE NIE from Simulation II"
    *4-1. SimulationII_parallel.R
              - Function used to estimate NDE/NIE under simulation setting II
    *4-2. SimulationII_parallel_boot.R
             -Function of the bootstrap for NDE/NIE under simulation setting II
    *4-3. SimulationII_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting II


3. Folder "NDE NIE from Simulation III"
    *5-1. SimulationIII_parallel.R
              - Function used to estimate NDE/NIE under simulation setting III
    *5-2. SimulationIII_parallel_boot.R
             -Function of the bootstrap for NDE/NIE under simulation setting III
    *5-3. SimulationIII_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting III


4. Folder "NDE NIE from Simulation IV"
    *6-1. SimulationIV_parallel.R
              - Function used to estimate NDE/NIE under simulation setting IV
    *6-2. SimulationIV_parallel_boot.R
             -Function of the bootstrap for NDE/NIE under simulation setting IV
    *6-3. SimulationIV_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting IV




*******************Part IV Real Data Analysis*************************

1. CPCRA study analysis.sas
    -Estimates the parameters under simulation setting I and II for CPCRA study

2. Folder "CPCRA NDE NIE Plot"

   * NDE_NIE_Trt.R 
      -Estimate NDE/NIE for treatment

   * NDE_NIE_CI_Trt.R
     -Boostrap of NDE/NIE for treatment

   * Summary_Trt.R
    -Summary the results and plot the estimates of NDE/NIE/TE with bootstraped 95% CI for treatment

  * NDE_NIE_CD4.R 
      -Estimate NDE/NIE for CD4

   * NDE_NIE_CI_CD4.R
     -Boostrap of NDE/NIE for CD4

   * Summary_CD4.R
    -Summary the results and plot the estimates of NDE/NIE/TE with bootstraped 95% CI for CD4



 