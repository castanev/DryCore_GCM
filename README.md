This model is based on the Held and Suarez (1994) benchmark and then configured by Wu & Reichler, 2018. Slight modifications were made to this model and some new files were added for purposes of automatic post-processing of the model outputs and generation of new experiments. From line 4 to 27 you can find the detailed description of the model that Wu & Reichler, 2018 provide to the public at (https://github.com/ZhengWinnieWu/WR_simpleGCM). From line 28 you will find the description of the files and the step-by-step to run the model successfully. 


# WR_simpleGCM
#### Improving GFDL's Idealized General Circulation Model

We start our experiments using an analytically determined equilibrium temperature (Teq) and Newtonian relaxation time scale (Tau) profile after Jucker et al. (2014). We then gradually optimize Teq by introducing zonal asymmetries into it using the iterative procedure following Chang (2006). The Teq we provide here (teq_ite31_0.9.nc) corresponds to iteration 31 of DRAG experiment D3 with the surface drag being 0.9 1/day. The file contains 12 monthly Teq fields, each varying in the three spatial dimensions. We also employ an actual orography with land and ocean. Given these modifications, the GCM produces temperatures and diabatic heating that are very similar to that of the reanalysis. A detailed description of the improved model and the outcomes of the various drag experiments are described in Wu and Reichler (2018a, 2018b). 

This code is based on **JFV-strat** (https://github.com/mjucker/JFV-strat). The following program and modules were modified:<br />
•	atmos_model<br />
•	hs_forcing_mod<br />
•	atmosphere_mod<br />
•	time_manager_mod	

The new code reads in twelve monthly Teq and Tau profiles during initialization and then uses linear interpolation to calculate daily fields. The interpolation is done only once during initialization to reduce the actual run time. Additional modification concerns the inclusion of a more flexible surface drag. The Rayleigh damping (surface drag) in the original code is defined in module “hs_forcing_mod”, which is a fixed number. In our modification, we introduce the Rayleigh damping as a two-dimensional matrix (longitude-latitude), and allow the code to read in the Rayleigh damping data files. Certain new variables are read in through the name-list: <br />
•	equilibrium_option = 'from_file', <br />
•	drag_file_name = 'name of drag data file.txt', <br />
•	temp_file = 'name of Teq data file.nc', <br />
•	tau_file = 'name of Tau data file.nc'. <br />

The main references are:<br />
Chang, E. K. M., 2006: An idealized nonlinear model of the Northern Hemisphere winter storm tracks. J. Atmos. Sci., 63, 1818–1839.<br />
Held, I. I. M., and M. M. J. Suarez, 1994: A proposal for the intercomparison of the dynamical cores of atmospheric general circulation models. Bull. Amer. Meteor. Soc., 75, 1825–1830.<br />
Jucker, M., S. Fueglistaler, and G. K. Vallis, 2014: Stratospheric sudden warmings in an idealized GCM. J. Geophys. Res. Atmos., 119, 11 054–11 064, doi:10.1002/2014JD022170.<br />
Wu, Z. and T. Reichler (2018a): Towards a More Earth-like Circulation in Idealized Models, *J. Adv. Model. Earth Sys.*, **30**(24), 10101-10116. <br />
Wu, Z. and T. Reichler (2018b): Surface Control of the Frequency of Stratospheric Sudden Warming Events, J. Climate (in preparation).


# DESCRIPTION OF FOLDERS AND FILES 
data: 
- teq_ite31_0.9.nc: Is the T_equilibrium used as input in the Wu & Reichler, 2018 model and described before.
- Mean_state_T_NCEP_asymmetrical_4D.nc: is the temporal mean zonal mean of the global temperature from NCEP. Used as T_equilibrium for one of the experiments. Structure: (time, lev, lat, lon)
- Mean_state_T_NCEP_symmetrical_4D.nc: is the temporal mean of the global temperature from NCEP. Used as T_equilibrium for one of the experiments. Structure: (time, lev, lat, lon)
- mean_state_NCEP.py: python code used to obtaining the two input files above for experiments.
- Figures: contains figures showing the distribution of the temperature in the T_equilibrium files.


fms_hs_winnie
- compile: all files needed for compilation. You must read the README_COMPILE.txt file and follow the step by step for compilation.
- Post-processing: automatically used at the end of each model run to decompress and organize the outputs in netcdf files. 
- src: is the source code. You must review **.../src/atmos_param/hs_forcing/hs_forcing.f90** to change parameters if needed.  


- vc_dry_job.sbatch: Is the file to run the model and post-processing is also done there. PATHS HAVE TO BE SPECIFIED. 
- vc_run_dry.sh: Is where vc_dry_job.sbatch is submited to the cluster using SLURM. PARAMETERS HAVE TO BE SPECIFIED. 


# STEP BY STEP TO RUN THE MODEL 
1. Change paths in **vc_dry_job.sbatch**, the name of the experiment and introduce your email in the SBATCH. In the variable **temp_file** you must define the file for the T_equilibrium.
2. Specify in **vc_run_dry.sh** the name of the experiment. You can change the parameters for SLURM.
3. Change paths in **input.nml** for the topog_file_name
4. Change the name of the experiment and the specified path at the beginning of the code **fms_hs_winnie/Post-processing/post_processing.py**.
5. Compile the model. To do that, go to **fms_hs_winnie** and follow the instructions in README_COMPILE. If you run the model in Brown, yo can skip the steps 1, 2 and 3. Otherwise, yo have to configurate the compilers for a different cluster. 
6. Run the model with **vc_run_dry.sh**. To to that, just run in the terminal the command [ ./vc_run_dry.sh ] in the folder where the file is located. 
7. Review the outputs located in the path specified into the sbatch file (**$archive/$name**). In the **post_processed** file will be the final outputs, each file corresponding to 365 days (in daily outputs) or 12 months (in monthly outputs). 
8. If you want to obtain a file for the temperature for all the simulated period, EDIT THE PATHS AND EXPERIMENT NAME of the code **.../fms_hs_winnie/Post_processing/Post_processing.py** and wun it using the sbatch code located in the same folder. 