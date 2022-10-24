# Surface_Charging

This is how you do a full surface charging calculation for studying electro-catalytic systems. 
Please feel free to add any thing else to this repository related to this method.

THEORY:
Regular CHE approach : 
Under the CHE scheme, the chemical potential of an electron−proton pair can be approximated in reference to the free energy of half a dihydrogen molecule, and it also allows the effect of the applied potential and pH to be easily included in the model.

Surface charging: 
Under the CHE approach, the electrochemical potential is assumed to affect only the chemical potential of the exchanged electrons. Thus, electronic energy is independent with respect to the potential and is taken from neutral systems. However, the free energy for the sites and the intermediates in the "neutral" situation is not the same as the U = 0.0 V situation.


STEPS ( I take these steps to reduce the computational time for the calculation) : 
1. First perform a full relaxation calculation of your system, you can use ase_geometry_opt.py for it. 
2. Symmetrize your system using the Symmetrize_simple.py ( NOTE: you might have to edit the symm.py based on your system but this should work most of the times)
    -- You will have to change the part you want to keep fixed in the slab. This would usually be the middle part of the system.
3. Run a relax + solvation calculation on the new system that you generated by the Symmetrize_simple.py this can be done with ase_symmetric_opt.py but change the paramaeters depending on the system.
    -- Use the "" file to run a relax(without LSOL) and then it will also create a folder for LSOL that will have the calcualtion with solvation
4. Once the solvation relation is done we then move on to do surface charging calculatation. Follow the followin steps: 
    -- create a folder SC inside your LSOL folder. 
    -- Copy the following files from the LSOL folder to SC folder : OUTCAR and CONTCAR 
    -- change the file name of CONTCAR to POSCAR inside the SC folder.
    -- put the ase_surface_charging.py in the SC folder and change its name to ase-bfgs.py ( ase python file that will do the vasp calcualtion, so change any input parameters accoding to the system) and the submission script as job.sh ( The names are important to be kept as the same)
5. Run the gen-scpot.py file inside the SC folder. This will create all the folders with different NELECT valules and also submit the calculation in each folder. You dont need to do anything else if your job.sh, ase-bfgs.py etc files are good .
6. These calculations do take considerable time so it can be taht you need to resubmit certain jobs in certain NELECT folder, you will have to do this step manually.
7. Once the converage criteria is reached in each folder, run the plot_sc.py, this will give you the paprabolic equation and also give out certain plots. To visualize if your plots are fitting to the points , please visualize the g-pot.png file before using the papabola to use in further analysis.


REFERENCES: 
