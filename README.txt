Modified model used in Gernon et al., 2024 - Nature Geoscience
Original model from Alcott et al., 2019 - Science
Coded by Lewis Alcott and Benjamin JW Mills // email lewis.alcott@bristol.ac.uk or b.mills@leeds.ac.uk
The only modifications to this code from the original are the routine to test multiple phosphorus inputs events, and a new constant for testing different rates of vertical mixing.

Contents:
Alcott_etal_2019_Science_front.m // Model frontend file // run this code to solve model and plot results // sets parameters and starting values, runs solver
Alcott_etal_2019_Science.m // Model equations file // do not run this code directly // contains flux and reservoir calculations

Reproduction of paper plots:
See line 87 in frontend file to modify number of model runs per forcing event. Currently set to 1 for demonstration but set to 1000 for runs shown in paper. 
*a full ensemble takes several hours on a single core and generates on the order of 1GB of data*

Notes:
Requires MATLAB. Tested in R2023b in Windows 10.
