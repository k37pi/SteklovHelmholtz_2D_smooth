1. Change inputs in inputs.m

2. steklov_HL.m computes Steklov-Laplace/Helmholtz spectrum and efns. 

3. helm_bvp.m computes Dirichlet, Neumann BVPs for given boundary fn for Laplace and Helmholtz equations. 

4. steklov_convergence.m computes self convergence and numerical orders

5. bayes_script.m does the shape optimization. Change the respective input end points depending on the curve from inputs.m. Shape opt can be done for up to 7 different scalings methods, see sh_bayes_fn.m. Add other scalings if required.

6. optimizers_bayes_plots.m shows the plots for shape opt experiments. Change the curve name in inputs.m and use the csvs from opt/files/.     

Coming up: Dirichlet/Neumann/Robin-Laplace and transmission EVPs in this setup.
Note: Can also compute Dirichlet-Laplace efn but need to know the appropriate wave number from D-L code. 
