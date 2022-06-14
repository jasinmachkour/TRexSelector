## Changes in tknock version 0.0.0.9000 2022-06-13)

* Initial pre-release version contains five major functions for executing the T-Rex selector (and/or executing its major building blocks separately):

	- tknock();
	- random_experiments();
	- lm_dummy();
	- add_dummies(); and
	- add_dummies_GVS().

* It includes two functions for computing the false discovery proportion (FDP) and the true positive proportion (TPP):

	- FDP(); and
	- TPP().
	
* It includes three non-exported helper functions:

  - fdp_hat();
  - Phi_prime_fun(); and
  - select_var_fun().
	
* Vignette illustrates how the tknock package is used.
