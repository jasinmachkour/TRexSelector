## Changes in trex version 0.0.1 2022-07-22)

* First release:

  - Polish pre-release version.

## Changes in trex version 0.0.0.9000 2022-06-13)

* Initial pre-release version contains five major functions for executing the T-Rex selector (and/or executing its major building blocks separately):

	- trex();
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
	
* Vignette illustrates how the trex package is used.
