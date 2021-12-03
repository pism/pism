Flowline ice shelf model (to test for driving stress scheme)
=================


Run a test to check for the effect of the driving stress scheme for a 50m ice thickness perturbations 
at the calving front in a flow line ice shelf setup and compare to analytical solution (van Der Veen). 

As in the flow line case no buttressing force exists, the SSA velocity should be unaffected by such perturbations. 
As surface gradients are low towards the calving front margin, the one-sided difference scheme with calving front 
boundary conditionproduced solutions comparable to the analytical solution, but it cannot reproduce the expected 
(absent) perturbation effects. 

     $ bash run_test.sh

provides a plot of anomalies in driving stress and SSA velcoity along the flow line ice shelf.
