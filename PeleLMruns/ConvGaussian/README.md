Convected Gaussian - README
================================

Description
-----------

The convection of a scalar gaussian perturbation is a classical test case for computational fluid dynamics solvers. It allows to directly evaluate the numerical characteristics of the solver by comparison with analytical solution. It consists in convecting a gaussian perturbation of a scalar (temperature or species mass fraction) across a periodic 2D box. 

Current implementation
----------------------

Running the convergence tests
-----------------------------

./multiRuns.py --test_name ConvTemp2SDC_x --input_file inputs.2d-convT_posx

./pprocConvOrder.py fcompare.gnu.ex --test_name ConvTemp2SDC_x

...produces Convergence_ConvTemp2SDC_x.png

