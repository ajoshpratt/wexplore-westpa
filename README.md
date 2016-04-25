# wexplore-westpa
A prototype of the WExplore method for [WESTPA](http://chong.chem.pitt.edu/WESTPA). This is a work-in-progress, 
but it nominally does what it's suppose to.

Based on:

**WExplore: Hierarchical Exploration of High-Dimensional Spaces Using the Weighted Ensemble Algorithm**<br \>
Alex Dickson and Charles L. Brooks, III<br \>
The Journal of Physical Chemistry B 2014 118 (13), 3532-3542<br \>

http://dx.doi.org/10.1021/jp411479c


This is a fork of the original repo, designed to work with steady state simulations.

INSTALLATION:

After cloning, checkout the following branches of WESTPA and WEST_TOOLS:

WESTPA: new_hooks

WEST_TOOLS: recursive_update

Or, after loading anaconda python:

`git clone https://github.com/westpa/westpa.git`

`git checkout new_hooks`

`./setup.sh`

`sed -i "s/checkout_remote/\#checkout_remote/g" setup.sh`

`cd lib/west_tools`

`git checkout recursive_update`

`cd ../../`

`./setup.sh`

In the CH4/CH4 example, symlink or copy in the main wexplore directory (from the root).

`./init.sh`

`./run.sh`
