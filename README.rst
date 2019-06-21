PeleProduction - README
=======================

`PeleProduction` is a collection of run folders for various Pele codes and processing. It includes git submodules for the dependent codes.


Getting started
---------------

To clone the repository and its dependent submodules, do ::

    git clone --recursive git@github.com:AMReX-Combustion/PeleProduction.git

To build one of the `PeleLM` cases, for example ::

    cd PeleProduction/PeleLMruns/ClosedChamber
    make -j 12

Note that the example will set the variables pointing to `PeleLM`, `IAMR`, `PelePhysics` and `AMReX` (to point to the respective folders in the Submodules area here) only if these variables have not been set in the user's environment.

Pulling with submodules
-----------------------

You can update the repository with fetch/pull like you would normally do. To pull everything including the submodules, use the `--recurse-submodules` and the `--remote` parameter in the git pull command ::

    # pull all changes in the repo including changes in the submodules
    git pull --recurse-submodules
    
    # pull all changes for the submodules
    git submodule update --remote

Updating the submodule commits being tracked
--------------------------------------------
The relevant state for the submodules are defined by the `PeleProduction` repository. If you commit here, the state of the submodules are also defined by this commit. The `git submodule update` command sets the Git repository of the submodule to that particular commit. The submodule repository tracks its own content which is nested into the main repository. The main repository refers to a commit of the nested submodule repository.

Use the `git submodule update` command to set the submodules to the commit specified by the main repository. This means that if you pull in new changes into the submodules, you need to create a new commit in your main repository in order to track the updates of the nested submodules ::

     # update submodule in the master branch
     # skip this if you use --recurse-submodules
     # and have the master branch checked out
     cd [submodule directory]
     git checkout master
     git pull

     # commit the change in main repo
     # to use the latest commit in master of the submodule
     cd ..
     git add [submodule directory]
     git commit -m "move submodule to latest commit in master"

     # share your changes
     git push

Another developer can get the update by pulling in the changes and running the submodules update command ::

     # another developer wants to get the changes
     git pull

     # this updates the submodule to the latest
     # commit in master as set in the last example
     git submodule update


Acknowledgment
--------------
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
