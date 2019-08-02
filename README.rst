PeleProduction - README
=======================

``PeleProduction`` is a collection of run folders for various Pele codes and processing. It includes git submodules for the dependent codes (such as ``PeleLM``, ``PeleC``, ``AMReX``, etc).  This organizational strategy is a work-in-progress attempt to manage the interactions between the various dependent repositories -- ie to keep them all compatible with each other.  So far, I have found that the system works reasonably well as long as you don't get too fancy and try to use the "recursive" clone thing that people try to do.  Instead, this note explains how I use it with some success.


Getting started
---------------

This repo consists of the working environment of several of our collaborators.  They are organized into git branches.  The submodules feature of git is used to track the specific commits of each of the dependent repositories that are compatible with the current state of each of those branches. So, the ``jbb`` branch might require one version of ``PeleLM`` while the ``UNSW`` branch needs another, and each depends on different branches of ``AMReX``.  In order to deal with this annoyance, there is a file in the root of this repository called ``.gitmodules`` that lists the local and remote locations of each dependent repo, and the specific commit for each of these repos is tucked away in ``PeleProduction``'s git history.

To clone all the sutff you need for a specific branch ::

    # Do a "shallow" clone of this repo
    git clone git@github.com:AMReX-Combustion/PeleProduction.git
    cd PeleProduction

Then change to the desired branch, e.g. ``UNSW`` ::

    git checkout UNSW

The first time you do this, you will need to tell git that there are submodules (a necessary evil of submodules).  Do ::

    git submodule init

Git will look at the ``.gitmodules`` file in this branch and use that to initialize stuff (note that a different branch may have a different version of that file, and therefore a different list of dependent repos). Finally, get the correct commits of the sub-repos set up for this branch ::

    git submodule update

Now, you are ready to to build one of the cases associated with this branch, for example ::

    cd PeleProduction/PeleLMruns/ClosedChamber
    make -j 12

As a side note, this "make" will assume that the usual environment variables pointing to ``PeleLM``, ``IAMR``, ``PelePhysics`` and ``AMReX`` have NOT been set, and will instead point to the local sub-repo folders.


Pulling/updating with submodules
--------------------------------

With subrepositories, there is an extra step to updating your local working files ::

    git pull
    git submodule update

The first command will grab commits for source in ``PeleProduction``. The second will do pulls in each the sub-repos, and set the commit to the one associated with the current ``PeleProduction`` branch.  If you make changes to ``PeleProduction`` files, ``git status`` will show the usual info.  If the commit of any of the sub-repos was changed in the ``PeleProduction`` commit,``git status`` will show that your sub-repo is out of date; the second command there will bring any of those in line. 

Updating the submodule commits being tracked
--------------------------------------------

The relevant state for the submodules are defined by the ``PeleProduction`` repository. If you commit here, the state of the submodules are also defined by this commit. The ``git submodule update`` command sets the Git repository of the submodule to that particular commit. The submodule repository tracks its own content separately. If you introduce new changes into any of the submodules, you will need to create a new commit in ``PeleProduction`` to update the pointer in git that selects the specific commit you want ::

     # e.g., update the master branch of a submodule
     cd [submodule directory]
     git checkout master
     git pull

     # commit the change in PeleProduction
     cd ..
     git add [submodule directory]
     git commit -m "move submodule to latest commit in master"

     # share your changes
     git push

Another developer can get the update by pulling in the changes and running the submodules update command ::

     # another developer wants to get the changes to PeleProduction, including new sub-repo commit pointers
     git pull

     # this updates all the submodules to the latest
     # commit in master as set in the last example
     git submodule update


More annoying aspects of submodules...sorry
--------------------------------------------

While the strategy outlined here is a way to manage the synchronization of multiple repos, it has some odd features.  One is that the ``git submodule update`` command restores the sub-repos to specific commits...not branches. If you cd into that sub-repo, you will be in a "detached head" state at that commit, even if that corresponds exactly to a branch, such as ``development``.  Moreover, if you edit in the local sub-repo and try to commit, the commit will be on top of a detached head, not ``development`` and almost certainly not what you intended.

Here's an example.  Say I've set the ``marc`` branch of ``PeleProduction`` to use the ``development`` branch of ``AMReX`` as a sub-repo.  If I move into my local clone of amrex, it will be in a detached head at the HEAD of ``development``.  If I commit a fix to ``AMReX`` there, it will not be applied to ``development`` but rather to the detached head.  I should rather ``git checkout development`` in my AMReX repo before committing any fixes there.  Then I can commit and push to AMReX as usual, and others will see it as a new commit to ``development``.  But I also have to move back to the ``marc`` branch of ``PeleProduction`` and checkin the pointer to this new commit of AMReX.  After doing so, when someone else checks out my branch of ``PeleProduction``, and does ::

     cd PeleProduction
     git checkout marc
     git pull
     git submodule update

my new commit to AMReX will be available at the remote site, and will be applied and set as the current version in the user's local submodules.

A second weird behavior was hinted at above, related to "recursive" clones.  I was tempted to include such a thing in these instructions, but here'sthe problem... ``PeleC`` was recently changed to itself contain submodules for ``PelePhysics`` and ``AMReX``.  If the ``master`` branch of ``PeleProduction`` contained a subrepo for ``PeleC``, a recursive clone would get multiple copies of things -- what a mess!  So, I opted for the manual approach, even if it is a little more verbose.

All this is a bit of torture, but in my experience it is still better than the alternative of manually keeping lists of branches and commits that are compatible.  Let me know if you stumble on something even better.

-M


Acknowledgment
--------------
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
