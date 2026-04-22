# Introduction

This is a python project dedicated to understanding Lefschetz fibrations, particularly through their monodromy representations. In particular the project is dedicated to exploring how the automorphism group of the (smooth) Lefschetz fibration maps to the smooth mapping class group of the base surface. At the moment it does not contain many features, but the author hopes to build upon these as they become relevant to her research. Any suggestions are welcome :)

## AI Usage

Recently, I have been experimenting with using AI (Claude Code) to optimize and rapidly iterate on the code contained here. Note that the original python project was written entirely without Claude. In particular, the computations for the 3 papers

1. "How large is the braid monodromy of low-genus Lefschetz fibrations"
2. "Torsion in the braid monodromy of elliptic fibrations"
3. "Universal braids for elliptic fibrations: from character varieties to Coxeter's factor groups of braid groups"

are contained in the PaperComputations folder and have not been touched by AI tools.

## For readers of my papers

If you have been directed here from my papers, self-contained versions of the code for those papers can be found in the PaperComputations subfolder. Most of these programs require sage to run, and can be run by inputting
```
sage [filename]
```
in the command line. The lib.sage file is an auxiliary file containing all the common code between the computations. There is one file which is run using gap, and can be run via
```
gap [filename]
```
in the command line. This file contains generators for Br(q_3), Br(q_4), and Br(q_5\[2\]) that are computed using the other files + the solomon package itself (for Br(q_5\[2\\)).

## Structure in general 

Currently the file [LefschetzDisk.ipynb](LefschetzDisk.ipynb) handles all of the calculations for genus 1 Lefschetz fibrations. It includes the ability to generate orbit graphs, compute generating sets, and more. These are generally included in the Setup Cells. The remainder of the sections are scratchpads which I use for my personal research, and contain no guarantee of being maintained.

An auxiliary program contained in [elliptic-subgroups.gap](elliptic-subgroups.gap) takes care of determining the presentations of certain finite index subgroups of interest.

These two files are not very well-documented/structured for general purpose use. I am currently trying to remove the dependencies on Sage--which is incompatible with the curver library--in order to create a package which works in higher genus. In the process I hope to document the code more thoroughly, and make it easier to use in general.

# Installation

TODO

# Dependencies

This package depends on
- curver
- flipper
- networkx
- distribute
- cysignals
- setuptools
We are particularly grateful to Mark Bell, the creator of curver/flipper, for his effort in creating computational tools for mapping class groups.

# License / Citing

This work is licensed under the [MIT License](LICENSE.md). However, if you find the package useful in your research, please consider citing it, after the package is officially published a suggested citation will be available.
