Codes for "Live-cell single RNA imaging reveals bursts of translational frameshifting"
=======

Kenneth Lyon, Luis U. Aguilera, Tatsuya Morisaki, Brian Munsky, and Timothy J. Stasevich.

Department of Biochemistry and Molecular Biology, Institute of Genome Architecture and Function, Colorado State University, Fort Collins, CO, 80523, USA
Department of Chemical and Biological Engineering and School of Biomedical Engineering, Colorado State University, Fort Collins, CO, 80523, USA
World Research Hub Initiative, Institute of Innovative Research, Tokyo Institute of Technology, Yokohama, Kanagawa, 226-8503, Japan

For questions about the codes, please contact:  Luis.aguilera@colostate.edu and brian.munsky@colostate.edu

---

This repository contains the codes necessary to reproduce figures from the above manuscript. All codes are implemented in Matlab 2018.

## Code organization

The codes are organized on three main cathegories.

  * **runner__FrameShifting.m** performs the stochastic simulation and reproduces the percentage of spots in shifted and non-shifted states, number of ribosomes, intensities and harringtonine assays.
  * **runner__Optimization.m** performs the parameter optimization using the pattern search algorithm. 
  * **runner__Parameter_Uncertanty.m** performs a random search to evaluate parameter uncertanty.

---

## Experimental data.

The experimental data used for the harringtonine assays
* **0and-1F.xls** Harringtonine data for the (1F and 0F) FSS MF tag.
* **0F.xls** (0F) Harringtonine data for the FSS MF tag.
* **FS_smHA_0F_only.xls** Harringtonine data for the (0F) HA MF tag, intensity for the FLAG proble (not used in the final fit).
* **FS_smHA_SunTag.xls** Harringtonine data for the (-1F) HA MF tag, intensity for the SunTag proble (not used in the final fit).
* **FS_smHA.xls** Harringtonine data for the (0F) HA MF tag, selecting only the frameshifting spots
* **nFS_smHA.xls** Harringtonine data for the (0F) HA MF tag , selecting only the non-frameshifting spots

---

## Gene Sequences.
* **Frame_0.txt** (0F) FSS MF tag
* **Frame_1.txt** (-1F) FSS MF tag
* **2X_FS_0F.txt** (0F) 2xFSS MF tag
* **2X_FS_1F.txt** (-1F) 2xFSS MF tag
* **0F_HA.txt** (0F) HA MF tag
* **1F_HA.txt** (-1F) HA MF tag

---  

## Code implementation.

The codes performing the optimization and uncertanty (runner__Optimization.m and runner__Parameter_Uncertanty.m), rely upon expensive computations that were performed on the W.M. Keck Compute Cluster. 
