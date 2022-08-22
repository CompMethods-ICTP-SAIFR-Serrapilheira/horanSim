# MoranSim: simulation of a well-mixed Moran process with constant selection
**Created by Lucas Alves de Melo Pontes**

# Objective

The objective of this repository is to execute simulations of Moran's model without spatial restrictions.

The implementation of Moran's model is based on:
* the Ricardo Martínez-Garcia's classes on "Stochastical Mathematical Modeling" during the Serrapilheira-ICTP QBio course;
* the description presented in the chapters of the book "NOVAK, M. A. Evolutionary dynamics: Exploring the equations of life. Belknap Press of Harvard University Press. 378 p., Cambridge, Massachusetts, and London, England, 2006".

# Structure of repository
```
endProject_compMethods/
├── data/
|     ├── raw/
|     └── processed/
├── docs/
├── graph/
├── R/
├── R_modules/
└── README.md
```
The order of execution for the scripts is given by the script filenames inside the "R/" directory.

"R_modules/" contains the modules required for the scripts in the "R/" directory.

"data/" contains the tables obtained through the simulation.

"docs/" contains the report file and its supporting files.

"graph/" contains the generated graphs.
