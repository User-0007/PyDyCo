# PyDYCO

**PyDYCO** is a thermodynamic solver primarily used for thermoelectricity, with the acronym standing for **Dynamique Couplée**. This solver leverages **NgSpice** and **Python** to calculate the coupling between energy current and electric current in 2D materials.

## Features

- **Coupling Calculation**: PyDYCO computes the interaction between energy and electric currents in 2D materials.
- **Temperature and Electrochemical Potential**: The solver helps find the local temperature and electrochemical potential in the material.
- **Current Distribution**: It calculates the different currents in the material, providing valuable insights into thermoelectric properties.

## Technologies Used

- **PySpice** Library and **NgSpice** for circuit simulation.
- **Python** for Network creation and computation.
- **ParaView** for data visualization and post-processing.
  
## Installation

To install PySpice, please refer to the official installation guide on [Fabrice Salvaire's website](https://pyspice.fabrice-salvaire.fr/).

## Usage
**PyDYCO** runs thermodynamic networks, and in this context, it's applied to thermoelectricity. The user defines the dimensions of the network, the size and properties of each element (thermoelectric dipole), the boundary conditions, and the operating mode.
<p align="center">
<img src="Pydyco_Network2.png" alt="7x11 element network simulated by PyDYCO with defined boundary conditions in terms of temperature and electrochemical potential." width="300" height="250">
</p>

( Usage instructions coming soon...)
