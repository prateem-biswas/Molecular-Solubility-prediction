# Molecular Solubility prediction App
## End to end project
![solubility-logo](https://user-images.githubusercontent.com/105977171/186254360-53e53dab-9887-4a32-95ff-5a4d7fcf29d1.jpg)


The solubility of various organic molecules is an important property in Drug discovery, design and development .The drug likeness property of a molecule is depends on
how easily it can dissolve and reach its target in the body, and solubility plays a big part in that!
Synthesizing compounds for potential drugs can be both time consuming and expensive, only to find out that the molecule is not compatible in terms of solubility! 
The following application tries to solve this issue building utilising machine learning to predict the solubility value of the molecule, along with some other parameters 
like log of partion cofficient between water and octanol (LogP), molecular weight, number of rotatable bonds and aromatic proportion. It needs the structure of the 
molecule in one dimensional form (SMILES notation) and returns the above values!


Click on the button below to use the app <br>
[![click](https://user-images.githubusercontent.com/105977171/186254747-e24a4395-747f-4d8b-8664-1cca7c354480.png)](https://prateem-biswas-molecular-solubility-predi-solubility-app-0gplt6.streamlitapp.com/)

## Project Description

The first step require using a dataset that contains the solubility values of various small organic molecules. This data was obtained from the paper published by John S Delaney : [ESOL: Estimating Aqueous Solubility Directly from Molecular Structure](https://pubs.acs.org/doi/10.1021/ci034243x). The data consists of more than 1100 rows containing Compound ID, their measured and predicted solubility values, and their SMILES notation.
![Capture](https://user-images.githubusercontent.com/105977171/186433157-06bf6970-f0e5-4224-ab95-c25311f6835e.PNG)

The **simplified molecular-input line-entry system (SMILES)** is a specification in the form of a line notation for describing the structure of chemical species using short ASCII strings. These SMILES values are used to calculate the Molecular descriptors of each molecule. A molecular descriptor is a structural or physicochemical property of a molecule . These descriptors are calculated using the **RDKit** library in Python. In his experiment, Delaney used 4 molecular descriptors : 
* **LogP (Log of octanol-water partition coefficient) :** Octanol-water partition coefficient is ratio of the concentration of a solute in a water-saturated octanolic phase to its concentration in an octanol-saturated aqueous phase. taking the Log of that value gives us LogP
* **MWt (Molecular weight):** The molecular weight of the compound
* **RB (Number of Rotatable Bonds) :** Number of  single non-ring bonds, attached to a non-terminal, non-hydrogen atom.
* **AP (Aromatic proportion) :** The ratio of number of aromatic atoms in the molecule to that of the total number of non hydrogen atoms
