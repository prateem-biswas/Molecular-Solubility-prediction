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
The project consitsts of two parts. Part one is aimed at recreating the work of John S Delany by building a simple linear regression model used for predicting the solubility of the molecules, and obtaining a formula that represents the relation between solubility and other features. The second part consists of analysing the dataset, comparing various algorithms and choosing one for building a better model, and using the model to deploy a solubility predicting app. 

## Part 1 :

The first step require using a dataset that contains the solubility values of various small organic molecules. This data was obtained from the paper published by John S Delaney : [ESOL: Estimating Aqueous Solubility Directly from Molecular Structure](https://pubs.acs.org/doi/10.1021/ci034243x). The data consists of more than 1100 rows containing Compound ID, their measured and predicted solubility values, and their SMILES notation.
![Capture](https://user-images.githubusercontent.com/105977171/186433157-06bf6970-f0e5-4224-ab95-c25311f6835e.PNG)

The **simplified molecular-input line-entry system (SMILES)** is a specification in the form of a line notation for describing the structure of chemical species using short ASCII strings. These SMILES values are used to calculate the Molecular descriptors of each molecule. A molecular descriptor is a structural or physicochemical property of a molecule . These descriptors are calculated using the **RDKit** library in Python. In his experiment, Delaney used 4 molecular descriptors : 
* **LogP (Log of octanol-water partition coefficient) :** Octanol-water partition coefficient is ratio of the concentration of a solute in a water-saturated octanolic phase to its concentration in an octanol-saturated aqueous phase. taking the Log of that value gives us LogP
* **MWt (Molecular weight):** The molecular weight of the compound
* **RB (Number of Rotatable Bonds) :** Number of  single non-ring bonds, attached to a non-terminal, non-hydrogen atom.
* **AP (Aromatic proportion) :** The ratio of number of aromatic atoms in the molecule to that of the total number of non hydrogen atoms

![descriptos](https://user-images.githubusercontent.com/105977171/186515195-329ccec0-cb1a-449e-b09f-12ec58b8b43f.PNG) <br>
The first 3 desriptors can be calculated by a built in method in rdkit, while aromatic proportion is calculated manually using functions !

In the first part of the project the aim is to recreate the work of Delany and obtain a fomula with intercept and coefficent values for various descriptors. So, a simple linear regression model was constructed! The four descriptors make up the feature set and the measured solubility is the target variable. Using 80 percent of the data for training gave us a model with Mean Square Error (MSE) of 1.01 and R<sup>2</sup> score of 0.77 . Using the entire dataset for training marginally improves the performance. The final formula obtained through the intercepts and coefficients is as folllows :

Equation using train dataset : <br>
LogS = **0.25** **- 0.73** LogP **- 0.0066** MW **+ 0.0050** RB **- 0.50** AP

Equation using the entire dataset :<br>
LogS = **0.26** **- 0.74** LogP **- 0.0066** MW **+ 0.0032** RB **- 0.42** AP

These formula compares to others as follows : 

Equation from Delany : <br>
LogS = **0.16 -  0.63** cLogP  **- 0.0062** MW **+ 0.066** RB **- 0.74** AP

The reproduction by Pat Walters  provided the following: <br>
LogS = **0.26 - 0.74** LogP **- 0.0066** MW + **0.0034** RB **- 0.42** AP

So, by using the entire dataset, we have built a model that comes close to replicating the work done by Pat Walters

![scatter](https://user-images.githubusercontent.com/105977171/186516508-59147d59-975c-477a-9d9e-b3e81d7cd7ea.png)

## Part 2 :

The r2 score for the in the previous model left a lot of scope for improvement . For this the dataset containing the features and target values is loaded and split into training and testing data. To understand the various features better, the training data is visualized via pair plots :

![pairplot](https://user-images.githubusercontent.com/105977171/186646497-95087088-0648-4d35-9d08-15cb7c36597b.png)

The various plots reveal the following :<br>
* Both the LogP and LogS values are normally distributed. LogP has a mean and standard deviation of **2.458** and **1.89** respectively, while for LogS the values are **-3.07** and **2.10** <br>
![hist](https://user-images.githubusercontent.com/105977171/186651670-ffd36145-9922-4445-8371-557accea4e3d.png)
* Both LogP and molecular weight are inversely proportional with the solubility values , while they are themselves linearly correlated <br>

![subplot](https://user-images.githubusercontent.com/105977171/186654585-4e8662da-30a8-402f-aa0a-00cf73200523.png)

Other two features do not show much correlation with the target variable

For building a better model, various algorithms were tested and their R<sup>2</sup> scores were used as the evaluation matrix. It yeilded the following results:<br>

![model_comp](https://user-images.githubusercontent.com/105977171/186655512-d1dba845-e3e5-48e7-8a0d-5e86471832bb.png)

As observed, Random Forest model yeilded the best result with R<sup>2</sup> score of 0.88 . So it was chosen as the model for building the application.

[View the app](https://prateem-biswas-molecular-solubility-predi-solubility-app-0gplt6.streamlitapp.com/)

[Notebook : Part 1](https://github.com/prateem-biswas/Molecular-Solubility-prediction/blob/main/Predicting_solubity_of_molecules.ipynb)

[Notebook : Part 2](https://github.com/prateem-biswas/Molecular-Solubility-prediction/blob/main/Analysis_and_Model_building.ipynb)

