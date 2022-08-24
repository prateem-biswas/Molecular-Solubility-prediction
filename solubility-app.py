# importing libraries

import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors

try :
	from rdkit.Chem import Draw
except :
	pass


# Calculate molecular descriptors
def AromaticProportion(m) :
	AromaticAtom = sum([m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())])
	HeavyAtom = Descriptors.HeavyAtomCount(m)
	AP = AromaticAtom/HeavyAtom
	return AP

def calculate(smiles) :
	moldata= []
	for elem in smiles :
		mol=Chem.MolFromSmiles(elem)
		moldata.append(mol)

	desc_MolLogP = []
	desc_MolWt = []
	desc_Rotatable_Bonds = []
	desc_Aromatic_Proportion = []    

	for mol in moldata :
		desc_MolLogP.append(Descriptors.MolLogP(mol))
		desc_MolWt.append(Descriptors.MolWt(mol))
		desc_Rotatable_Bonds.append(Descriptors.NumRotatableBonds(mol))
		desc_Aromatic_Proportion.append(AromaticProportion(mol))

	table = {"MolLogP" : desc_MolLogP , "MolWt" : desc_MolWt , "RotatableBonds": desc_Rotatable_Bonds, "AromaticProportion" : desc_Aromatic_Proportion }
	descriptors = pd.DataFrame(table)

	return descriptors


st.markdown("""
# Molecular Solubility Predictor App
""")

image = Image.open('solubility-logo.jpg')
st.image(image, use_column_width = True)

st.write("""
#### Use this app to get the predicted **Solubility(LogS)** values of molecules!
""")

steps ="Step 1 : Enter your SMILES code on the box on left <br> Step 2 : Click on the (Get Values) button <br> step 3 : Get the values "


st.sidebar.header('User Input Features')
# Smiles Input
SMILES_input = "c1ccncc1\nCC=C(C)O\nO=C(O)c1ccccc1"

with st.sidebar.form("Input"):
	SMILES = st.text_area("Enter SMILES code", SMILES_input)
	click = st.form_submit_button("Get Values")

if click : #going into the click function


	#SMILES = "C\n" + SMILES 		#Adds C as a dummy, first item
	SMILES = SMILES.split('\n')

	## Remove blank items from list to prevent error 
	for item in SMILES:
		if item == "" :
			SMILES.remove(item)

	#st.header('SMILES values :')
	#SMILES[1:] # Skips the dummy first item


	## Calculate molecular descriptors

	X = calculate(SMILES)
	#X[1:] # Skips the dummy first item

	#Loading the pickle file
	load_model = pickle.load(open('solubility_model.pkl', 'rb'))
	# Apply model to make predictions
	prediction = load_model.predict(X)
	#prediction_proba = load_model.predict_proba(X)

	st.header('Computed Molecular Descriptors')

	i = 0
	for item in SMILES :
		st.header(f"Molecule {i+1}")
	
		mol=Chem.MolFromSmiles(item)

		try :		#Using try except to bypass the problem Chem.Draw problem online 
			im = Draw.MolToImage(mol)
			st.markdown("""**Stucture :**""")
			st.image(im)
		except :
			pass

		st.markdown(f""" **SMILES** : {item} <br>
			**LogP** : {X.iloc[i,0]:.2f} <br>
			**Mol_Wt** : {X.iloc[i,1] :.2f} <br>
			**Rotatable Bonds** : {X.iloc[i,2]:.0f} <br> 
			**Aromatic Proportion** : {X.iloc[i,3]:.2f}<br>
			<h3> **LogS** : {prediction[i]:.2f} </h3>
			""",True)

		i += 1

	
	steps =""


st.sidebar.write(""" ### Examples of some SMILES codes """)
st.sidebar.markdown("""
	Benzene : c1ccccc1 <br>
	Ethanol : CCO <br>
	Methane : C <br>
	Octanol : CCCCCCCCO <br>
	Methyl amine : CN <br>
	Anisole : COc1ccccc1 <br>
	Aspirin : O=C(C)Oc1ccccc1C(=O)O <br>
	""",True)

st.markdown(f"""{steps}""", True)
	
st.write("""
Data obtained from the John S. Delaney. [ESOL: Estimating Aqueous Solubility Directly from Molecular Structure](https://pubs.acs.org/doi/10.1021/ci034243x). ***J. Chem. Inf. Comput. Sci.*** 2004, 44, 3, 1000-1005.
***
	""")
