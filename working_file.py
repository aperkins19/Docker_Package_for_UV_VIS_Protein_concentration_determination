import pandas as pd

from Bio.SeqUtils.ProtParam import ProteinAnalysis




protein_sequences_df = pd.read_csv("proteinsequences.csv")

protein_sequences_df  = protein_sequences_df.iloc[:,1:]

reduced=[]
oxidised = []
for i, row in protein_sequences_df.iterrows():

    reduced.append(ProteinAnalysis(row["Sequence"]).molar_extinction_coefficient()[0])
    oxidised.append(ProteinAnalysis(row["Sequence"]).molar_extinction_coefficient()[1])

protein_sequences_df["reduced_extinction_coefficient"] = reduced
protein_sequences_df["oxidised_extinction_coefficient"] = oxidised

print(protein_sequences_df)

protein_sequences_df.to_excel("protein_sequences_extinction_coeff.xlsx")