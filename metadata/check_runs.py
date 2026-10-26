import pandas as pd

df = pd.read_csv("metadata/SraRunTable.csv")   

print("Unique samples (GSM):", df["GEO_Accession (exp)"].nunique())   
print(df.groupby("GEO_Accession (exp)").size())                       
print(df["Assay Type"].value_counts())                                
print(df["treatment"].value_counts())                                 

samples = (df.groupby("GEO_Accession (exp)")
             .agg(treatment=("treatment", "first"),
                  runs=("Run", lambda x: ",".join(x)))
             .reset_index())
print(samples)
samples.to_csv("metadata/sample_metadata_draft.csv", index=False)
