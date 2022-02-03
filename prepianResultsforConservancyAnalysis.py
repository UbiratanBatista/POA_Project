#To create epitope files (fasta) for Epitope Conservancy Analysis
#Files separated by organism species 
def Epitope_EptConsAnalysis(dataframe, path):
    #To list species in Dataframe
    Specieslist = list(dataframe["Specie"].unique())
    for sp in Specieslist:
        #Select predicted epitopes (by species) 
        df_Specie = dataframe.loc[(dataframe['Specie'] == sp)]
        df_Specie = df_Specie.reset_index(drop = True)
        #Create fasta files (by species) with epitopes
        with open (f'{path}/Conservancy Analysis/{sp}_epitopes.fasta', 'w') as epitopes:
            for linha in range(len(df_Specie["Peptide Sequence"])):
                id_epitope = f">{df_Specie['Specie'][linha]}_{df_Specie['Protein'][linha]}_{df_Specie['Method'][linha]}_{df_Specie['Initial Position'][linha]}_{df_Specie['Final Position'][linha]}"
                epitopes.write(f"{id_epitope}\n")
                epitopes.write(f"{df_Specie['Peptide Sequence'][linha]}\n")
                    
#Build the fasta files for next step (Epitope Conservancy Analysis)
def filesforEptConsAnalysis(FileofProteins, EpitopeDataframe, directory):
    #Creating a folder to receive the fasta files
    import os
    Analysis_directory = False
    for folders, subfolders, files in os.walk(directory):
        if folders == directory:
            if "Conservancy Analysis" in subfolders:
                Analysis_directory = True
        else:
            continue
    if Analysis_directory is not True:
        dir = f'{directory}/Conservancy Analysis'       
        os.mkdir(dir)
    #Creating the fasta files
    Epitope_EptConsAnalysis(EpitopeDataframe, directory)
