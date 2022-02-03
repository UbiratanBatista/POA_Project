#Get and organize epitopes from fasta file
def fasta_epitopes(fasta, lenght_min, lenght_max):
    import pandas as pd
    #Create a new data frame (table for the results)
    results_df = pd.DataFrame()
    colunas = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    results_df = pd.DataFrame(columns=colunas)
    #Open prediction file in fasta format
    with open (fasta, 'r') as Input:
        index = 0
        for line in Input:
            if not line.isspace():
                line = line.upper(); line = line.strip()
                if line[0] == '>': #Get sequence information from the epitope header
                    ID_data = line.replace('>', '')
                    ID_data = ID_data.split('_')
                    sp = ID_data[1]
                    prot = ID_data[0]
                    method = ID_data[2]
                    if ID_data[3] == 'NP':
                        idSeq = f"{ID_data[3]}_{ID_data[4]}"
                    else:
                        idSeq = ID_data[3]
                    InitialPos = ID_data[-2]
                    FinalPos = ID_data[-1]
                else: #Get sequence of the epitopes
                    Epitope = line
                    line = [
                        method,
                        sp,
                        prot,
                        idSeq,
                        InitialPos,
                        FinalPos,
                        Epitope
                        ]
                    
                    results_df.loc[index] = line
                    index += 1
    #To set the minimum and maximum length of epitopes
    if lenght_min == 0 and lenght_max == 0 :
        results_df = results_df
    elif lenght_min == 0:
        results_df = results_df[results_df['Peptide_sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    elif lenght_max == 0:
        results_df = results_df[results_df['Peptide_sequence'].str.len() >= lenght_min]
        results_df = results_df.reset_index(drop = True)
    else:
        results_df = results_df[results_df['Peptide_sequence'].str.len() >= lenght_min]
        results_df = results_df[results_df['Peptide_sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    return results_df