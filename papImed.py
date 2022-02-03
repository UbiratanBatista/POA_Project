#Get Epitopes predictes from IMED Server
def PAPepitopes(file, lenght_min, lenght_max):
    import pandas as pd
    #Building a new dataframe with predicting antigenic peptides server (IMED) output information
    with open (file, 'r') as txt:
        data = txt.readlines()
        IMED_df = pd.DataFrame()
        col = ['ID', 'n', 'Start Position', 'Sequence', 'End Position']
        IMED_df = pd.DataFrame(columns=col)
        x = 0
        for line in data:
            if not line.isspace():
                line = line.strip()
                if line[0] == '>':
                    id = line.replace('>', '')
                    id = id.upper()
                elif line[0] == 'n':
                    continue
                elif line in '\t\n ':
                    continue
                else:
                    IMED_df_line = f"{id}\t{line}"
                    IMED_df_line = IMED_df_line.replace('\n', '')
                    IMED_df_line = IMED_df_line.split('\t')
                    IMED_df.loc[x] = IMED_df_line
            x += 1
        IMED_df = IMED_df.reset_index(drop = True)
    #Organizing the results in a new dataframe
    results_df = pd.DataFrame()
    col = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    results_df = pd.DataFrame(columns=col)
    index = 0
    #Select informations as specie, protein, id sequence, etc.
    for line in range(len(IMED_df['Sequence'])):
        imed_id = str(IMED_df['ID'][line]).split('_')
        sp = imed_id[1]
        protein = imed_id[0]
        if len(imed_id) > 2:
            if imed_id[2] == "NP":
                idSeq = imed_id[2] + "_" + imed_id[3]
            else:
                idSeq = imed_id[2]
        else:
            idSeq = '-'
        #Building the final dataframe with the results
        results_line = ['PAP/IMED',
                        sp,
                        protein,
                        idSeq,
                        IMED_df['Start Position'][line],
                        IMED_df['End Position'][line],
                        IMED_df['Sequence'][line]]      
        results_df.loc[index] = results_line
        index += 1
    #To set the minimum and maximum length of epitopes 
    if lenght_min == 0 and lenght_max == 0 :
        pass
    elif lenght_min == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    elif lenght_max == 0:
        results_df = results_df[results_df['Peptide Sequence'].str.len() >= lenght_min]
        results_df = results_df.reset_index(drop = True)
    else:
        results_df = results_df[results_df['Peptide Sequence'].str.len() >= lenght_min]
        results_df = results_df[results_df['Peptide Sequence'].str.len() <= lenght_max]
        results_df = results_df.reset_index(drop = True)
    return results_df