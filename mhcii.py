#To convert HTML file to pandas DataFrame
def MHCIIHTMLconverter(file):
    import pandas as pd
    with open (file, 'r') as Rawdata:
        MHCII_df = pd.DataFrame()
        x = 0
        for line in Rawdata:
            line = line.replace('\n', '')
            if line[:6] == 'allele':
                col = line.split('\t')
                MHCII_df = pd.DataFrame(columns=col)
            elif line[0:1] == 'H':
                mhcii_df_line = line.split('\t')
                if len(mhcii_df_line) != len(col):
                    print(f'ValueError: cannot set a row with mismatched columns in {file}')
                    print(f"Line: {mhcii_df_line}")
                    new_cells = len(col) - len(mhcii_df_line)
                    for cell in range(new_cells):
                        mhcii_df_line.append('[NaN]')
                MHCII_df.loc[x] = mhcii_df_line
                x += 1
    return(MHCII_df)
#Cutting Dataframe by HLA type of interest
def HLA_df_slice(df, HLAlist, HLAtype):
    delete_alleles = []
    for line in HLAlist:
        if f'HLA-{HLAtype}' not in line:
            delete_alleles.append(line)
    for index in range(len(delete_alleles)):
        df_remove = df.loc[(df["allele"] == delete_alleles[index])]
        final_df = df.drop(df_remove.index)
        df = final_df
    return final_df
#Get ID infos like specie and protein in the filename
def files_map_MHCIIBD(directory):
    import os
    MHCII_files = []
    specie = []
    protein = []
    for folders, subfolders, files in os.walk(directory):
        for file in files:
            path_file = (os.path.join(folders, file))
            path_file = str(path_file)
            if '.html' in file:
                MHCII_files.append(path_file)
                filename = file.replace('.html', '')   
                filename = filename.upper()
                filename = filename.replace('_PROTEINA', '')
                filename = filename.replace('_PROTEIN', '')
                filename = filename.replace('_PROT', '')
                infos = filename.split('_')
                if infos[0] > infos[1]:
                    specie.append(infos[0])
                    protein.append(infos[1])
                else:
                    specie.append(infos[1])
                    protein.append(infos[0])
    return (specie, protein, MHCII_files)
#To organize epitopes predicted by MHCII - Binding (IEDB)
# We use the results of the NN_align algorithm as default 
def MHCIIAntigenEpitopes(file, specie, protein, allele_class, ic50):
    import pandas as pd
    df = MHCIIHTMLconverter(file)
    #DataFrame Slice for NN_align Selection results
    NN_df = df[[ 'allele', 'seq_num', 'start', 'end', 'method', 'peptide', 'smm_align_ic50','nn_align_ic50', 'nn_align_rank', 'nn_align_adjusted_rank' ]]
    sp = [f'{specie}'] * (len(NN_df['peptide']))
    NN_df = NN_df.assign(specie = sp)#Add Species Column
    prot = [f'{protein}'] * (len(NN_df['peptide']))
    NN_df = NN_df.assign(protein = prot)#Add Proteins Column
    NN_df_remove = NN_df.loc[(NN_df['nn_align_ic50'] == '-')]
    new_NN_df = NN_df.drop(NN_df_remove.index)
    #Prediction Label for strong and weak bindings
    predNN = []
    for line in new_NN_df['nn_align_ic50']:
        line = float(line)
        if line < 50:
            predNN.append('SB')
        elif line < 500:
            predNN.append('WB')
        else:
            predNN.append('-')

    new_NN_df['Prediction_NN'] = predNN

    #To list HLA types in Dataframe
    HLAlist = list(new_NN_df['allele'].unique())
    
    #DataFrame Selection by HLA type
    if allele_class == 'DR':
        new_NN_df = HLA_df_slice(new_NN_df, HLAlist, allele_class)
    elif allele_class == 'DP':
        new_NN_df = HLA_df_slice(new_NN_df, HLAlist, allele_class)
    elif allele_class == 'DQ':
        new_NN_df = HLA_df_slice(new_NN_df, HLAlist, allele_class)
    else:
        print('Allele List Complete')

    #DataFrame Selection by IC50 Quartile(method NN_Align)
    if ic50 == 50:
        new_NN_df['nn_align_ic50'] = pd.to_numeric(new_NN_df['nn_align_ic50'])
        cut_NNSB_remove = new_NN_df.loc[(new_NN_df['nn_align_ic50'] > 50)]
        cut_NNSB = new_NN_df.drop(cut_NNSB_remove.index)
        cut_NNSB = cut_NNSB.reset_index(drop = True)
    elif ic50 <= 0:
        print("IC50 Value doesn't exist !")
    else:
        new_NN_df['nn_align_ic50'] = pd.to_numeric(new_NN_df['nn_align_ic50'])
        cut_NNSB_remove = new_NN_df.loc[(new_NN_df['nn_align_ic50'] > ic50)]
        cut_NNSB = new_NN_df.drop(cut_NNSB_remove.index)
        cut_NNSB = cut_NNSB.reset_index(drop = True)
    #Building the Final DataFrame with the data of interest
    df = cut_NNSB[['specie', 'protein', 'allele', 'start', 'end', 'peptide']]
    df = df.assign(method = 'MHC-II Binding')#Add Method Column
    df = df.assign(ID_Sequence = '-')#Add ID_Sequence Column
    df = df[['method', 'specie', 'protein', 'allele', 'ID_Sequence', 'start', 'end', 'peptide']]
    #standardizing column names with dataframes for other methods 
    df = df.rename(columns={
        'method': 'Method', 
        'specie': 'Specie', 
        'protein': 'Protein',
        'allele': 'Allele',
        'start': 'Initial Position',
        'end': 'Final Position',
        'peptide': 'Peptide Sequence'
        })
    
    return df