#Get epitopes predictes from NetCTL-1.2 Server
def netctlAntigenEpitopes(file):
    import pandas as pd
    import html
    NetCTL_df = pd.DataFrame()
    #Create a DataFrame with columns of the same name as those in the NetCTL result table
    columns = 'Residue_number ID Protein_identifier pep Peptide_sequence aff Predicted_MHC_binding_affinity aff_rescale Rescale_binding_affinity cle C_terminal_cleavage_affinity tap TAP_transport_efficiency COMB Prediction_score Identified_MHC_ligands'
    columns = columns.split(' ')
    NetCTL_df = pd.DataFrame(columns=columns)
    with open(file,'r') as RawData:
        text = RawData.readlines()
    x = 0
    #Filling in the DataFrame with netctl results (html)
    for line in text:
        line = html.unescape(line) #Converting html named and numeric character references (e.g. &lt) in Unicode characters
        if not line.isspace():
            line = line.strip()
            if (line[0]) in '0123456789':
                line = line.replace('   ',' '); line = line.replace('  ',' ');
                netctl_df_line = line.split(' ')
                if line [-1] == 'E':
                    NetCTL_df.loc[x] = netctl_df_line
                else:
                    netctl_df_line.append('-')
                    NetCTL_df.loc[x] = netctl_df_line
                x += 1
    #to Select the predicted epitopes by NetCTL method
    NetCTL_epitopes = NetCTL_df[NetCTL_df['Identified_MHC_ligands'] == '<-E']
    NetCTL_epitopes = NetCTL_epitopes.reset_index(drop = True)
    #Set data of interest in the DataFrame
    NetCTL_epitopes_slice = NetCTL_epitopes[['Protein_identifier', 'Residue_number', 'Peptide_sequence']]
    #Insert columns for method, id, specie, and protein into the DataFrame
    NetCTL_epitopes_slice = NetCTL_epitopes_slice.assign(Method='NetCTL')
    NetCTL_epitopes_slice = NetCTL_epitopes_slice.assign(ID_Sequence='-') #NetCTL does not save the ID information, but the column has been included for standardization
    sp = []; prot = []
    for line in NetCTL_epitopes_slice['Protein_identifier']:
        line = line.split('_')
        sp.append(line[1])
        prot.append(line[0])
    NetCTL_epitopes_slice['Specie'] = sp
    NetCTL_epitopes_slice['Protein'] = prot
    #From the length, creating a new column for the final position of the amino acid residue in the epitope         
    FinalPosition_aa = []
    for line in range(len(NetCTL_epitopes_slice['Peptide_sequence'])):
        final_Pos = int(NetCTL_epitopes_slice['Residue_number'][line]) + (len(NetCTL_epitopes_slice['Peptide_sequence'][line]) - 1)
        FinalPosition_aa.append(final_Pos)
    NetCTL_epitopes_slice = NetCTL_epitopes_slice.rename(columns={'Residue_number': 'Initial Position', 'Peptide_sequence': 'Peptide Sequence'})
    NetCTL_epitopes_slice['Final Position'] = FinalPosition_aa
    #Organizing information in the final DataFrame
    results_df = NetCTL_epitopes_slice[['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']]
    return (results_df)