#Verification of the viability of epitopes to be analyzed by conservancy analysis 
def checkIntegrity(dataframe, column_epitopes, column_ID):
    indice = 0
    for line in dataframe[column_epitopes]:
        line = str(line)
        line = line.lower()
        if 'x' in line:
            raise Exception(f"ERROR: There is an epitope in sequence {dataframe.iloc[indice][column_ID]} containing an invalid character: 'x'")
        else:
            indice += 1

def runningPrEpiAn(args):
    import pandas as pd
    import mhcii
    import netctl
    import papImed
    import bepipred
    import othersPred
    import os

    predictions = []
    #Analysis of epitopes predicted by Bepipred
    #Parsing the Json file
    if (args.b) != '':
        print("Collecting Bepipred-2.0 epitope data ...\n")
        if (args.bmin) == 0 and (args.bmax) == 0:
            bp_data = bepipred.BpJsonAnalysis(args.b)
            Bp_df = bepipred.bpAntigenEpitopes(bp_data, 0, 0)
        else:
            bp_data = bepipred.BpJsonAnalysis(args.b)
            Bp_df = bepipred.bpAntigenEpitopes(bp_data, args.bmin, args.bmax)
        checkIntegrity(Bp_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(Bp_df)
        
    #Analysis of epitopes predicted by IMED (Predicting Antigenic Peptides)
    #Parsing the txt file
    if (args.p) != '':
        print("Collecting PAP epitope data ...\n")
        if (args.pmin) == 0 and (args.pmax) == 0:
            PAP_df = papImed.PAPepitopes(args.p, 0, 0)
        else:
            PAP_df = papImed.PAPepitopes(args.p, args.pmin, args.pmax)
        checkIntegrity(PAP_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(PAP_df)
            
    #Analysis of epitopes predicted by NetCTL 1.2
    #Analyzing the file in html format
    if (args.n) != '':
        print("Collecting NetCTL epitope data ...\n")
        NetCTL_df = netctl.netctlAntigenEpitopes(args.n)
        checkIntegrity(NetCTL_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(NetCTL_df)
        
    #Analysis of epitopes predicted by MHCII - Binding (IEDB)
    #Data collected from html files in a folder
    if (args.m) != '':
        print("Collecting MHCII epitope data ...\n")
        directory = rf"{args.m}"
        specie, protein, MHCII_files = mhcii.files_map_MHCIIBD(directory)
        
        MHCII_df = pd.DataFrame()
        columns = ['Method', 'Specie', 'Protein', 'Allele', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
        MHCII_df = pd.DataFrame(columns=columns)

        MHCII_df_list = []
        for file in MHCII_files:
            df = mhcii.MHCIIHTMLconverter(file)
            MHCII_df_list.append(df)
        for indice in range(len(MHCII_df_list)):
            df = mhcii.MHCIIAntigenEpitopes(MHCII_df_list[indice], specie[indice], protein[indice], args.mhla, args.mic)
            MHCII_df = pd.concat([MHCII_df, df])
        checkIntegrity(MHCII_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(MHCII_df)

    #Analysis of epitopes predicted by others methods
    #Parsing the fasta file
    if (args.x) != '':
        filename = str(os.path.basename(rf"{args.x}"))
        filename = filename.split('.')
        print("Analyzing data from other predictors ...\n")
        if (args.xmin) == 0 and (args.xmax) == 0:
            epitopes_df = othersPred.fasta_epitopes(args.x, 0, 0)
        else:
            epitopes_df = othersPred.fasta_epitopes(args.x, args.xmin, args.xmax)
        checkIntegrity(epitopes_df, 'Peptide Sequence', 'ID_Sequence')
        predictions.append(epitopes_df)
    
    #Build a single dataframe for the result of all predictions 
    df_predictions = pd.DataFrame()
    columns = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    df_predictions = pd.DataFrame(columns=columns)
    for prediction in predictions:
        df_predictions = pd.concat([df_predictions, prediction], join="inner")
    
    #result files (by method) in .xlsx format
    if (args.e).lower() == 'y':
        if args.b != '':
            Bp_df.to_excel(rf'{args.d}/Bepipred_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.p != '':
            PAP_df.to_excel(rf'{args.d}/PAP_IMED_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.n != '':
            NetCTL_df.to_excel(rf'{args.d}/NetCTL_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.m != '':
            MHCII_df.to_excel(rf"{args.d}/MHCII-Binding_Epitopes.xlsx", sheet_name = 'Prediction Results', startcol=0, index=False)
        if args.x != '':
            epitopes_df.to_excel(rf'{args.d}/{filename[0]}_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)

    return(df_predictions)