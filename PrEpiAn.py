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
        predictions.append(Bp_df)
        
    #Analysis of epitopes predicted by IMED (Predicting Antigenic Peptides)
    #Parsing the txt file
    if (args.p) != '':
        print("Collecting PAP epitope data ...\n")
        if (args.pmin) == 0 and (args.pmax) == 0:
            PAP_df = papImed.PAPepitopes(args.p, 0, 0)
        else:
            PAP_df = papImed.PAPepitopes(args.p, args.pmin, args.pmax)
        predictions.append(PAP_df)
            
    #Analysis of epitopes predicted by NetCTL 1.2
    #Analyzing the file in html format
    if (args.n) != '':
        print("Collecting NetCTL epitope data ...\n")
        NetCTL_df = netctl.netctlAntigenEpitopes(args.n)
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
        for indice in range(len(MHCII_files)):
            df = mhcii.MHCIIAntigenEpitopes(MHCII_files[indice], specie[indice], protein[indice], args.mhla, args.mic)
            MHCII_df = pd.concat([MHCII_df, df])
        predictions.append(MHCII_df)

    #Analysis of epitopes predicted by others methods
    #Parsing the fasta file
    if (args.x) != '':
        filename = os.path.basename(rf"{args.x}")
        print("Analyzing data from other predictors ...\n")
        if (args.xmin) == 0 and (args.xmax) == 0:
            epitopes_df = othersPred.fasta_epitopes(args.x, 0, 0)
        else:
            epitopes_df = othersPred.fasta_epitopes(args.x, args.xmin, args.xmax)
        predictions.append(epitopes_df)
    
    #Build a single dataframe for the result of all predictions 
    df_predictions = pd.DataFrame()
    columns = ['Method', 'Specie', 'Protein', 'ID_Sequence', 'Initial Position', 'Final Position', 'Peptide Sequence']
    df_predictions = pd.DataFrame(columns=columns)
    for prediction in predictions:
        df_predictions = pd.concat([df_predictions, prediction], join="inner")
    
    #result files (by method) in .xlsx format
    if (args.e).lower() == 'y':
        #
        Bp_df.to_excel(rf'{args.d}/Bepipred_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        #
        PAP_df.to_excel(rf'{args.d}/PAP_IMED_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        #
        NetCTL_df.to_excel(rf'{args.d}/NetCTL_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)
        #
        MHCII_df.to_excel(rf"{args.d}/MHCII-Binding_Epitopes.xlsx", sheet_name = 'Prediction Results', startcol=0, index=False)
        #
        epitopes_df.to_excel(rf'{args.d}/{filename}_Epitopes.xlsx', sheet_name = 'Prediction Results', startcol=0, index=False)

    return(df_predictions)