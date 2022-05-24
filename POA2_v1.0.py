#writing the fasta file of the results 
def getfastafile(path,portion_epitopes, dataframe):
        with open (f"{path}.fasta", "w") as results_fasta:
            for indice_line in range(len(dataframe["Epitope sequence"])):
                idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
                epitope = dataframe.iloc[indice_line]["Epitope sequence"]
                column = ''
                #Definition of which topology region will be referenced in the fasta file 
                if portion_epitopes == 0:#all regions 
                    results_fasta.write(f'>{idt_epitope}')
                    results_fasta.write(f'\n{epitope}\n')
                elif portion_epitopes == 1: #outer regions
                    column = "Portion_Outside"
                elif portion_epitopes == 2: #transmembrane regions 
                    column = "Portion_TM"
                else: #inner regions
                    column = "Portion_Inside"
                if column != '':
                    if dataframe.iloc[indice_line][column] == 1: #To check if the cell contains a specific text
                        results_fasta.write(f'>{idt_epitope}')
                        results_fasta.write(f'\n{epitope}\n')
                 
def main():
    import pandas as pd
    import argparse
    import conservancyAnalysis
    import PepiTMHMM
    import sys

    parser = argparse.ArgumentParser(add_help=False, description = 'POA2 - Conservancy Analysis and TMHMM')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-g", help= "-g: greater than or equal to (>=) in Sequence identity threshold from your analysis", type=str)
    mutually_exclusive_group.add_argument("-l", help= "-l: less-than sign (<) in Sequence identity threshold from your analysis", type=str)
    required.add_argument("-h", "--help", action ="help", default=argparse.SUPPRESS, help= "Print this help message.")
    required.add_argument("-d", help= "Directory for Results of the Epitope Conservancy Analysis in CSV format", type=str, default='', required = "True")
    required.add_argument("-t", help= "Sequence identity threshold", type=int, default='', required = "True")
    optional.add_argument("-r", help= "Directory for analysis results", type=str, default='')
    optional.add_argument("-imin", help= "Minimum identity (%%) in Conservancy Analysis", type=int, default='60')
    optional.add_argument("-imax", help= "Maximum identity (%%) in Conservancy Analysis", type=int, default='100')
    optional.add_argument("-m", help= "Percent of protein sequence matches at identity", type=int, default='60')
    required.add_argument("-f", help= " Protein sequence(s) fasta for prediction of transmembrane helices (TMHMM)", required = "True")
    optional.add_argument("-rf", help= '''Analysis results in fasta format: [0]all epitopes 
                                    [1]epitopes in Outside portion(TMHMM)
                                    [2]epitopes in Transmembrane portion(TMHMM)
                                    [3]epitopes in Inside portion(TMHMM)''', type=int, default=None)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    
    print("\n" + 'POA2 - Analysis Conservancy' + "\n")
    if (args.t == '') or (args.d == ''):
        raise Exception(f'ERROR: The following arguments are required: -t, -d')
    
    if args.rf != None:
        if (args.rf < 0) or (args.rf > 3):
            raise Exception(f'ERROR: Argument -rf is invalid')

    if args.r != '':
        path = f"{args.r}/POA2_analysis"
    else:
        path = f"POA2_analysis"

    if args.g:
        type_symbol = '>='
    else:
        type_symbol = '<'

    #organize the results of conservancy analysis 
    ConservancyAnalysis_DF = conservancyAnalysis.EpitConservAnalysis(args.t, type_symbol, args.m, args.imax, args.imin, args.d)   
    #applying the TMHMM prediction on the epitopes
    POA2_df = PepiTMHMM.tmhmmAnalysis(args, ConservancyAnalysis_DF)
    #assembling the spreadsheet with the final result 
    POA2_df.to_excel(f"{path}_{args.t}.xlsx", sheet_name = f'{type_symbol}{args.t}', startcol=0, index=False)

    if (args.rf) != None:
        getfastafile(path, args.rf, POA2_df)

    print("Analysis of Results Completed!\n")

if __name__ == '__main__':
    main()