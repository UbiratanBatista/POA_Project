def getfastafile(path,portion_epitopes, dataframe):
        with open (f"{path}.fasta", "w") as results_fasta:
            for indice_line in range(len(dataframe["Epitope sequence"])):
                idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
                epitope = dataframe.iloc[indice_line]["Epitope sequence"]
                column = ''
                if portion_epitopes == 0:
                    results_fasta.write(f'>{idt_epitope}')
                    results_fasta.write(f'\n{epitope}\n')
                elif portion_epitopes == 1:
                    column = "Portion_Outside"
                elif portion_epitopes == 2:
                    column = "Portion_TM"
                else:
                    column = "Portion_Inside"
                if column != '':
                    if dataframe.iloc[indice_line][column] == 1:
                        results_fasta.write(f'>{idt_epitope}')
                        results_fasta.write(f'\n{epitope}\n')
                 
def main():
    import pandas as pd
    import argparse
    import conservancyAnalysis
    import PepiTMHMM

    parser = argparse.ArgumentParser(add_help=False, description = 'POA2 - Analysis Conservancy and TMHMM')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-g", help= "-g: greater than or equal to (>=) in Sequence identity threshold from your analysis", type=str)
    mutually_exclusive_group.add_argument("-l", help= "-l: less-than sign (<) in Sequence identity threshold from your analysis", type=str)
    required.add_argument("-h", "--help", action ="help", default=argparse.SUPPRESS, help= "Print this help message.")
    required.add_argument("-d", help= "Directory for Results of the Epitope Conservancy analysis in CSV format", type=str, default='', required = "True")
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
 

    args = parser.parse_args()
    
    print('POA2 - Analysis Conservancy' + "\n")
    if (args.t == '') or (args.d == ''):
        raise Exception(f'ERROR: The following arguments are required: -t, -d')
    
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
    
    ConservancyAnalysis_DF = conservancyAnalysis.EpitConservAnalysis(args.t, type_symbol, args.m, args.imax, args.imin, args.d)   
    
    POA2_df = PepiTMHMM.tmhmmAnalysis(args, ConservancyAnalysis_DF)

    POA2_df.to_excel(f"{path}_{args.t}.xlsx", sheet_name = f'{type_symbol}{args.t}', startcol=0, index=False)

    if (args.rf) != None:
        getfastafile(path, args.rf, POA2_df)

    print("Analysis of Results Completed!\n")

if __name__ == '__main__':
    main()