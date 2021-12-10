import pandas as pd
import argparse
import prepianResultsforConservancyAnalysis
import os
import PrEpiAn

#Write epitope analysis report within the pipeline 
def writereport(dataframe, path):
    import math
    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    with open (f"{path}/Analysis_report.txt", "w") as report:
        report.write("Analysis Report: Pipeline XXX\n\n")
        report.write(dt_string+"\n\n")
        organisms = list(dataframe["Specie"].unique())#Number of species
        report.write(f"#Species present in the analysis:\n{len(organisms)}\n")
        proteinsInAnalysis = 0 #Number of proteins
        percentGreater = 0; percentfewer = 0
        Prot_greaterEpit = []; greater_Epitopes = -(math.inf)
        Prot_fewerEpit = []; fewer_Epitopes = (math.inf)
        Total_epitopes = len(dataframe["Peptide Sequence"])
        for organism in organisms:
            df_organism = dataframe.loc[(dataframe["Specie"] == organism)] #dataframe analysis for each species
            organismProteins = list(df_organism["Protein"].unique())
            proteinsInAnalysis += len(organismProteins)
            for protein in organismProteins:
                df_protein = df_organism.loc[(df_organism["Protein"] == protein)] #dataframe analysis by protein 
                epitopesinProtein = len(df_protein["Peptide Sequence"])
                percent_epitopesinProtein = (epitopesinProtein*100)/Total_epitopes
                #proteins with fewer epitopes 
                if epitopesinProtein < fewer_Epitopes:
                    fewer_Epitopes = len(df_protein["Peptide Sequence"])
                    percentfewer = percent_epitopesinProtein
                    Prot_fewerEpit = []
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == fewer_Epitopes:
                    Prot_fewerEpit.append(f"{protein}/{organism}")
                else:
                    pass
                #proteins with more epitopes 
                if epitopesinProtein > greater_Epitopes:
                    greater_Epitopes = len(df_protein["Peptide Sequence"])
                    percentGreater = percent_epitopesinProtein
                    Prot_greaterEpit = []
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                elif epitopesinProtein == greater_Epitopes:
                    Prot_greaterEpit.append(f"{protein}/{organism}")
                else:
                    pass
        report.write(f"#Proteins present in the analysis:\n{proteinsInAnalysis}\n")
        report.write("#Epitopes\n\n")
        methods = list(dataframe["Method"].unique())
        for method in methods:
            df_method = dataframe.loc[(dataframe["Method"] == method)]
            epitopesInMethod = len(df_method["Peptide Sequence"])
            report.write(f"{method}:\t{epitopesInMethod}\n")
        report.write(f"\nTotal:\t{Total_epitopes}\n\n")
        report.write("#Proteins with the greatest amount of epitopes\n")
        report.write(f"{Prot_greaterEpit}: {greater_Epitopes} ({percentGreater:.2f}%)\n\n")
        report.write("#Proteins with the least amount of epitopes\n")
        report.write(f"{Prot_fewerEpit}: {fewer_Epitopes} ({percentfewer:.2f}%)\n\n")  
        average_length = 0
        for epitope in dataframe["Peptide Sequence"]:
            average_length += len(epitope)
        average_length /= Total_epitopes
        report.write("#Average length of epitopes:\n")
        report.write(f"{int(average_length)} aa.")
    
def main():
    
    parser = argparse.ArgumentParser(add_help=False)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments (at least one prediction method)')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "Print this help message.")
    required.add_argument("-n", help= "NetCTL prediction in HTML format", type=str, default='')
    required.add_argument("-m", help= "Directory for MHC-II Binding Predictions in HTML format", type=str, default='')
    optional.add_argument("-mhla", help="HLA-type allele of the predicted MHC-II Binding Predictions Epitopes  (default = DR)", type=str, default = 'DR')
    optional.add_argument("-mic", help="IC50 threshold - NN_align 2.3 (default = 50), High binding peptides", type=int, default = 50) 
    required.add_argument("-b", help= "Bepipred-2.0 prediction in JSON format", type=str, default='')
    optional.add_argument("-bmin", help="Min. Length (MERS) of the predicted Bepipred Epitopes  (default = None)", type=int, default=0)
    optional.add_argument("-bmax", help="Max. Length (MERS) of the predicted Bepipred Epitopes  (default = None)", type=int, default=0)
    required.add_argument("-p", help= "PAP-IMED prediction in txt format", type=str, default='')
    optional.add_argument("-pmin", help="Min. Length (MERS) of the predicted PAP-IMED Epitopes  (default = None)", type=int, default=0)
    optional.add_argument("-pmax", help="Max. Length (MERS) of the predicted PAP-IMED Epitopes  (default = None)", type=int, default=0)
    required.add_argument("-x", help= "Prediction in others web servers (in fasta format)", type=str, default='')
    optional.add_argument("-xmin", help="Min. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    optional.add_argument("-xmax", help="Max. Length (MERS) of the results predicted in others web servers  (default = None)", type=int, default=0)
    required.add_argument("-d", help= "Directory for analysis results",required = "True", type=str, default='')
    required.add_argument("-f", help= "Polyproteins/proteins used in predictions analysis in fasta format", required = "True", type=str, default='')
    optional.add_argument("-e", help= "Analysis results in .xlsx format (y/n)", type=str, default='n')
    
    args = parser.parse_args()
    print("PrEpiAn: Predicted Epitopes Analyzer v1.0" + "\n")
    
    #At least one argument must be selected.
    if (args.b == '') and (args.p == '') and (args.n == '') and (args.m == '') and (args.x == ''):
        raise Exception(f'ERROR: Some the following arguments are required: -n, -m, -b, -p')
    #organize prediction results
    results_df = PrEpiAn.runningPrEpiAn(args)
  
    #Analysis Report
    writereport(results_df, args.d)

    #Organize data for Epitope Conservancy Analysis
    prepianResultsforConservancyAnalysis.filesforEptConsAnalysis(args.f, results_df, args.d)

    print("Data collected and analyzed...\n")
    print(".\n.\n.\n.\n.\n.\n.\n")
    print("Finished")


if __name__ == '__main__':
    main()