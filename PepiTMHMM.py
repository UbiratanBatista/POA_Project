#checking that entries are not empty
def notEmptyValidate(file):
    import os
    result = os.stat(file).st_size==0
    return result

#Check if epitope position values are possible 
def rangeVerify(valor, lenghtProt):
    if int(valor) == 0 or int(valor) < 0:
        return False
    elif int(valor) > int(lenghtProt):
        return False
    else:
        return True

def tmhmmAnalysis(args, dataframe):
    import TMHMM
    import warnings
    import pandas as pd
    import Bio
    from Bio import SeqIO

    print("TMHMM - Epitope Analysis v1.0" + "\n")
    #check the entries
    if notEmptyValidate(args.f):
        raise Exception("failed because protein fasta is empty.")
        
    #Printing fasta protein sequences 
    for seq_record in SeqIO.parse(args.f, "fasta"):
        print("Completed Sequence Analysis:"+ seq_record.id)       

    #Update columns in dataframe of the Conservancy Analysis results     
    dataframe = dataframe.assign(Portion_Outside = "-")
    dataframe = dataframe.assign(Portion_TM = "-")
    dataframe = dataframe.assign(Portion_Inside = "-")

    #TMHMM prediction for fasta protein sequences 
    ids_TMHMM_list, seqs_TMHMM_list = TMHMM.pyTMHMMpredict(args.f)

    x = 0 #Counter for dataframe rows 

    Out_portion_list = []; TM_portion_list = []; Ins_portion_list = []

    #Analyzing each of the DF epitopes 
    for indice_line in range(len(dataframe["Epitope sequence"])):
        #Organize and separate epitope information 
        idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
        idt = idt_epitope.upper()
        idt = idt.split('_')
        if len(idt) == 5:
            EpitopeVirus = idt[0]
            init_pos = idt[-2]
            fin_pos = idt[-1]
        else:
            raise Exception("failed because id epitope " + idt_epitope + " is not correctly formatted. {} infos were found".format(str(idt)))
        #Get the epitope sequence
        epitope = dataframe.iloc[indice_line]["Epitope sequence"]
        epitope = epitope.lower() #lowercase standardization
        
        #Searching for epitope and protein matches 
        for seq_record in SeqIO.parse(args.f, "fasta"):
            #format the protein identification and sequence 
            seq_polyprot = str(seq_record.seq)
            seq_polyprot = seq_polyprot.lower() #lowercase standardization
            seq_polyprot_id = seq_record.id.upper() #uppercase standardization
            #Checking if the sequence in question belongs to the virus of the analyzed epitope 
            if EpitopeVirus not in seq_polyprot_id:
                continue
            #the sequence virus species is the same as the epitope virus species 
            else:
                #Ensuring epitope sequence and protein sequence matching 
                if epitope not in seq_polyprot:
                    warnings.warn("Epitope "+ idt_epitope + " are not in the protein(s) fasta file, check your sequences output")
                    continue
                else:
                    #verify feasibility of epitope position data 
                    start_pos_verif = rangeVerify(init_pos, len(seq_record))
                    end_pos_verif = rangeVerify(fin_pos, len(seq_record))
                    if (start_pos_verif != True) or (end_pos_verif != True):
                        raise Exception("failed because epitope position in" + idt_epitope + "does not exist.")
                    #evaluate the correspondence between the epitope and the protein after the prediction of the TMHMM
                    for x in range(len(seqs_TMHMM_list)):
                        if seq_polyprot_id == ids_TMHMM_list[x]:
                            init = seq_polyprot.find(epitope)
                            if init != -1:
                                epit_lenght = len(epitope)
                                #Calculate the percentage of amino acid residues in each of the membrane topology classifications (outer, transmembrane, and inner) 
                                Out_portion, TM_portion, Ins_portion = TMHMM.epitTMHMMcaract (seqs_TMHMM_list[x], init, epit_lenght)
                    
                    Out_portion_list.append(Out_portion); TM_portion_list.append(TM_portion); Ins_portion_list.append(Ins_portion)
    
    dataframe["Portion_Outside"] = Out_portion_list
    dataframe["Portion_TM"] = TM_portion_list
    dataframe["Portion_Inside"] = Ins_portion_list
    
    return (dataframe)