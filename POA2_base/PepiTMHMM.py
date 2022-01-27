#verificando se as entradas não estão vazias
def notEmptyValidate(file):
    import os
    result = os.stat(file).st_size==0
    return result

#Verificar se os valores da posição do epitopo são possíveis
def rangeVerify(valor, lenghtProt):
    if int(valor) == 0 or int(valor) < 0:
        return False
    elif int(valor) > int(lenghtProt):
        return False
    else:
        return True

#OBS. O nome da proteina deve aparecer no cabeçalho do fasta proteína da mesma que forma que no cabeçalho dos epitopos
#EX: Se aparece capsid no fasta da proteina, deve aparecer capsid no fasta dos epitopos, não proteina C
def tmhmmAnalysis(args, dataframe):
    import TMHMM
    import warnings
    import pandas as pd
    import Bio
    from Bio import SeqIO

    print("TMHMM - Epitope Analysis v1.0" + "\n")
    #Verificando as entradas
    if notEmptyValidate(args.f):
        raise Exception("failed because protein fasta is empty.")
        
    #Reconhecimento das sequencias proteicas do fasta
    for seq_record in SeqIO.parse(args.f, "fasta"):
        print("Completed Sequence Analysis:"+ seq_record.id)       

    #Atualizando a planilha do Analysis Conservancy    
    dataframe = dataframe.assign(Portion_Outside = "-")
    dataframe = dataframe.assign(Portion_TM = "-")
    dataframe = dataframe.assign(Portion_Inside = "-")

    #fazendo a predição do TMHMM para as sequencias das proteínas no fasta
    ids_TMHMM_list, seqs_TMHMM_list = TMHMM.pyTMHMMpredict(args.f)

    x = 0 #Contador para as linhas da planilha

    Out_portion_list = []; TM_portion_list = []; Ins_portion_list = []

    #Lendo cada um dos epitopos do DF
    for indice_line in range(len(dataframe["Epitope sequence"])):
        #Organizando e separando informações do epitopo
        idt_epitope = dataframe.iloc[indice_line]["Epitope name"]
        idt = idt_epitope.upper()
        idt = idt.split('_')
        if len(idt) == 5:
            EpitopeVirus = idt[0]
            EpitopeProtein = idt[1]
            method = idt [2]
            init_pos = idt[-2]
            fin_pos = idt[-1]
        else:
            raise Exception("failed because id epitope " + idt_epitope + " is not correctly formatted. {} infos were found".format(str(idt)))
        #Recuperando a sequencia do Epitopo
        epitope = dataframe.iloc[indice_line]["Epitope sequence"]
        epitope = epitope.lower() #padronizando em minusculo
        
        #Fazendo a leitura por proteina (math do epitopo com alguma das proteinas)
        for seq_record in SeqIO.parse(args.f, "fasta"):
            #formatando a sequencia da proteina
            seq_polyprot = str(seq_record.seq)
            seq_polyprot = seq_polyprot.lower() #padronizando em minusculo
            seq_polyprot_id = seq_record.id.upper() #padronizando o ID em maiusculo
            #Conferindo se a sequencia em questão pertence ao virus do epitopo analizado
            if EpitopeVirus not in seq_polyprot_id:
                continue
            #virus é o mesmo do epitopo
            else:
                #Conferindo o match do epitopo com a sequencia da proteina
                if epitope not in seq_polyprot:
                    warnings.warn("Epitope "+ idt_epitope + " are not in the protein(s) fasta file, check your sequences output")
                    continue
                else:
                     #testando a viabilidade dos dados de posição do epitopo
                    start_pos_verif = rangeVerify(init_pos, len(seq_record))
                    end_pos_verif = rangeVerify(fin_pos, len(seq_record))
                    if (start_pos_verif != True) or (end_pos_verif != True):
                        raise Exception("failed because epitope position in" + idt_epitope + "does not exist.")
                    
                    for x in range(len(seqs_TMHMM_list)):
                        if seq_polyprot_id == ids_TMHMM_list[x]:
                            init = seq_polyprot.find(epitope)
                            if init != -1:
                                epit_lenght = len(epitope)
                                Out_portion, TM_portion, Ins_portion = TMHMM.epitTMHMMcaract (seqs_TMHMM_list[x], init, epit_lenght)
                    
                    Out_portion_list.append(Out_portion); TM_portion_list.append(TM_portion); Ins_portion_list.append(Ins_portion)
    
    dataframe["Portion_Outside"] = Out_portion_list
    dataframe["Portion_TM"] = TM_portion_list
    dataframe["Portion_Inside"] = Ins_portion_list
    
    return (dataframe)