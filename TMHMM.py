#OBS: o pyTMHMM só recebe a cadeia da proteina(sem o cabeçalho do fasta)
def pyTMHMMpredict(file):
    import pyTMHMM
    from Bio import SeqIO
    id_seqs_list = []
    prediction_TMHMM_list = []
    for seq_record in SeqIO.parse(file, "fasta"):
        id_seq = (seq_record.id).upper()
        seq = str(seq_record.seq)
        #Apenas a sequencia de aminoácidos para a predição        
        annotation = pyTMHMM.predict(seq, compute_posterior=False)
        annotation = annotation.lower()
        id_seqs_list.append(id_seq)
        prediction_TMHMM_list.append(annotation)
    return(id_seqs_list, prediction_TMHMM_list)

def epitTMHMMcaract (seq_TMHMM, initial_pos, lenght): #Contabilizando as posições dos aminoácidos
    epitope_TMHMM = seq_TMHMM[initial_pos:(initial_pos+lenght)]
    out_count = epitope_TMHMM.count('o')
    tm_count = epitope_TMHMM.count('m')
    ins_count = epitope_TMHMM.count('i')
    Outporct = out_count/lenght
    TMporct = tm_count/lenght
    Insporct = ins_count/lenght
    return (round(Outporct, 4), round(TMporct, 4), round(Insporct, 4))