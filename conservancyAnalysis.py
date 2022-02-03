#Get the list of Conservancy Analysis results filenames
def map_EpConservFiles(directory):
    import os
    EpConservFiles = []
    for folders, subfolders, files in os.walk(directory):
        for file in files:
            path_file = (os.path.join(folders, file))
            path_file = str(path_file)
            if '.csv' in file:
                EpConservFiles.append(path_file)
    return (EpConservFiles)
#Data Screening of Epitope Conservancy Analysis results 
def EpitConservAnalysis(ID_threshold, symbol, seq_match, max_ID, min_ID, directory):
    import pandas as pd
    filesList = map_EpConservFiles(directory)
    allSp_ConservDF = pd.DataFrame() #create a unique dataframe for results
    data = ["Epitope #","Epitope name","Epitope sequence","Epitope length","Percent of protein sequence matches at identity <= 100%","Minimum identity","Maximum identity","View details"]
    allSp_ConservDF = pd.DataFrame(columns=data)
    for file in filesList:
        Sp_ConservDF = pd.read_csv(file) #open conservancy analysis result (csv) as dataframe
        allSp_ConservDF = pd.merge(allSp_ConservDF, Sp_ConservDF, how = 'outer') #merge conservancy analysis result dataframe with general dataframe for all results 
    allSp_ConservDF['Percent of protein sequence matches at identity <= 100%'] = allSp_ConservDF['Percent of protein sequence matches at identity <= 100%'].map(lambda x: x.replace('%', ''))
    allSp_ConservDF['Minimum identity'] = allSp_ConservDF['Minimum identity'].map(lambda x: x.lstrip('').rstrip('%')) #removing the '%' character from the line 
    allSp_ConservDF['Maximum identity'] = allSp_ConservDF['Maximum identity'].map(lambda x: x.lstrip('').rstrip('%'))
    allSp_ConservDF['Minimum identity'] = allSp_ConservDF['Minimum identity'].apply(float) #converting line to float type
    allSp_ConservDF['Maximum identity'] = allSp_ConservDF['Maximum identity'].apply(float)
    #applying filters 
    allSp_ConservDF_idMin = allSp_ConservDF[allSp_ConservDF['Minimum identity'] >= min_ID]
    if allSp_ConservDF_idMin.empty:#checking if the dataframe is empty 
        return(allSp_ConservDF_idMin)
    allSp_ConservDF_idMax =  allSp_ConservDF_idMin[ allSp_ConservDF_idMin['Maximum identity'] <= max_ID]
    if allSp_ConservDF_idMax.empty:#checking if the dataframe is empty 
        return(allSp_ConservDF_idMax)
    seqMatchSplit = allSp_ConservDF_idMax["Percent of protein sequence matches at identity <= 100%"].str.split(" ", n = 1, expand = True)
    allSp_ConservDF_idMax["Percent"]= seqMatchSplit[0] 
    allSp_ConservDF_idMax["Reason"]= seqMatchSplit[1] 
    allSp_ConservDF_idMax["Percent"] = allSp_ConservDF_idMax["Percent"].apply(float)
    allSp_ConservDF_final = allSp_ConservDF_idMax[allSp_ConservDF_idMax["Percent"] >= seq_match] 

    new_colNames = {
        "Percent of protein sequence matches at identity <= 100%":f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}",
        "Minimum identity": "Minimum identity(%)",
        "Maximum identity": "Maximum identity(%)"
        }
        
    allSp_ConservDF_final = allSp_ConservDF_final.rename(columns = new_colNames, inplace = False)        
    final_Conserv_df = allSp_ConservDF_final[["Epitope name","Epitope sequence","Epitope length",f"Protein(s) sequence match(es) at Sequence identity threshold {symbol}{ID_threshold}","Minimum identity(%)","Maximum identity(%)"]]
        
    return(final_Conserv_df)