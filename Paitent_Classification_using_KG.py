import pandas as pd
import numpy as np

# read CTD knowlege graph 
KG_ctd_gene_dis_required_genes = pd.read_csv("KG.csv", header = 0, sep=",")

KG_ctd_gene_dis_required_genes.head()

# Use Respiratory viral infection in DiseaseName (test data)
# Ampligraph baseline

import pickle
import os

def save_pickle(filename,data):
    output = open(filename, 'wb')
    # Pickle dictionary using protocol 0.
    pickle.dump(data, output)
    output.close()

def load_pickle(filename):
    pkl_file = open(filename, 'rb')
    data1 = pickle.load(pkl_file)
    pkl_file.close()
    return data1

if not os.path.exists('./pickle/unique_genes.pkl'):
    unique_genes = set(KG_ctd_gene_dis_required_genes['GeneSymbol'])
    save_pickle('./pickle/unique_genes.pkl',unique_genes)
else :
    unique_genes = load_pickle('./pickle/unique_genes.pkl')

# number of unique genes symbols
len(unique_genes)

if not os.path.exists('./pickle/unique_disease.pkl'):
    unique_disease = set(KG_ctd_gene_dis_required_genes['DiseaseName'])
    save_pickle('./pickle/unique_disease.pkl',unique_disease)
else :
    unique_disease = load_pickle('./pickle/unique_disease.pkl')

# number of unique diseases
len(unique_disease)

# read Gene expression data GSE73072 
Gene_expression_data_full = pd.read_csv("gene_expression_data.csv", header = 0, sep=",")

Gene_expression_data_full.head()

set(Gene_expression_data_full['Super_Subject_ID'])

# selecting only 15 subjects data which are infected
# Although the file already has required entries
infected_sub_ids = Gene_expression_data_full.loc[Gene_expression_data_full['Label'] == 1, ["Super_Subject_ID"]]
Gene_expression_data_full = Gene_expression_data_full.merge(infected_sub_ids)

infected_sub_ids.head()

Gene_expression_data_full.head()

j = 0
for i in range(0,len(Gene_expression_data_full),2):
    Gene_expression_data_full.loc[i,['Super_Subject_ID']] = j
    Gene_expression_data_full.loc[i+1,['Super_Subject_ID']] = j
    j += 1

Gene_expression_data_full.head()

Gene_expression_data_full.iloc[0:9,0:15]

Gene_expression_data_full.shape

print('Number of columns in the data',len(Gene_expression_data_full.columns))
print('Number of unique columns in the data',len(set(Gene_expression_data_full.columns)))

# removing duplicates
Gene_expression_data_full_no_duplicate = Gene_expression_data_full.loc[:,~Gene_expression_data_full.columns.duplicated()]

Gene_expression_data_full_no_duplicate.iloc[0:9,0:9]

Gene_expression_data_full_no_duplicate.shape

Gene_expression_data_full_no_duplicate.loc[Gene_expression_data_full_no_duplicate.shape[0] -1, "Super_Subject_ID"]

# Total number of subject
s = Gene_expression_data_full_no_duplicate.loc[Gene_expression_data_full_no_duplicate.shape[0] -1, "Super_Subject_ID"] + 1
print(s)

s_index = 9
e_index = Gene_expression_data_full_no_duplicate.shape[1]

def normalize(list_val):
    numerator =  list_val - min(list_val)
    denominator = max(list_val) - min(list_val)
    return numerator/denominator

if not os.path.exists('./pickle/Gene_expression_data_full_no_duplicate.pkl'):    
    for i in range(len(Gene_expression_data_full_no_duplicate)):
        Gene_expression_data_full_no_duplicate.iloc[i, s_index:e_index] = normalize(Gene_expression_data_full_no_duplicate.iloc[i, s_index:e_index])
    save_pickle('./pickle/Gene_expression_data_full_no_duplicate.pkl',Gene_expression_data_full_no_duplicate)
else :
    Gene_expression_data_full_no_duplicate = load_pickle('./pickle/Gene_expression_data_full_no_duplicate.pkl')

Gene_expression_data_full_no_duplicate.iloc[0:10,0:14]

# Taking out just gene expression columns
Gene_expression_data_inicial_columns = Gene_expression_data_full_no_duplicate.iloc[ :, 0:s_index]

Gene_expression_data_inicial_columns.head()

# taking only genes 
Gene_expression_data_only_genes = Gene_expression_data_full_no_duplicate.iloc[:, s_index:]

Gene_expression_data_only_genes.head()

# find common genes in KG and gene expression data
common_gene_KG_and_exp_data = Gene_expression_data_only_genes.loc[:,list(set(KG_ctd_gene_dis_required_genes.loc[:,'GeneSymbol']) 
                                   & set(Gene_expression_data_only_genes.columns))]

# print dimensions
common_gene_KG_and_exp_data.shape

# combine first 6 column again which has information except genes like label, title, time point, subject id etc.
Gene_expression_data_full_no_duplicate = pd.concat([Gene_expression_data_inicial_columns.reset_index(drop=True), common_gene_KG_and_exp_data], axis=1)
print(Gene_expression_data_full_no_duplicate.shape)

# print the first subject's data only 10 rows and 10 columns
Gene_expression_data_full_no_duplicate.iloc[0:10,0:15]

Gene_expression_data_full_no_duplicate.to_csv("Gene_expression_data_full_no_duplicate.csv")

pd.read_csv("Data_Unique_Disease.csv")

Data_Unique_Disease = pd.read_csv("Data_Unique_Disease.csv")

Data_Unique_Disease

import copy
Data_Unique_Disease_initial_weights = copy.deepcopy(Data_Unique_Disease)

############## User define input for the rang of top P genes and bottom Q genes ################# 
P_start = 200
P_end = 100 
P_step = 50
Q_start = 50
Q_end = 100
Q_step = 50

len(range(P_start,P_end,P_step))

PS = len(list(range(P_start,P_end,P_step))) + 1 
QS = len(list(range(Q_start,Q_end,Q_step))) + 1

list_fifty = [0 for i in range(PS*QS)]

list_fifty

Accuracy_matrix = pd.DataFrame({"P":list_fifty, "Q":list_fifty, "Acc_Top_1_Dis":list_fifty,
                 "Acc_Top_2_Dis":list_fifty, "Acc_Top_3_Dis":list_fifty,
                 "Acc_Top_4_Dis":list_fifty, "Acc_Top_5_Dis":list_fifty, 
                 "Acc_Top_10_Dis":list_fifty, "Max_Score":list_fifty, "Min_Score":list_fifty})

Accuracy_matrix

Acc_index = 0

Gene_expression_data_full_no_duplicate.shape[0]

s_index

predicted_test = Gene_expression_data_full_no_duplicate.iloc[:,:s_index]
predicted_test.head()

# # Ti - T1
# Gene_Transition_Matrix = Gene_expression_data_sub_l.iloc[1, ] - Gene_expression_data_sub_l.iloc[0, ]

# Gene_Transition_Matrix_top_p_Genes = Gene_Transition_Matrix.apply(lambda x:abs(x)).sort_values(ascending=False)
# Gene_Transition_Matrix_top_p_Genes = Gene_Transition_Matrix_top_p_Genes.iloc[0:50]


# print("show values of top 5 Gene_Transition_Matrix_top_p_Genes:")
# print(Gene_Transition_Matrix_top_p_Genes.iloc[0:5])

KG_ctd_gene_dis_required_genes['Disease_ID']

# abs(Gene_Transition_Matrix_top_p_Genes[0])

# Data_Unique_Disease.loc[Data_Unique_Disease["Disease_ID"] ==  1935, "Disease_Weight"] = abs(Gene_Transition_Matrix_top_p_Genes[0])

Data_Unique_Disease.loc[Data_Unique_Disease["Disease_ID"] ==  1935, "Disease_Weight"]

# Gene_Transition_Matrix_bottom_q_Genes = Gene_Transition_Matrix.apply(lambda x:abs(x)).sort_values(ascending=True)
# Gene_Transition_Matrix_bottom_q_Genes = Gene_Transition_Matrix_bottom_q_Genes.iloc[0:q+1]

# q

KG_ctd_gene_dis_required_genes.loc[KG_ctd_gene_dis_required_genes["GeneSymbol"] == 'SLC14A2', "Disease_ID" ]

s

len(Gene_expression_data_full_no_duplicate)

Data_Unique_Disease.loc[:,"Disease_Weight"] = [0 for i in range(len(Data_Unique_Disease))]

Data_Unique_Disease.head()

# pd.DataFrame(Gene_Transition_Matrix).sort_values(ascending=True, key = lambda x: abs(x))

# (Gene_expression_data_sub_l.iloc[1, ] - Gene_expression_data_sub_l.iloc[0, ]).sort_index(ascending=True, key=lambda x:abs(x))

# Gene_Transition_Matrix_top_p_Genes

Data_Unique_Disease.sort_values(by="Disease_Weight", ascending=False).loc[:20,:]

import time

Gene_Data_All_ti_prediction["Label"].values[k] == All_Sub_temp_prediction.iloc[k,i]

Commented out IPython magic to ensure Python compatibility.
%%time
for p in range(P_start,P_end,P_step):
    for q in range(Q_start,Q_end,Q_step):
        G_Max_Score = 0
        G_Min_Score = 0
        
        # total number of subjects
        s = Gene_expression_data_full_no_duplicate.loc[Gene_expression_data_full_no_duplicate.shape[0]-1, 
                                                       "Super_Subject_ID"] + 1
        
        predicted_info = Gene_expression_data_full_no_duplicate.iloc[:,:s_index]
        predicted_info['predicted_label_top_1'] = list(range(0,2*s))
        predicted_info['predicted_label_top_2'] = list(range(0,2*s))
        predicted_info['predicted_label_top_3'] = list(range(0,2*s))
        predicted_info['predicted_label_top_4'] = list(range(0,2*s))
        predicted_info['predicted_label_top_5'] = list(range(0,2*s))
        predicted_info['predicted_label_top_10'] = list(range(0,2*s))
        
        Gene_Data_All_ti_prediction = Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate["Time_Point_Adjusted"] == 72].iloc[:,:s_index]
        
        Gene_Data_All_ti_prediction['predicted_label_top_1'] = list(range(0,s))
        Gene_Data_All_ti_prediction['predicted_label_top_2'] = list(range(0,s))
        Gene_Data_All_ti_prediction['predicted_label_top_3'] = list(range(0,s))
        Gene_Data_All_ti_prediction['predicted_label_top_4'] = list(range(0,s))
        Gene_Data_All_ti_prediction['predicted_label_top_5'] = list(range(0,s))
        Gene_Data_All_ti_prediction['predicted_label_top_10'] = list(range(0,s))
        
        All_Sub_temp_prediction = pd.DataFrame({"Top_1":list(range(s)), "Top_2":list(range(s)),
                                                "Top_3":list(range(s)), "Top_4":list(range(s)),
                                                "Top_5":list(range(s)), "Top_10":list(range(s))})
        #debug
#         print(All_Sub_temp_prediction)
        for l in range(s):            
            Gene_expression_data_sub_l = Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate["Super_Subject_ID"] == l]
            Gene_expression_data_sub_l = Gene_expression_data_sub_l.iloc[:,s_index:]
            
            print("########################## New subject computation start from here ############################")
            print("Subject id is:")
            print(l)
            
            # for each time point initialize again
#             Data_Unique_Disease = Data_Unique_Disease_initial_weights
            
            Data_Unique_Disease.loc[:,"Disease_Weight"] = [0 for i in range(len(Data_Unique_Disease))]
            
            # Ti - T1
            Gene_Transition_Matrix = Gene_expression_data_sub_l.iloc[1, ] - Gene_expression_data_sub_l.iloc[0, ]
            
            Gene_Transition_Matrix_top_p_Genes = Gene_Transition_Matrix.loc[list(Gene_Transition_Matrix.apply(lambda x: abs(x)).sort_values(ascending=False).index)]
#             Gene_Transition_Matrix_top_p_Genes = Gene_Transition_Matrix.sort_values(ascending=False, key=lambda x:abs(x))
            Gene_Transition_Matrix_top_p_Genes = Gene_Transition_Matrix_top_p_Genes.iloc[0:p]


            print("show values of top 5 Gene_Transition_Matrix_top_p_Genes:")
            print(Gene_Transition_Matrix_top_p_Genes.iloc[0:5])
            
            # loop of j for number of genes
            for j in range(p):
                Disease_IDs = KG_ctd_gene_dis_required_genes[KG_ctd_gene_dis_required_genes["GeneSymbol"] == Gene_Transition_Matrix_top_p_Genes.keys()[j]]["Disease_ID"]
                # loop for every disease id
                for k in range(len(Disease_IDs)):
                    Data_Unique_Disease.loc[Data_Unique_Disease["Disease_ID"] ==  Disease_IDs.values[k], "Disease_Weight"] += abs(Gene_Transition_Matrix_top_p_Genes[j])
            
            # Penalty
            Gene_Transition_Matrix_bottom_q_Genes = Gene_Transition_Matrix.loc[list(Gene_Transition_Matrix.apply(lambda x: abs(x)).sort_values(ascending=True).index)]
            
#             Gene_Transition_Matrix_bottom_q_Genes = Gene_Transition_Matrix.sort_values(ascending=True, key=lambda x:abs(x))
            Gene_Transition_Matrix_bottom_q_Genes = Gene_Transition_Matrix_bottom_q_Genes.iloc[0:q]
            
            print("show values of Gene_Transition_Matrix_bottom_q_Genes:")
            print(Gene_Transition_Matrix_bottom_q_Genes[0:5])
            
            
            # loop of j for number of bottom genes
            for j in range(q):
                Disease_IDs = KG_ctd_gene_dis_required_genes.loc[KG_ctd_gene_dis_required_genes["GeneSymbol"] == Gene_Transition_Matrix_bottom_q_Genes.keys()[j], "Disease_ID" ]
                
                for k in range(len(Disease_IDs)):
                    Data_Unique_Disease.loc[Data_Unique_Disease["Disease_ID"] ==  Disease_IDs.values[k], "Disease_Weight"] -= abs(Gene_Transition_Matrix_top_p_Genes[0])
                    Data_Unique_Disease.loc[Data_Unique_Disease["Disease_ID"] ==  Disease_IDs.values[k], "Disease_Weight"] -= abs(Gene_Transition_Matrix_bottom_q_Genes[j])
            
            Max_Score = Data_Unique_Disease.loc[:,'Disease_Weight'].sort_values(ascending=False)[0]
            Min_Score = Data_Unique_Disease.loc[:,'Disease_Weight'].sort_values(ascending=True)[0]
            
            # create file name to write data into csv file
            file_name = "Disease_Weight_Sub_" + str(l) + "_p_" + str(p) + "_q_" + str(q) + ".csv"
            
            # write data into csv file
            Data_Unique_Disease.loc[:,['Disease_Weight']].sort_values(ascending=False, by="Disease_Weight").to_csv(file_name)
            
            print("Subject id is:")
            print(l)
  
            print("Value of p is :")
            print(p)
  
            print("Value of q is :")
            print(q)
  
            print("This subject at this time point has following label:")
            print(Gene_expression_data_full_no_duplicate.loc[Gene_expression_data_full_no_duplicate["Super_Subject_ID"] == l, 'Label'].reset_index(drop=True)[1])
            
#             disease_idx = list(Data_Unique_Disease.loc[:,'Disease_Weight'].apply(lambda x:abs(x)).sort_values(ascending=False).index)
            # loop for how many top disease you want to look for Acc calc
            for i in range(6):
                if i < 6 :
#                     Top_Disease_Names = Data_Unique_Disease.loc[disease_idx[:i], "Disease_Name"]
#                     Top_Disease_Names = Data_Unique_Disease.loc[:i, ["Disease_Name", "Disease_Weight"]].sort_values(by= "Disease_Weight",ascending=False)
                    Top_Disease_Names = Data_Unique_Disease.sort_values(by= "Disease_Weight",ascending=False).loc[:i,:]
                else:
#                     Top_Disease_Names = Data_Unique_Disease.loc[disease_idx[:10], "Disease_Name"]
#                     Top_Disease_Names = Data_Unique_Disease.loc[:10, ["Disease_Name", "Disease_Weight"]].sort_values(by= "Disease_Weight",ascending=False)
                    Top_Disease_Names = Data_Unique_Disease.sort_values(by= "Disease_Weight",ascending=False).loc[:10,:]
                
                if "Respiratory Viral Infection" in Top_Disease_Names.loc[:,"Disease_Name"].values:
#                     check index
                    predicted_info.loc[predicted_info["Super_Subject_ID"] == l, predicted_info.columns[s_index + i]] = 1
#                     check df
                    Gene_Data_All_ti_prediction.loc[Gene_Data_All_ti_prediction["Super_Subject_ID"] == l , 
                                                    s_index + i] = 1
                    All_Sub_temp_prediction.iloc[l,i] = 1
                else:
                    predicted_info.loc[predicted_info["Super_Subject_ID"] == l, predicted_info.columns[s_index + i]] = 0
                    Gene_Data_All_ti_prediction.loc[Gene_Data_All_ti_prediction["Super_Subject_ID"] == l , 
                                                    s_index + i] = 0
                    All_Sub_temp_prediction.iloc[l,i] = 0
                
                print("Predicted label using top ", i, " disease:")
                print(All_Sub_temp_prediction.iloc[l,i])
                
            print("Top 20 Diseases are:")

            print(Data_Unique_Disease.sort_values(by= "Disease_Weight",ascending=False).loc[:20,:])
            
#             print(Data_Unique_Disease.loc[disease_idx[:20],:])

            print("Top 5 Genes are:")

            print(Gene_Transition_Matrix_top_p_Genes[:5])
            
            print("############################################################################")
            
            if G_Min_Score > Min_Score:
                G_Min_Score = Min_Score
            
            else:
                G_Max_Score = Max_Score
                
        file_name_1 = "predicted_info_60hr_Sub_"+ str(l) + "_p_" + str(p) + "_q_" + str(q) + ".csv"
        
        Gene_Data_All_ti_prediction.to_csv(file_name_1)
        
        
        for i in range(6):
            
            hit = 0
            
            for k in range(All_Sub_temp_prediction.shape[0]):
                # debug
                print(Gene_Data_All_ti_prediction["Label"].values[k])
                print(All_Sub_temp_prediction.iloc[k,i])
                if Gene_Data_All_ti_prediction["Label"].values[k] == All_Sub_temp_prediction.iloc[k,i]:
                    hit += 1
            # if hit is zero: index issue
            # debug
            print(hit)       
            Acc = hit/All_Sub_temp_prediction.shape[0]
            print("Accuracy using top ", i, "disease is:")
            print(Acc)
            # debug
            Accuracy_matrix.iloc[Acc_index, 2+i] = Acc
            
        Accuracy_matrix.loc[Acc_index, "P"] = p
        Accuracy_matrix.loc[Acc_index, "Q"] = q
        Accuracy_matrix.loc[Acc_index, "Max_Score"] = G_Max_Score
        Accuracy_matrix.loc[Acc_index,"Min_Score"] = G_Min_Score
        
        Acc_index += 1
        
        # create file name to write data into csv file
        file_name_1 = "Accuracy_matrix_60hr_Sub"+ str(l) + "_p_"+ "p" + "_q_" + "q" + ".csv"
        
        Accuracy_matrix.to_csv(file_name_1)
        
        print("Accuracy Matrix:")
        print(Accuracy_matrix.iloc[1:Acc_index,:])      
        
       
