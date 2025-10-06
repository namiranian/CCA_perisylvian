import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
import csv

def get_files(directory_path):
    """
    Returns a list of CSV file names in the specified directory using os.listdir().
    """
    csv_files = []
    for filename in os.listdir(directory_path):
        if filename.endswith(".csv"):
            csv_files.append(filename)
    return csv_files

def savefeatures(feat, name, y=None, directoy='C:/SF_data/'):
    pheno = y.astype(float)  # phenosmoother(y, interval=1, mode='reg')
    x = regoutfunc(feat, pheno)
    filesave(x, directoy + name + '.csv')

def regoutfunc(x, pheno):
    """  Regout the pheno variable from x """
    paramsx = np.zeros((2, x.shape[1]))
    pvalsx = np.zeros((2, x.shape[1]))
    exog = sm.add_constant(pheno)
    x_res=np.zeros(x.shape)
    for i in range(x.shape[1]):
        model = sm.OLS(x[:, i], exog)
        results = model.fit()
        paramsx[:, i] = results.params
        pvalsx[:, i] = results.pvalues
        # for attributeIndex in range (0, 1):
        #         pvalsx[:,i] = results.pvalues[attributeIndex]
        x_res[:, i] = results.resid

    return x_res

def filesave(x, name):
    try:
        row = np.zeros((1, x.shape[1]))
    except:
        row = np.zeros((1, 1))
        x = x.reshape((x.shape[0], 1))
    x = np.vstack([row, x])
    with open(name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        try:
            writer.writerows(x)
        except:
            writer.writerows(map(lambda i: [i], x))


#read individual metrics
main_dir=os.getcwd()+'/input_data/'
file_list=get_files(main_dir+'/metrics/')

subjct_list = pd.read_excel(main_dir+'subject_list.xlsx')
sub_age=subjct_list['scan_age'].astype(float)

for file in file_list:
    feature=pd.read_csv(main_dir+'/metrics/'+file).to_numpy()
    feat_regout = regoutfunc(x=feature, pheno=sub_age)
    filesave(x=feat_regout, name=main_dir+'/age_reg_out_metrics/'+file)

