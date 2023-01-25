import streamlit as st
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
from tensorflow import keras
from tensorflow.keras import layers
from keras import regularizers
import numpy as np
import pandas as pd
from statistics import mean
#from st_radial import st_radial
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors
#import datetime
#from streamlit_option_menu import option_menu
#import streamlit.components.v1 as components
#import psycopg2
#from st_aggrid import AgGrid,GridUpdateMode
#from st_aggrid.grid_options_builder import GridOptionsBuilder

#import psycopg2.extras
import pandas as pd

#################################### Functions #########################
class MACCS:
    def __init__(self, smiles):
        self.smiles = smiles
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]

    def compute_MACCS(self, name):
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for i in range(len(self.mols)):
            ds = list(MACCSkeys.GenMACCSKeys(self.mols[i]).ToBitString())
            MACCS_list.append(ds)
        #df = pd.DataFrame(MACCS_list,columns=header)
        #df.insert(loc=0, column='smiles', value=self.smiles)
        #df.to_csv(name[:-4]+'_MACCS.csv', index=False)
        return MACCS_list

def extract_MACCSds(column,  from_smiles=True):

    maccs_descriptor = MACCS(column)
    feature_MACCS = maccs_descriptor.compute_MACCS('MACCSds.csv')
    return np.array(feature_MACCS)

def get_fingerprint(smiles):
    '''generates fingerprint dataframe given the smiles dataframe with 'SMILES' column containing the smiles. Also saves bit information needed for recovering the substructures later'''
    bit_infos = []
    rdkit_molecules = []
    for x in smiles['SMILES']:
        x = x.replace(' ', '')
        rdkit_molecules.append(Chem.MolFromSmiles(x))
    rdkit_fingerprint = []
    i = 0
    for mol in rdkit_molecules:
        print(i)
        bit_info = {}
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=4, nBits=10000, bitInfo=bit_info).ToBitString()
        bit_infos.append(bit_info)
        rdkit_fingerprint.append(fp)
        i = i + 1
    fingerprint_df = pd.DataFrame([np.array(list(x)).astype(int) for x in rdkit_fingerprint])
    return fingerprint_df, bit_infos
def predict_tox (string_smile):
    beta = 0.1
    neurons = 100
    hidden_layers = 2
    LR = 0.001
    epochs = 100

    inputs = keras.Input(shape=(10000,), name='digits')
    x = layers.Dense(neurons, activation='relu', kernel_initializer='random_uniform', bias_initializer='zeros',
                     kernel_regularizer=regularizers.l2(beta), name='dense_1')(inputs)
    x = layers.Dense(neurons, activation='relu', kernel_initializer='random_uniform', bias_initializer='zeros',
                     kernel_regularizer=regularizers.l2(beta), name='dense_2')(x)
    outputs = layers.Dense(2, activation='softmax', name='predictions')(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    ######
    # The new smile to predict
    ######
    smiles = list()
    smiles.append(string_smile) #'C1=C(N=CS1)CC(C(=O)O)N'
    df = pd.DataFrame(smiles, columns=['SMILES'])
    # generate fingeprint
    fingerprint, bit_infos = get_fingerprint(df)
    fingerprint = np.array(fingerprint)
    ########
    prefix = 'tox '  # xenobiotic1010 model10
    dir_models = 'models_tox'
    prefix = dir_models + '/' + prefix

    list_xenobiotic = list()
    list_nonXenobiotic = list()
    list_prob = list()
    #for f in range(10):
    model.load_weights(prefix + "model.h5")
    model.compile(optimizer=keras.optimizers.Adam(lr=LR),
                  loss='sparse_categorical_crossentropy',
                  metrics=['sparse_categorical_accuracy'])
    # print("Model loaded==>"+str(f+1))

    y_pred = model.predict(fingerprint)
    #st.write(str(f + 1), '===>', y_pred)
    list_prob.append(y_pred[0])
    if y_pred[0][0] >= 0.50:
        list_xenobiotic.append(y_pred[0][0])
    else:
        list_nonXenobiotic.append(y_pred[0][1])
   # if y_pred[0][0] > y_pred[0][1]:

    #st.write(len(list_nonXenobiotic))
    #if len(list_xenobiotic) >= len(list_nonXenobiotic):
        #formatted_string = "{:.2%}".format(mean(list_xenobiotic))
        #st.write(smiles[0], '\nis Xenobiotics with a probability of ==>', formatted_string, '\nThis is the decision of ',
              #len(list_xenobiotic), ' out of 10 classifiers')
        #st_radial('You selected', mean(list_xenobiotic))
    #else:
        #formatted_string = "{:.2%}".format(mean(list_nonXenobiotic))
        #st.write(smiles[0], '\nNON-Xenobiotics with a probability of ==>',formatted_string,
              #'\nThis is the Decision of ', len(list_nonXenobiotic), ' out of 10 classifiers')
    df_prob = pd.DataFrame(list_prob, columns=['Toxic', 'Non-Toxic'])
    df_prob = df_prob.style.format({
        'Toxic': '{:,.2%}'.format,
        'Non-Toxic': '{:,.2%}'.format,
    })
    return df_prob

def predict_xeno (string_smile, batch):
    beta = 0.1
    neurons = 100
    hidden_layers = 2
    LR = 0.001
    epochs = 100

    inputs = keras.Input(shape=(10167,), name='digits')
    x = layers.Dense(neurons, activation='relu', kernel_initializer='random_uniform', bias_initializer='zeros',
                     kernel_regularizer=regularizers.l2(beta), name='dense_1')(inputs)
    x = layers.Dense(neurons, activation='relu', kernel_initializer='random_uniform', bias_initializer='zeros',
                     kernel_regularizer=regularizers.l2(beta), name='dense_2')(x)
    outputs = layers.Dense(2, activation='softmax', name='predictions')(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    ######
    # The new smile to predict
    ######
    smiles = list()
    smiles.append(string_smile) #'C1=C(N=CS1)CC(C(=O)O)N'
    df = pd.DataFrame(smiles, columns=['SMILES'])
    # generate fingeprint
    fingerprint, bit_infos = get_fingerprint(df)
    MACC = extract_MACCSds(df['SMILES'])
    #st.write(MACC.shape, fingerprint.shape)
    MACC = MACC[:,0:167]
    #st.write(MACC.shape)
    fingerprint_MACC = np.concatenate((fingerprint, MACC), axis=1)
    #st.write(fingerprint_MACC)
    fingerprint_MACC = np.array(fingerprint_MACC)
    #st.write(fingerprint_MACC.shape)
    ########
    prefix = 'check10pointw'  # xenobiotic1010 model10
    dir_models = 'models_xeno'
    prefix = dir_models + '/' + prefix

    list_xenobiotic = list()
    list_nonXenobiotic = list()
    list_prob = list()
    #st.write('Predicting...')
    for f in range(10):
        #st.write(prefix + str(f + 1) + "_.h5")
        model.load_weights(prefix + str(f + 1) + "_.h5")
        model.compile(optimizer=keras.optimizers.Adam(lr=LR),
                      loss='sparse_categorical_crossentropy',
                      metrics=['sparse_categorical_accuracy'])
        # print("Model loaded==>"+str(f+1))

        y_pred = model.predict(fingerprint_MACC)
        #st.write(y_pred)
        y_score = np.array(y_pred)[:, 1]
        #st.write(y_score)
        #st.write(str(f + 1), '===>', y_pred)
        list_prob.append(y_pred[0])
        if y_pred[0][1] >= 0.50:
            list_xenobiotic.append(y_pred[0][1])
        else:
            list_nonXenobiotic.append(y_pred[0][0])
    #st.write(len(list_nonXenobiotic))

    if len(list_xenobiotic) >= len(list_nonXenobiotic):
        formatted_string = "{:.2%}".format(mean(list_xenobiotic))
        yes_xeno = 'Xenobiotic'
        if batch == 0:
            st.write(smiles[0], '\nis Xenobiotics with a probability of ==>', formatted_string, '\nThis is the decision of ',
                  len(list_xenobiotic), ' out of 10 classifiers')
        #st_radial('You selected', mean(list_xenobiotic))

    else:
        formatted_string = "{:.2%}".format(mean(list_nonXenobiotic))
        yes_xeno = 'Non-Xenobiotic'
        if batch == 0:
            st.write(smiles[0], '\nNON-Xenobiotics with a probability of ==>',formatted_string,
              '\nThis is the Decision of ', len(list_nonXenobiotic), ' out of 10 classifiers')

    ddf_prob = pd.DataFrame(list_prob, columns=['Non-Xenobiotic', 'Xenobiotic'])
    df_prob = ddf_prob.style.format({
        'Non-Xenobiotic': '{:,.2%}'.format,
        'Xenobiotic': '{:,.2%}'.format})
    return df_prob, ddf_prob, yes_xeno, formatted_string

def convert_df(df):
   return df.to_csv().encode('utf-8')
def isNaN(string):
    return string != string
