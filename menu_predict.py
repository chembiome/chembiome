import streamlit as st
from functions import predict_tox
from functions import predict_xeno
from functions import isNaN
import pandas as pd

import sys
import os
from prepare_input_file import canonicalise_smile
from prepare_input_file import smi_tokenizer
#sys.path.append('./MetaTrans')
#from MetaTrans.translate import *
import prepare_input_file

def perform_toxicity_prediction():
    st.info("Toxic or Non-Toxic")
    title = st.text_input('SMILE', 'C1=C(N=CS1)CC(C(=O)O)N')
    if st.button('Predict'):
        df_prob = predict_tox(title.strip())
        st.write('Probabilities of the classifier:')
        st.write(df_prob)

def perform_xenobiotic_prediction():
    st.header("Xenobiotics vs Non-xenobiotics")
    choice = st.radio("Predict", ('SMILE', 'Upload CSV'))
    if choice == "SMILE":
        st.subheader("SMILE")
        title = st.text_input('SMILE', 'C1=C(N=CS1)CC(C(=O)O)N')
        if st.button('Predict'):
            df_prob, ddf_prob, yes_xeno, formatted_string = predict_xeno(title.strip(), 0)
            st.write('Probabilities of 10 classifiers')
            st.write(df_prob)
    else:
        st.subheader("File Upload :")
        data_file = st.file_uploader(
            "Please upload a CSV with a column that contains the smiles structures. The column title should be (SMILE)",
            type=["csv"])
        st.write(data_file)
        if data_file is not None:
            file_details = {"filename": data_file.name, "filetype": data_file.type,
                            "filesize": data_file.size}

            # st.write(file_details)
            df = pd.read_csv(data_file, sep='\t')
            st.dataframe(df)
            if 'SMILE' not in df.columns:
                st.write('Please be sure that your file has a column with the name SMILE (all capital letters)')
            else:
                if st.button('Predict'):
                    with st.spinner('Predicting, please wait for it...'):

                        my_smiles = df['SMILE'].tolist()
                        ok = 1
                        for sm in my_smiles:
                            if isNaN(sm):
                                ok = 0
                                break
                        if ok == 0:
                            st.write('We cannot process your file, there is a compound without SMILE structure')
                        elif ok == 1:
                            list_yes_xeno, list_formatted_string = list(), list()
                            df_prob, ddf_prob, yes_xeno, formatted_string = predict_xeno(my_smiles[0].strip(), 1)
                            list_yes_xeno.append(yes_xeno)
                            list_formatted_string.append(formatted_string)
                            for it in range(1, len(my_smiles)):
                                my_smile = my_smiles[it].strip()
                                df_prob, ddf_prob, yes_xeno, formatted_string = predict_xeno(my_smiles[it].strip(), 1)
                                list_yes_xeno.append(yes_xeno)
                                list_formatted_string.append(formatted_string)

                                # break
                            df_yes_xeno = pd.DataFrame(list_yes_xeno, columns=['Status'])
                            df_formatted_string = pd.DataFrame(list_formatted_string, columns=['Probability'])
                            df_xeno_prediction = pd.concat([df_yes_xeno, df], axis=1)
                            df_xeno_prediction = pd.concat([df_formatted_string, df_xeno_prediction], axis=1)
                            st.write(df_xeno_prediction)
                            df_xeno_prediction = df_xeno_prediction.to_csv(index=False).encode('utf-8')
                            st.download_button(
                                label="Press to Download",
                                data=df_xeno_prediction,
                                file_name='prediction_xeno.csv',
                                mime="text/csv"

                            )
                            # print ('hi')

def perform_biotransformation_prediction():
    st.header("Biotransformation Prediction")
    smile = st.text_input('SMILE', 'C1=C(N=CS1)CC(C(=O)O)N')
    if st.button('Predict'):
        #smiles = canonicalise_smile(smile)
        #smiles_tok = smi_tokenizer(smile)
        df = pd.DataFrame(smile)
        df.to_csv('smile_biotransformation.csv', index=True)
        prepare_input_file.main ()


        df_prob, ddf_prob, yes_xeno, formatted_string = predict_xeno(title.strip(), 0)
        st.write('Probabilities of 10 classifiers')
        st.write(df_prob)