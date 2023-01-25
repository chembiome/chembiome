import streamlit as st
import hydralit_components as hc
import tensorflow.compat.v1 as tf
#######################
# https://mcsorkun-chemplot-web-web-app-chemplot-jrrecy.streamlitapp.com/
# https://github.com/Mariewelt/OpenChem/tree/master/devtools
# https://ai-dd.eu/sites/default/files/school-1/esben.pdf
################
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
import datetime
from streamlit_option_menu import option_menu
import streamlit.components.v1 as components

import sys, os

import pandas as pd

#from functions import predict_xeno
#from functions import isNaN
from menu_predict import perform_toxicity_prediction
from menu_predict import perform_xenobiotic_prediction
from menu_chembiome import core_microbiome_organisms
from menu_chembiome import core_microbiome_reactions
from menu_chembiome import core_microbiome_pathways
from menu_chembiome import core_microbiome_enzymes
#C:\Users\masun\Anaconda3\envs\enzyme_predict\chembiome

#make it look nice from the start
st.set_page_config(layout='wide',initial_sidebar_state='collapsed',)

# specify the primary menu definition
menu_data = [
    {'icon': "far fa-copy", 'label':"About"},
    #{'id':'Copy','icon':"🐙",'label':"Copy"},
    {'icon': "fa-solid fa-radar",'label':"CHEMBIOME", 'submenu':[{'id':' Organisms','icon': "fa fa-database", 'label':"Organisms"},{'id':'Reactions','icon': "fa fa-database", 'label':"Reactions"},{'id':'Pathways','icon': "fa fa-database", 'label':"Pathways"}, {'id':'enzymes','icon': "fa fa-database", 'label':"Enzymes"}, {'id':'metabolites','icon': "fa fa-database", 'label':"Metabolites"}]},
    {'icon': "fa-solid fa-radar",'label':"METAGENOMES", 'submenu':[{'id':' subid21','icon': "🧬️", 'label':"Diseases"},{'id':'subid22','icon': "🧬️", 'label':"Studies"},{'id':'subid23','icon': "🧬️", 'label':"Pathways"}, {'id':'subid13','icon': "🧬️", 'label':"Genes"}]},
    {'icon': "fa-solid fa-radar",'label':"PREDICT", 'submenu':[{'id':' subid21','icon': "❁", 'label':"Enzymes and Organisms"},{'id':'xeno','icon': "❁", 'label':"Xenobiotics or Non-xenobiotics"},{'id':'tox','icon': "❁", 'label':"Toxic and Non-toxic"}, {'id':'subid13','icon': "❁", 'label':"Xenobiotic Biotransformation"}]},
    {'icon': "fa-solid fa-radar",'label':"SEARCH & DOWNLOAD", 'submenu':[{'id':' subid21','icon': "🔎", 'label':"CHEMBIOME"},{'id':'subid22','icon': "🔎", 'label':"FoodDB"},{'id':'subid23','icon': "🔎", 'label':"FlavorDB"}, {'id':'subid13','icon': "🔎", 'label':"ToxicDB"}]},
    {'icon': "fa-thin fa-diagram-project",'label':"VISUALIZE", 'submenu':[{'id':' subid21','icon': "⌘", 'label':"Chemical compounds"},{'id':'subid22','icon': "far fa-chart-bar", 'label':"Visualize1"},{'id':'subid23','icon': "far fa-chart-bar", 'label':"Visualiz2"}, {'id':'subid13','icon': "far fa-chart-bar", 'label':"Visualiz3"}]},
    {'icon': "✉", 'label':"Contact"},#no tooltip message
    #{'id':' Crazy return value 💀','icon': "💀", 'label':"Calendar"},
    #{'icon': "fas fa-tachometer-alt", 'label':"Dashboard",'ttip':"I'm the Dashboard tooltip!"}, #can add a tooltip message
    #{'icon': "far fa-copy", 'label':"Right End"},
    #{'icon': "fa-solid fa-radar",'label':"Dropdown2", 'submenu':[{'label':"Sub-item 1", 'icon': "fa fa-meh"},{'label':"Sub-item 2"},{'icon':'🙉','label':"Sub-item 3",}]},
]
over_theme = {'txc_inactive': 'white','menu_background':'purple','txc_active':'purple','option_active':'white'}
#over_theme = {'txc_inactive': '#FFFFFF'}
menu_id = hc.nav_bar(
    menu_definition=menu_data,
    override_theme=over_theme,
    home_name='Home',
    #login_name='Logout',
    hide_streamlit_markers=True, #will show the st hamburger as well as the navbar now!
    sticky_nav=False, #at the top or not
    sticky_mode='sticky', #jumpy or not-jumpy, but sticky or pinned
)

#get the id of the menu item clicked
#st.info(f"{menu_id}")


#with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
if menu_id.strip() == 'tox':
    perform_toxicity_prediction()
elif menu_id.strip() == 'xeno':
    perform_xenobiotic_prediction()
elif menu_id.strip()== 'Organisms':
    core_microbiome_organisms()
elif menu_id.strip() == 'Reactions':
    core_microbiome_reactions()
elif menu_id.strip() == 'Pathways':
    core_microbiome_pathways()
elif menu_id.strip() == 'enzymes':
    core_microbiome_enzymes()

        #HtmlFile = open("homes.html")
        #source_code_2 = HtmlFile.read()
        #components.html(HtmlFile.read())




#if st.button('click me'):
    #st.info('You clicked at: {}'.format(datetime.datetime.now()))


    #components.html(source_code_2, height=700)
#if st.sidebar.button('click me too'):
  #st.info('You clicked at: {}'.format(datetime.datetime.now()))

