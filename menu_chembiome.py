import streamlit as st
import psycopg2
from st_aggrid import AgGrid,GridUpdateMode
from st_aggrid.grid_options_builder import GridOptionsBuilder
import psycopg2.extras
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Draw import IPythonConsole
from PIL import Image
#from functions import string_2_html

########################### Connect to the database #################
hostname = 'localhost'
database = 'chembiome'
username = 'postgres'
pwd = '1'
port_id = 5432
conn = None
################################# Functions
def sql_executor(raw_code, val):
    #with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
    conn = psycopg2.connect(host=hostname, port=port_id, dbname=database, user=username,
                            password=pwd)
    cur = conn.cursor()
    cur.execute(raw_code, (val,))
    data = cur.fetchall()
    cur.close()
    return data

def string_2_html(my_string, field):
    if '-' in my_string:
        return '-'
    my_string = my_string.replace('__ ', '__')
    list_op = my_string.split('__')
    html_text = ''
    name = ''
    for op in list_op:
        #st.write(op)
        ko = 0
        if '-' in op:
            continue
        else:
            if field == 'enzymes':
                query_result = sql_executor(
                    'SELECT * FROM "kegg_enzymes" where "EC_Number" =%s', op)
                #st.write(query_result)
                if query_result != None  :
                    name = query_result[0][1].split(';')[0]
            elif field == 'pathways':
                if 'rn' in op:
                    op = op.replace('rn','map')
                elif 'ko' in op :
                    op = op.replace('ko', 'map')
                    ko=1
                query_result = sql_executor(
                    'SELECT * FROM "kegg_pathways" where "pathway_id" =%s', op)
                #st.write(query_result)
                if query_result != None:
                    name = query_result[0][1].split(';')[0]
                    if ko==1:
                        op = op.replace('map', 'ko')
            elif field == 'modules':
                query_result = sql_executor(
                    'SELECT * FROM "kegg_modules" where "module_id" =%s', op)
                if query_result != None:
                    name = query_result[0][1]
            elif field == 'orthology':
                query_result = sql_executor(
                    'SELECT * FROM "kegg_orthology" where "orthology_id" =%s', op)
                if query_result != None:
                    name = query_result[0][1].split(';')[1]
                    #name = name.split('[')[0]
            elif field == 'Rclass':
                op = op.split(' ')[0]
                query_result = sql_executor(
                    'SELECT * FROM "kegg_Rclass" where "rclass_id" =%s', op)
                if query_result != None:
                    name = query_result[0][1]
            elif field == 'disease':
                op = op.split(' ')[0]
                query_result = sql_executor(
                    'SELECT * FROM "kegg_diseases" where "disease_id" =%s', op)
                if query_result != None:
                    name = query_result[0][1]
            html_text = html_text + "<a href='https://www.kegg.jp/entry/"+op+"'>"+op+"</a>\t"+name+"</br>"
    #st.write(query_result[0])
    return html_text

def core_microbiome_organisms ():
    try:
        with psycopg2.connect(
                host=hostname,
                dbname=database,
                user=username,
                password=pwd,
                port=port_id) as conn:
            #st.header("Core Gut Microbiome")
            # st.info(f"{menu_id.strip()}")
            # cur.execute('Select * FROM metadata_genomes, core_genome WHERE core_genome.genome_name = metadata_genomes.Genome')
            # for record in cur.fetchall():
            # st.write(record['Domain'], record['Phylum'])
            # raw_code = 'Select * FROM metadata_genomes, core_genome WHERE core_genome.genome_name = metadata_genomes.Genome'
            # raw_code = 'Select * FROM metadata_genomes'
            # query_results = sql_executor(raw_code)
            # df = pd.DataFrame(query_results)
            # st.write(df.head())
            st.header("Core Gut Microbiome's Organisms")
            with st.spinner('Searching, please wait for it...'):
                query_results = sql_executor(
                    'SELECT * FROM "core_genome"  LEFT JOIN  "metadata_genomes" on "genome_name" = "Genome";',
                    None)
                df = pd.DataFrame(query_results)
                #st.write(df.head(5))
                #df_metadata_genomes = pd.read_sql_query('Select * FROM metadata_genomes', con=conn)

                #df_core_genome = pd.read_sql_query('Select * FROM core_genome  ', con=conn)
                #df = pd.merge(df_core_genome, df_metadata_genomes, left_on=['genome_name'], right_on=['Genome'], how='left')
                df.rename(columns={df.columns[0]: 'Genome', df.columns[21]: 'Domain', df.columns[22]: 'Phylum',
                                   df.columns[23]: 'Class', df.columns[24]: 'Ordr', df.columns[25]: 'Family',
                                   df.columns[26]: 'Genus', df.columns[27]: 'Specie', df.columns[18]: 'Country',
                                   df.columns[19]: 'Continent'},
                          inplace=True)
                df = df[['Genome', 'Domain', 'Phylum', 'Class', 'Ordr', 'Family', 'Genus', 'Specie', 'Country','Continent']]
                # df = df.astype(str)
                # st.dataframe(df)
                st.cache()
                gd = GridOptionsBuilder.from_dataframe(df)
                gd.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
                gd.configure_default_column(editable=False, groupable=True)
                gd.configure_pagination(enabled=True)
                gridoptions = gd.build()
                grid1 = AgGrid(df, gridOptions=gridoptions,
                               update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True)

                sel_row = grid1["selected_rows"]  # Type -> List
                # st.write(sel_row)

                df_sel = pd.DataFrame(sel_row)
                if len(df_sel) > 0:
                    st.write(df_sel.iat[0,0])
                    query_sel = sql_executor(
                        'SELECT * FROM "metadata_genomes" where "Species_rep" =%s', df_sel.iat[0, 0])
                    df_sel1 = pd.DataFrame(query_sel)
                    if len(df_sel1)>0 :
                        #st.write(df_sel)
                        df_sel1.rename(columns={df_sel1.columns[0]: 'Genome', df_sel1.columns[20]: 'Domain', df_sel1.columns[21]: 'Phylum',
                                           df_sel1.columns[22]: 'Class', df_sel1.columns[23]: 'Ordr', df_sel1.columns[24]: 'Family',
                                           df_sel1.columns[25]: 'Genus', df_sel1.columns[26]: 'Specie', df_sel1.columns[1]: 'Genome_type', df_sel1.columns[17]: 'Country',
                                           df_sel1.columns[18]: 'Continent', df_sel1.columns[15]: 'Sample_accession', df_sel1.columns[12]: 'Genome_accession'},
                                  inplace=True)

                        #df_sel = df_metadata_genomes[df_metadata_genomes.Species_rep == df_sel.iat[0, 0]]  # Convert list to dataframe
                        df_sel1 = df_sel1[
                            ['Genome', 'Genus', 'Specie', 'Genome_type', 'Country', 'Continent', 'Genome_accession', 'Sample_accession']]
                        df_sel1 = df_sel1[df_sel1['Genome'] != df_sel.iat[0, 0]]
                        st.subheader("Extended Genomes")
                        # Grid 2 - Highlight
                        gd2 = GridOptionsBuilder.from_dataframe(df_sel1)
                        gd2.configure_selection(selection_mode='single', use_checkbox=False)  # multiple
                        gd2.configure_default_column(editable=False, groupable=True)
                        gd2.configure_pagination(enabled=True)
                        gridoptions2 = gd2.build()
                        grid2 = AgGrid(df_sel1, gridOptions=gridoptions2,
                                       update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True, height=250)
                        # sel_row2 = grid2["selected_rows"]  # Type -> List
                        # df_sel2 = pd.DataFrame(sel_row2)
                        if len(df_sel) > 0:
                            # st.write(df_sel.iat[0,0])
                            df_metadata_genome_detailes = pd.read_sql_query('Select * FROM metadata_genome_detailes', con=conn)
                            df_sel2 = df_metadata_genome_detailes[
                                df_metadata_genome_detailes.genome_detail_id.str.contains(df_sel.iat[0, 0])]
                            # df_sel2 = df_metadata_genome_detailes.str.contains(df_sel.iat[0, 0], regex=False)
                            st.subheader("Genomes Functional Annotation")
                            # Grid 2 - Highlight
                            gd3 = GridOptionsBuilder.from_dataframe(df_sel2)
                            gd3.configure_selection(selection_mode='single', use_checkbox=False)  # multiple
                            gd3.configure_default_column(editable=False, groupable=True)
                            gd3.configure_pagination(enabled=True)
                            gridoptions3 = gd3.build()
                            grid3 = AgGrid(df_sel2, gridOptions=gridoptions3,
                                           update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True)
    except Exception as error:
        st.write(error)
    finally:
        if conn is not None:
            conn.close()

def core_microbiome_reactions():
    try:
         with psycopg2.connect(
                 host=hostname,
                 dbname=database,
                 user=username,
                 password=pwd,
                 port=port_id) as conn:
            st.header("Core Gut Microbiome's Reactions")
            query_results = sql_executor('SELECT * FROM "chembiome_unique_reactions" LEFT JOIN "kegg_reactions_catalogue" on "KEGG_Reaction" = "reaction_ID";', None)
            df = pd.DataFrame(query_results)
            df.fillna("-", inplace=True)
            df.rename(columns={df.columns[0]: 'KEGG Reaction', df.columns[2]: 'Name', df.columns[3]: 'Definition', df.columns[4]: 'Equation', df.columns[5]: 'Enzyme', df.columns[6]: 'RClass', df.columns[7]: 'Pathways', df.columns[8]: 'Modules', df.columns[9]: 'Orthology'}, inplace=True)
            #st.write(df.columns)
            df_reactions = df[['KEGG Reaction', 'Name', 'Definition', 'Equation']]
            gd2 = GridOptionsBuilder.from_dataframe(df_reactions)
            gd2.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
            gd2.configure_default_column(editable=False, groupable=False)
            gd2.configure_pagination(enabled=True)
            gridoptions2 = gd2.build()
            grid2 = AgGrid(df_reactions, gridOptions=gridoptions2,
                           update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True, height=250 )

            sel_row = grid2["selected_rows"]  # Type -> List
            if len(sel_row) > 0:
                #mol = Chem.MolFromSmiles('c1ccncc1')
                rxn = Reactions.ReactionFromSmarts('CC(=O)C>>CC(O)C', useSmiles=True)

                #Draw.MolToImage(mol)
                drawer = rdMolDraw2D.MolDraw2DCairo(800, 200)
                drawer.SetFontSize(1.0)
                drawer.DrawReaction(rxn)
                drawer.FinishDrawing()
                drawer.WriteDrawingText('test2.png')
                im = Image.open('test2.png')

                #st.write(sel_row)
                #st.write(sel_row[0])
                keggReaction = sel_row[0]["KEGG Reaction"]
                #st.write(keggReaction)
                query_genome_reactions = sql_executor(
                    'SELECT "RGenome", "Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "count_sequences" AS "sequences in the Genome", count(*) as "Reaction Abundance"  FROM "chembiome_Genomes_Reactions" LEFT OUTER JOIN "metadata_genomes"  ON "RGenome"="Genome" LEFT OUTER JOIN "chembiome_genome_count_sequences" ON "RGenome" = "genomeid" where "KEGG_Reaction"= %s group by "RGenome","Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "sequences in the Genome" order by "Reaction Abundance" DESC ', keggReaction)
                df_query_genome_reactions = pd.DataFrame(query_genome_reactions)
                if len(df_query_genome_reactions)>0:
                    df_query_genome_reactions.rename(columns={df_query_genome_reactions.columns[0]: 'Genome ID', df_query_genome_reactions.columns[1]: 'Domain',
                                                df_query_genome_reactions.columns[2]: 'Phylum', df_query_genome_reactions.columns[3]: 'Class',
                                                df_query_genome_reactions.columns[4]: 'Ordr', df_query_genome_reactions.columns[5]: 'Family',
                                                df_query_genome_reactions.columns[6]: 'Genus', df_query_genome_reactions.columns[7]: 'Specie',
                                                df_query_genome_reactions.columns[8]: 'Total Sequences',
                                                df_query_genome_reactions.columns[9]: 'Reaction Qty'},
                                    inplace=True)
                    #df_query_genome_reactions = df_query_genome_reactions[['Genome ID', 'Domain', 'Phylum', 'Class', 'Ordr', 'Family', 'Genus', 'Specie']]
                    #query_genome_count = sql_executor('SELECT * FROM "chembiome_genome_count_sequences" ', '')
                    #df_query_genome_count = pd.DataFrame(query_genome_count, columns=['Genome_ID', 'countt'])
                    #df_genome_reactions_count = pd.merge(df_query_genome_reactions, df_query_genome_count,
                                                               #left_on=['Genome ID'],
                                                               #right_on=['Genome_ID'], how='left')
                    #st.write(df_query_genome_reactions)
                    #df_query_genome_reactions ['count']
                    #st.write(df_query_genome_reactions)
                    df_query_genome_reactions['Reaction Abundance'] = df_query_genome_reactions['Reaction Qty'] / df_query_genome_reactions['Total Sequences']
                    gd3 = GridOptionsBuilder.from_dataframe(df_query_genome_reactions)
                    gd3.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
                    gd3.configure_default_column(editable=False, groupable=False)
                    gd3.configure_pagination(enabled=True)
                    gridoptions = gd3.build()
                    grid2_gdgr = AgGrid(df_query_genome_reactions, gridOptions=gridoptions,
                                   update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True, height=250)
                    sel_row_gr = grid2_gdgr["selected_rows"]


                ########## KEGG Metadata ##########
                st.header("KEGG Metadata")
                query_links = sql_executor(
                    'SELECT * FROM "chembiome_reactions_links" where kegg_reaction =%s', keggReaction)
                df_query_links = pd.DataFrame(query_links)

                df_links = pd.DataFrame(df_query_links)
                links = list()
                for c in df_links.columns:
                    if not df_links.iloc[0][c]:
                        links.append('')
                    else:
                        links.append(df_links.iloc[0][c])
                #st.write(links)
                df_links.rename(columns={df_links.columns[0]: 'KEGG Reaction', df_links.columns[1]: 'Metanetx ID', df_links.columns[2]: 'MetaCyc ID', df_links.columns[3]: 'In Pathways'}, inplace=True)
                #st.write(df_links)
                #st.write(df_links.columns)
                df_more_details = df [df['KEGG Reaction']==links[0]]
                df_more_details.reset_index(drop=True)
                #st.write(df_more_details.iloc[0]["Pathways"])
                if sel_row[0]["Name"] == '-':
                    st.write('There is no more data regarding this reaction:\t'+links[0])
                else:
                    html_string = "<table>" \
                                    "<tr><td>KEGG reaction ID : </td><td><a href='https://www.kegg.jp/entry/"+links[0]+"'>"+links[0]+"</a></td></tr>" \
                                    "<tr><td>Name : </td><td>" + sel_row[0]["Name"] + "</td></tr>" \
                                    "<tr><td>Equation : </td><td>" + df_more_details.iloc[0]["Equation"] + "</td></tr>" \
                                    "<tr><td>Enzymes : </td><td>" + string_2_html (df_more_details.iloc[0]["Enzyme"], 'enzymes') + "</td></tr>" \
                                    "<tr><td>Pathways : </td><td>" + string_2_html (df_more_details.iloc[0]["Pathways"], 'pathways') + "</td></tr>" \
                                    "<tr><td>Modules:</td><td>" + string_2_html (df_more_details.iloc[0]["Modules"], 'modules') + "</td></tr>" \
                                    "<tr><td>Orthology:</td><td>" + string_2_html (df_more_details.iloc[0]["Orthology"], 'orthology') + "</td></tr>" \
                                    "<tr><td>RClass:</td><td>" + string_2_html (df_more_details.iloc[0]["RClass"], 'Rclass') + "</td></tr>" \
                                     "<tr><td>Metanetx ID:</td><td><a href='https://www.metanetx.org/equa_info/"+links[1]+"'>"+links[1]+"</a></td></tr>" \
                                    "<td>MetaCyc ID:</td><td><a href='https://biocyc.org/META/NEW-IMAGE?type=REACTION&object="+links[2]+"'>"+links[2]+"</a></td></tr>" \
                                    "</table)"
                    st.markdown(html_string, unsafe_allow_html=True)
                    #if len(df_sel) > 0:

    except Exception as error:
        st.write(error)
    finally:
        if conn is not None:
            conn.close()

def core_microbiome_pathways():
    try:
         with psycopg2.connect(
                 host=hostname,
                 dbname=database,
                 user=username,
                 password=pwd,
                 port=port_id) as conn:
            st.header("Core Gut Microbiome's pathways")
            query_results = sql_executor('SELECT * FROM "chembiome_unique_pathways" LEFT JOIN "kegg_pathways_ko_catalogue" on "KEGG_map_code" = "pathway_ID";', None)
            df = pd.DataFrame(query_results)
            df.fillna("-", inplace=True)
            df.rename(columns={df.columns[1]: 'KEGG Pathway ID', df.columns[2]: 'KEGG Ko Code', df.columns[4]: 'Name', df.columns[5]: 'Class', df.columns[6]: 'PATHWAY_MAP', df.columns[7]: 'DISEASE', df.columns[8]: 'DRUG', df.columns[9]: 'Modules', df.columns[10]: 'DBLINKS', df.columns[11]: 'Orthology', df.columns[12]: 'Compounds', df.columns[13]: 'Related Pathways'}, inplace=True)
            #st.write(df)
            df_pathway = df[['KEGG Pathway ID', 'Name', 'Class']]
            gd2 = GridOptionsBuilder.from_dataframe(df_pathway)
            gd2.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
            gd2.configure_default_column(editable=False, groupable=False)
            gd2.configure_pagination(enabled=True)
            gridoptions2 = gd2.build()
            grid2 = AgGrid(df_pathway, gridOptions=gridoptions2,
                           update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True, height=250 )

            sel_row = grid2["selected_rows"]  # Type -> List
            if len(sel_row) > 0:
                keggPath = sel_row[0]["KEGG Pathway ID"]
                # st.write(keggReaction)
                query_genome_Path = sql_executor(
                    'SELECT "PGenome", "Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "count_sequences" AS "sequences in the Genome", count(*) as "Pathway Abundance"  FROM "chembiome_Genomes_Pathways" LEFT OUTER JOIN "metadata_genomes"  ON "PGenome"="Genome" LEFT OUTER JOIN "chembiome_genome_count_sequences" ON "PGenome" = "genomeid" where "KEGG_Pathway"= %s group by "PGenome","Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "sequences in the Genome" order by "Pathway Abundance" DESC ',
                    keggPath)
                df_query_genome_Path = pd.DataFrame(query_genome_Path)
                st.subheader('Organisms')
                if len(df_query_genome_Path) > 0:
                    df_query_genome_Path.rename(columns={df_query_genome_Path.columns[0]: 'Genome ID',
                                                              df_query_genome_Path.columns[1]: 'Domain',
                                                              df_query_genome_Path.columns[2]: 'Phylum',
                                                              df_query_genome_Path.columns[3]: 'Class',
                                                              df_query_genome_Path.columns[4]: 'Ordr',
                                                              df_query_genome_Path.columns[5]: 'Family',
                                                              df_query_genome_Path.columns[6]: 'Genus',
                                                              df_query_genome_Path.columns[7]: 'Specie',
                                                              df_query_genome_Path.columns[8]: 'Total Sequences',
                                                              df_query_genome_Path.columns[9]: 'Pathway Qty'},
                                                     inplace=True)

                    df_query_genome_Path['Pathway Abundance'] = df_query_genome_Path['Pathway Qty'] / \
                                                                      df_query_genome_Path['Total Sequences']
                    gd3 = GridOptionsBuilder.from_dataframe(df_query_genome_Path)
                    gd3.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
                    gd3.configure_default_column(editable=False, groupable=False)
                    gd3.configure_pagination(enabled=True)
                    gridoptions = gd3.build()
                    grid2_gdgr = AgGrid(df_query_genome_Path, gridOptions=gridoptions,
                                        update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True,
                                        height=250)
                    sel_row_gr = grid2_gdgr["selected_rows"]


                st.subheader('KEGG Metadata')
                df_more_details = df [df['KEGG Pathway ID']==sel_row[0]["KEGG Pathway ID"]]
                df_more_details.reset_index(drop=True)
                st.write(df_more_details.iloc[0]["KEGG Pathway ID"])
                ID_map = df_more_details.iloc[0]["KEGG Pathway ID"]
                ID_ko = df_more_details.iloc[0]["KEGG Ko Code"]
                if df_more_details.iloc[0]["Name"] == '-':
                    st.write('No data regarding this pathway :\t'+ID_map)
                else:
                    with st.spinner('Searching, please wait for it...'):
                        html_string = "<table>" \
                                        "<tr><td>KEGG pathway ID : </td><td><a href='https://www.kegg.jp/entry/"+ID_map+"'>"+ID_map+"</a>&emsp;<a href='https://www.kegg.jp/entry/"+ID_ko+"'>"+ID_ko+"</a></td></tr>" \
                                        "<tr><td>Name : </td><td>" + sel_row[0]["Name"] + "</td></tr>" \
                                        "<tr><td>Class : </td><td>" + df_more_details.iloc[0]["Class"] + "</td></tr>" \
                                        "<tr><td>PATHWAY_MAP : </td><td>" + string_2_html (df_more_details.iloc[0]["PATHWAY_MAP"], 'pathways') + "</td></tr>" \
                                        "<tr><td>Modules :</td><td>" + string_2_html (df_more_details.iloc[0]["Modules"], 'modules') + "</td></tr>" \
                                        "<tr><td>Disease :</td><td>" + string_2_html(df_more_details.iloc[0]["DISEASE"], 'disease') + "</td></tr>" \
                                        "<tr><td>Orthology:</td><td>" + string_2_html (df_more_details.iloc[0]["Orthology"], 'orthology') + "</td></tr>" \
                                        "<tr><td>Orthology:</td><td>" + string_2_html(df_more_details.iloc[0]["Related Pathways"], 'pathways') + "</td></tr>" \
                                         "</table)"
                        st.markdown(html_string, unsafe_allow_html=True)
                #if len(df_sel) > 0:

    except Exception as error:
        st.write(error)
    finally:
        if conn is not None:
            conn.close()

def core_microbiome_enzymes():
    try:
         with psycopg2.connect(
                 host=hostname,
                 dbname=database,
                 user=username,
                 password=pwd,
                 port=port_id) as conn:
            st.header("Core Gut Microbiome's Enzymes")
            query_results = sql_executor('SELECT * FROM "chembiome_unique_enzymes" LEFT JOIN "kegg_enzymes" on "KEGG_EC" = "EC_Number";', None)
            df = pd.DataFrame(query_results)
            df.fillna("-", inplace=True)
            df.rename(columns={df.columns[0]: 'EC Number', df.columns[2]: 'Names'}, inplace=True)
            #st.write(df.columns)
            df_ec = df[['EC Number', 'Names']]
            gd2 = GridOptionsBuilder.from_dataframe(df_ec)
            gd2.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
            gd2.configure_default_column(editable=False, groupable=False)
            gd2.configure_pagination(enabled=True)
            gridoptions2 = gd2.build()
            grid2 = AgGrid(df_ec, gridOptions=gridoptions2,
                           update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True, height=250 )

            sel_row = grid2["selected_rows"]
            if len(sel_row) > 0:
                keggPath = sel_row[0]["EC Number"]
                # st.write(keggReaction)
                query_genome_EC = sql_executor(
                    'SELECT "RGenome", "Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "count_sequences" AS "sequences in the Genome", count(*) as "Enzyme Abundance"  FROM "chembiome_Genome_Enzymes" LEFT OUTER JOIN "metadata_genomes"  ON "RGenome"="Genome" LEFT OUTER JOIN "chembiome_genome_count_sequences" ON "RGenome" = "genomeid" where "EC"= %s group by "RGenome","Domain", "Phylum", "Class", "Ordr", "Family", "Genus", "Specie", "sequences in the Genome" order by "Enzyme Abundance" DESC ',
                    keggPath)
                df_query_genome_EC = pd.DataFrame(query_genome_EC)
                st.subheader('Organisms')
                if len(df_query_genome_EC) > 0:
                    df_query_genome_EC.rename(columns={df_query_genome_EC.columns[0]: 'Genome ID',
                                                              df_query_genome_EC.columns[1]: 'Domain',
                                                              df_query_genome_EC.columns[2]: 'Phylum',
                                                              df_query_genome_EC.columns[3]: 'Class',
                                                              df_query_genome_EC.columns[4]: 'Ordr',
                                                              df_query_genome_EC.columns[5]: 'Family',
                                                              df_query_genome_EC.columns[6]: 'Genus',
                                                              df_query_genome_EC.columns[7]: 'Specie',
                                                              df_query_genome_EC.columns[8]: 'Total Sequences',
                                                              df_query_genome_EC.columns[9]: 'EC Qty'},
                                                     inplace=True)

                    df_query_genome_EC['Enzyme Abundance'] = df_query_genome_EC['EC Qty'] / \
                                                                      df_query_genome_EC['Total Sequences']
                    gd3 = GridOptionsBuilder.from_dataframe(df_query_genome_EC)
                    gd3.configure_selection(selection_mode='single', use_checkbox=True)  # multiple
                    gd3.configure_default_column(editable=False, groupable=False)
                    gd3.configure_pagination(enabled=True)
                    gridoptions = gd3.build()
                    grid2_gdgr = AgGrid(df_query_genome_EC, gridOptions=gridoptions,
                                        update_mode=GridUpdateMode.SELECTION_CHANGED, enable_enterprise_modules=True,
                                        height=250)
                    sel_row_gr = grid2_gdgr["selected_rows"]
    except Exception as error:
        st.write(error)
    finally:
        if conn is not None:
            conn.close()