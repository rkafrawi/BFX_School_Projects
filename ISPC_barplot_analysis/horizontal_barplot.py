#!/usr/bin/env python
import cgitb
import cgi
import plotly.express as px
import plotly.io as pio
import pymysql
import pandas as pd
import plotly.graph_objects as go

cgitb.enable()

# connect to database
connection = pymysql.connect(
    host="localhost",
    user="root",
    password="Irsyadian0!BU2023",
    database="Team_J"
)

# define SQL query
cell_q = '1'
pw_q = 'kegg%'
ind_q = '1'


ind_1_pw_Kg_ct_Mc = """select distinct pathway, pval, padj, ES, NES, nMoreExtreme, ed.size, leadingEdge
                  from expression_data ed join individual i using (cell_id) join cell_type ct using (cell_id)
                  where pathway like '{pw_q}'
                  and cell_id = '{cell_q}'
                  and ed.indiv_id = '{ind_q}'
                  and padj <=0.25;""".format(pw_q=pw_q, cell_q=cell_q, ind_q=ind_q)

# fetch data
result = pd.read_sql_query(ind_1_pw_Kg_ct_Mc, connection)
result['NES'] = pd.to_numeric(result['NES'])
print(result.head(10))
print(len(result))
# if len(result == 0):
#     print("No rows returned")
# else:
duplicates=result['pathway'].value_counts()
print(duplicates)
# print("cocka")

# extra filtration
size_90 = result['size'].quantile(0.9)

# imaginary vars for graph
pathway = 'kegg'
ind = 'ind 2'
cell = 'microglia'

#use plotly to visualize pathways
fig = px.bar(data_frame= result,
            x = 'NES',
            y='pathway', 
            orientation = 'h',
            custom_data=['pval','leadingEdge']
            )

fig.update_layout(xaxis=dict(range=[min(result['NES']) - 0.1, max(result['NES']) + 0.1]))
fig.update_layout(
                xaxis_title='NES', 
                yaxis_title='Pathway',
                title='FGSEA Results for {cell} in {ind} ({pathway})'.format(cell=cell,ind=ind,pathway=pathway))

fig.update_xaxes(tickfont=dict(size=8))
fig.update_yaxes(tickfont=dict(size=8), showticklabels=True)

#check plotly default hovertemplate
print("plotly express hovertemplate:", fig.data[0].hovertemplate)

fig.update_traces(marker_line_width=0)

fig.update_traces(hovertemplate='P-val: %{customdata[0]}<br>Leading Edge Genes: %{customdata[1]}')

print("usr def hovertemp:", fig.data[0].hovertemplate)
fig.show()

# Define the file name based on the input variables
file_name = "C:/Users/rkafr/Desktop/BF768 Project/Figures/fgsea_results_{cell}_{ind}_{pathway}.jpeg".format(cell=cell,ind=ind,pathway=pathway)

# Save the figure as a JPEG image with the dynamic file name
pio.write_image(fig, file_name)

# close connection
connection.close()
