# import
from ncnetworkviz import Jaal
from ncnetworkviz.datasets import load_got
import pandas as pd

# load the data
edge_df = pd.read_csv('ectoderm_grn.csv')

# define vis options
vis_opts = {'height': '600px', # change height
            'interaction':{'hover': True}, # turn on-off the hover 
            'physics':{'stabilization':{'iterations': 100}}} # define the convergence iteration of network

# init Jaal and run server (with opts)
#Jaal(edge_df, node_df).plot(host="neucrest.herokuapp.com", port="8050", vis_opts=vis_opts)

# init Jaal and run server (with gunicorn)
app = Jaal(edge_df).create()
server = app.server
