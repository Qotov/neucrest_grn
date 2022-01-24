# import
from jaal import Jaal
from jaal.datasets import load_got
import pandas as pd

# load the data
edge_df = pd.read_csv('viz_dfs/nc_grn.csv')

# define vis options
vis_opts = {'height': '600px', # change height
            'interaction':{'hover': True}, # turn on-off the hover 
            'physics':{'stabilization':{'iterations': 100}}} # define the convergence iteration of network

# init Jaal and run server (with gunicorn)
#app = Jaal(edge_df).create()
#server = app.server

Jaal(edge_df).plot(host="https://neuralcrest.herokuapp.com", vis_opts=vis_opts)
