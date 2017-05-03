# coding: utf-8
get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
import stat_it as si
it = si.it_method
dct = si.io.loadmat('../mat/ddom_20x20x24_210.mat')
A, b = [dct[key] for key in ['A','b']]
