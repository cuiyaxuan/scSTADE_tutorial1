Tutorial 2: 10X Visium Breast dataset
=====================================

.. raw:: html

    <div style="font-size: 15px;">In this tutorial, we show how to apply scSTADE to identify spatial domains on 10X Visium data. As a example, we analyse the Breast dataset.</div>

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    <div style="font-size: 15px;">The source code package is freely available at https://github.com/cuiyaxuan/scSTADE/tree/master. The datasets used in this study can be found at https://drive.google.com/drive/folders/1H-ymfCqlDR1wpMRX-bCewAjG5nOrIF51?usp=sharing.</div>

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    #Low-resolution spatial transcriptomics data simplified version.
    from scSTADE import scSTADE
    import os
    import torch
    import pandas as pd
    import scanpy as sc
    from sklearn import metrics
    import multiprocessing as mp
    
    def setup_seed(seed=41):
        import torch
        import os
        import numpy as np
        import random
        torch.manual_seed(seed)  
        np.random.seed(seed)  # Numpy module.
        random.seed(seed)  # Python random module.
        if torch.cuda.is_available():
            # torch.backends.cudnn.benchmark = False
            torch.backends.cudnn.deterministic = True
            torch.cuda.manual_seed(seed)  
            torch.cuda.manual_seed_all(seed) 
            #os.environ['PYTHONHASHSEED'] = str(seed)
    
    setup_seed(41)
    
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    
    
    n_clusters = 30  ###### the number of spatial domains.
    file_fold = '/home/cuiyaxuan/spatialLIBD/3.Human_Breast_Cancer' #### to your path
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True) #### project name
    adata.var_names_make_unique()
    model = scSTADE(adata,device=device,n_top_genes=5000)
    adata = model.train()
    radius = 50
    tool = 'mclust' # mclust, leiden, and louvain
    from utils import clustering
    
    if tool == 'mclust':
       clustering(adata, n_clusters, radius=radius, method=tool, refinement=True)
    elif tool in ['leiden', 'louvain']:
       clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)
    
    adata.obs['domain']
    adata.obs['domain'].to_csv("label.csv")


.. parsed-literal::

    Begin to train ST data...


.. parsed-literal::

    
      0%|                                                   | 0/500 [00:00<?, ?it/s]

.. parsed-literal::

    0


.. parsed-literal::

    
      0%|                                           | 1/500 [00:03<29:53,  3.59s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
      0%|▏                                          | 2/500 [00:07<30:10,  3.64s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
      1%|▎                                          | 3/500 [00:10<29:51,  3.61s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
      1%|▎                                          | 4/500 [00:14<29:26,  3.56s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
      1%|▍                                          | 5/500 [00:17<27:51,  3.38s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
      1%|▌                                          | 6/500 [00:20<28:20,  3.44s/it]

.. parsed-literal::

    0


.. parsed-literal::

    
    100%|████████████████████████████████████████▉| 499/500 [29:26<00:03,  3.67s/it]

.. parsed-literal::

    0


.. parsed-literal::

    100%|█████████████████████████████████████████| 500/500 [29:30<00:00,  3.54s/it]


.. parsed-literal::

    Optimization finished for ST data!
    fitting ...
      |======================================================================| 100%


.. code:: ipython3

    import matplotlib as mpl
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import warnings
    import visual
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    warnings.filterwarnings('ignore')
    file_fold = '/home/cuiyaxuan/spatialLIBD/3.Human_Breast_Cancer' #### to your path
    adata = sc.read_visium(file_fold, count_file='filtered_feature_bc_matrix.h5', load_images=True) #### project name
    df_label=pd.read_csv('./label.csv', index_col=0) 
    visual.visual(adata,df_label)



.. parsed-literal::

    #cells after MT filter: 3798



.. image:: 2_Example_BreastSlice_test_files/2_Example_BreastSlice_test_4_1.png
   :width: 362px
   :height: 337px


.. parsed-literal::

    WARNING: saving figure to file figures/showvisualdomainplot_plot.pdf



.. image:: 2_Example_BreastSlice_test_files/2_Example_BreastSlice_test_4_3.png
   :width: 362px
   :height: 337px


