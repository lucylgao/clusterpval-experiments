# clusterpval-experiments 

Code to reproduce simulation results, figures, and real data analysis results from the paper "Selective inference for hierarchical clustering" by Lucy L. Gao, Jacob Bien, and Daniela Witten. Preprint available at https://arxiv.org/abs/2012.02936. Download instructions for our R package "clusterpval" can be found at http://www.lucylgao.com/clusterpval/. 

## Organization

### simulation-code  

naive_type1.R produces the results displayed in Figures 1-2. 

Running our_type1.R with this.sim.id between 1-6 produces the results displayed  in Figure 4.  

Running cond.R with this.sim.id between 1-3 produces the results displayed  in Figure 5.

Running power.R with this.sim.id between 1-9 produces the results displayed  in Figures 6-7. 

### simulation-results  

Contains the results from running the code in the simulation-code folder as described above. 

### figures  

Produces Figures 1-8 in the paper. Figures 1-2, 4, 5, and 6-7 depend on the simulation-results folder. 

### figures-code  

Contains the results from calling the code in the figures-code folder.

### real-data-code  

The code to run the real data analysis in Section 6.1 is in penguins.R. Instructions on how to download the data can be found [here](https://allisonhorst.github.io/palmerpenguins/articles/download.html). 

The code to run the real data analysis in Section 6.2 is in zheng.R. The T-cell data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/memory_t), the B-cell data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/b_cells) and the monocyte data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/cd14_monocytes).  


