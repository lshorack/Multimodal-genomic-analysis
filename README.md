# Multimodal-genomic-analysis
This project is based on the [Open Problems in Multimodal Single-Cell Integration](https://www.kaggle.com/competitions/open-problems-multimodal/overview) Kaggle competition. This competition is made up of two parts, and I have focused only on one of them. My work here is about the data acquired by the cite-seq technology. It takes measurements of the quantity of various types of RNA within a cell and the amount of multiple kinds of proteins embed in the cell's surface.

## Problem Statement
This competition aims to develop methodologies for relating multi-modal cellular data. This means finding the relationship between two varieties of data from the same cell. In this case, I am looking at a cell's transcriptome and the presence of specific proteins on the cell's outer membrane. The transcriptome is the RNA released from the nucleus that is ready for ribosomes to build proteins from its sequence.

My goal is to find ways to isolate signals that relate the transcriptome of a cell to the proteins on that cell's surface that are being considered.

## Data
This data consists of two main sets that share the same index, the cells the data was taken from. 10X genomics have built the tools and procedures to simultaneously sample a cell's transcriptome and specific proteins on the cell's surface. In this case, I have a library of 140 different proteins that 10X genomics was able to capture. In general, cells can have thousands of different types of proteins embedded in their outer membrane, not to mention the still more within the cytoplasm. This method captures only a small subset of the proteins this cell produces.

Further descriptions of this data can be found on the competition's data [page](https://www.kaggle.com/competitions/open-problems-multimodal/data).
I also include a blog post describing the process of acquiring this type of data in more detail.

### Note on Data location and associated Kaggle contract
This data can be downloaded from the Kaggle competition [page](https://www.kaggle.com/competitions/open-problems-multimodal/data).

### Some further notes on this data
These cells were taken from three healthy donors who, during medical treatment, were provided a drug that stimulated various progenitor or stem cells from their bone marrow to move into the bloodstream, where they were collected. An explanation for the nature of these cells can be found [here](https://allcells.com/research-grade-tissue-products/mobilized-leukopak/).

Another important detail is that these cells were cultured in the lab, and samples were taken at several points over ten days. These samples were processed via the technology from 10X genomics.

Finally, before posting this data, the Kaggle Competition hosts applied their own normalizing schema to the data. They used different types of normalization for the two modalities, the transcriptome, and the protein samples. To the transcriptome data, they applied a library-size normalization from [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_per_cell.html) followed by a log1p transformation from [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html). To the proteins, they applied a different schema; they used a method called dsb-normalization from [muon](https://muon.readthedocs.io/en/latest/omics/citeseq.html#dsb). (These libraries and methodologies appear to be specific to genomic data.) Rather than getting too deep into the weeds of these normalization methods, I applied a standard scaler to the data, both to the transcriptome and the proteins.

### Data Dictionary
There are three primary datasets I employed for this project. The first is a table of the proteins measured on the cell's membrane called targets. The next is a table of the RNA by amount present within the cell (the transcriptome) called inputs. Finally, there is a metadata table for each cell called meta.


|Feature|Type|Dataset|Description|
|---|---|---|---|
|cell_id|code|all|This is the index used for all tables.|
|Day|code|meta|This is the day the cell was assayed using the 10X Genomic methodology. It is measured from the time the cell was thawed and placed into culture. (This occurred simultaneously for all cells involved.)|
|Donor|code|meta|The ID code for the donor the cell was taken from. (All donors are human.)|
|Cell-type|string|meta|The cell's cell type. The cell types under consideration are all blood progenitor cells; they are cells that generate blood cells of various kinds. It is important to note that these labels were not derived experimentally but from this data and supplied by the Kaggle hosts. They did this to help support EDA. As a result, it does not add any signal not found in the input data.|
|Technology|string|meta|The 10X Genomic methodology was used to assay the cells. This is not relevant to this project. The only cells I consider use the CiteSeq assay that measures transcriptome and proteins in the cellular membrane. The other assay, Multiome, is also developed by 10X Genomics and measures the accessibility of the cell's chromatin and transcriptome. (Again, I do not yet have an analysis of the data associated with the Multiome assay.)|
|Genes|string|inputs|This is the amount of RNA present for a particular gene in the cell’s transcriptome. These values have been normalized as described above using library-size normalization followed by a log1p transformation. An example label is ENSG00000268895_A1BG-AS1; the first component is the gene's unique Ensembl code. Searching for this label on any internet search engine will provide several websites like [GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene=A1BG-AS1) or [Ensembl](https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000268895;r=19:58859117-58866549) that contain additional information about that gene. There are about twenty-two thousand different genes in this table. Together, these columns make up all the inputs to apply to models to predict the proteins, the targets.
|Proteins|string|targets| This is the amount of protein of a particular type on the cell’s outer surface. These values have been normalized as described above using dsb-normalization. An example label is CD86. Most of the labels consist of no more than four characters, although there are some with a few more. There are one hundred and forty such proteins in this table. Together, these columns make up all the targets upon which to apply models.

## Notebook Summary

My Notebooks are best read in the following order:

1. EDA
1. train test split to use across notebooks
1. NVIDIA Rapids correlation data
1. Build Initial Models

### EDA
In this notebook, I look for the obscured signals within this data to build models with these insights.

### train test split to use across notebooks
I applied several methods that took some time to compute, so I needed to run them in google colab. So I decided to train-test split my data in a separate notebook to use across all the other ones.

### NVIDIA Rapids correlation data
#### Dimension Reduction Strategy using correlations and PCA
The models I decided to use were regression models. Since Regression models have a closed form, I thought it would be easier and quicker to produce the 140 models I would need to fit all my targets. So I decided to find the genes that most correlate to each protein, pick the top representatives, and fit my model on those.

Given that this was intensive computationally, I decided to try NVIDIA Rapids, a library similar to pandas, except it is optimized to use a GPU to compute these correlations.

### Build Initial Models
I conclude by building a method for individually fitting all 140 targets to regression models. 

I found that PCA was a valuable tool in this case to augment my dimension reduction strategy involving the correlations found in the previous notebook. 

To reduce the dimension for each model, I use the table of correlations to pick the top one thousand genes most correlated to my target protein. I then run PCA upon all these one thousand predictors picked out of twenty-two possible thousand. After running PCA, I reduce it further by setting the PCA n_components parameter to fifty. This reduces my selected predictors to the top fifty components from the PCA. This produces reasonably robust fits with very low overfitting.


## Conclusions
I have found several options to derive signals from this data. I can use the EDA graphs I created to filter cells based on their behavior by considering how many overall proteins are present in the cell’s membrane compared to the volume of the cell’s transcriptome. This way, I have cells that are more able to drive insight into the underlying biological phenomena. I have also fit all my targets to regression models using several strategies and found a way to significantly reduce the number of dimensions, from twenty-two thousand to fifty, while maintaining a significant amount of signal relevant to many of the target proteins, as demonstrated by my R-Squared values.

## Next Steps
Some of my best conclusions involve filtering cells to reduce noise to produce better models from which new insights can be derived.

In the future, I would like to try more methods to reduce the number of dimensions to consider, expanding on my work in the t-SNE notebook. This includes Multi-Dimensional scaling and some other models under Sci-Kit Learn's manifolds module.
