# ENCORE
The periphery and the core properties explain the omnigenic model in the human interactome



ENCORE.ENCORE.Pattern_detection is used to find the peripheral and core properties of a disease neighborhood in the Human Interactome.
We develop pEriphery aNd COre pRopErties (ENCORE) model to disclose component property of a disease neighborhood. 
ENCORE reveal that highly perturbed genes connect into disease-specific cores that drive complex pathological processes, while genes in peripheries are more likely to be shared across diseases.
ENCORE inspires a deeper understanding of interconnectivity of phenotypically related genes, revealing the characteristics of core and periphery.
At present, ENCORE.Pattern_detection has one useful function whose usages are listed below.

-------------------------------------------------------------------------------------------------------------------------------------------------------
ENCORE.ENCORE.ENCORE.Pattern_detection(dir)
Basic: write all results into output files without returning value
Parameters:
			dir:The complete file path of differential expression Genes for a phenotype.
				The input file contains 3 columns: GeneSymbol,P-values and log(fold change).
							Gene	P.Value		logFC
							CCK		7.59E-09	-1.177045455
							CST1	1.12E-08	-3.434015152
							CTSC	2.74E-08	-0.763611111
							...		...			...
returning value:(DiseaseNeighborhood,Core)
			DiseaseNeighborhood: The disease neighborhood stored as a graph by NetworkX in Python. 
			Core:The core component of the disease neighborhood stored as a subgraph by NetworkX in Python.
Output files:
			disease neighborhood.txt: The edge lists of diseas module.
			core and periphery genes.txt:The Genes in core and periphery components.
-------------------------------------------------------------------------------------------------------------------------------------------------------
example one:
import ENCORE.ENCORE as ee
(DiseaseNeighborhood,Core)=ee.Pattern_detection(r'examples\Pulmonary hypertension GSE703.txt')
********************************************************************************
warning: each str referring to the path must be wrapped with r'' just as in the above example
********************************************************************************
