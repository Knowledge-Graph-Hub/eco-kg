eco-KG
================================================
-- UNDER DEVELOPMENT --

Based on the KGH [template](https://github.com/Knowledge-Graph-Hub/kg-dtm-template/generate). 

Documentation
------------------------------------------------

Starting with Planteome - hope to include Monarch Initiative, EOL TraitBank, and more

Currently only contains data from Arabidopsis thaliana, Zea mays, Oryza sativa, Sorghum bicolor, and Populus trichocarpa.

**Components**

- Download: The [download.yaml](download.yaml) contains all the URLs for the source data.
- Transform: The [transform_utils](project_name/transform_utils) contains the code relevant to trnsformations of the raw data into nodes and edges (tsv format)
- Merge: Implementation of the 'merge' function from [KGX](https://github.com/biolink/kgx) using [merge.yaml](merge.yaml) as a source file.

**Utilities**

The code for these are found in the [utils](project_name/utils) folder.

- [ROBOT](https://github.com/ontodev/robot) for transforming ontology.OWL to ontology.JSON
