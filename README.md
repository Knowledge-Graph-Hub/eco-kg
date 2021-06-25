eco-KG
================================================
-- UNDER DEVELOPMENT --

Based on the KGH [template](https://github.com/Knowledge-Graph-Hub/kg-dtm-template/generate). 

To Use
------------------------------------------------

```
  #Clone the repository
  git clone https://github.com/diatomsRcool/eco-kg
  cd eco-kg
  
  #Set up a virtual environment and install the requirements
  python -m venv kg-env
  source kg-env/bin/activate
  pip install -r requirements.txt
  python setup.py install
  
  #Download the data - Use -o to set a path for the download if you don't like the default
  python run.py download
  
  #Transform the data - Use -o and -i to change the default path if you want
  python run.py transform
  
  #Merge the transformed data into a single KG
  python run.py merge
```

Documentation
------------------------------------------------

Starting with Planteome - hope to include Monarch Initiative, EOL TraitBank, and more

Currently only contains data from Arabidopsis thaliana, Zea mays, Oryza sativa, Sorghum bicolor, and Populus trichocarpa.

Currently unable to transform PO and TO fully within this workflow. I have to separately use ROBOT to transform the OWL file to JSON and change the OBONamespace.

Outside of the virtual environment started above, use robot to convert the owl file (in the data/raw diretory) to json.

```
./robot convert --input path/to.owl --format json --output path/to.json
```
Then, using find and replace change plant_anatomy to anatomical_entity. Change plant_structure_development_stage to life_stage. Change plant_trait_ontology to phenotypic_feature.

**Components**

- Download: The [download.yaml](download.yaml) contains all the URLs for the source data.
- Transform: The [transform_utils](project_name/transform_utils) contains the code relevant to trnsformations of the raw data into nodes and edges (tsv format)
- Merge: Implementation of the 'merge' function from [KGX](https://github.com/biolink/kgx) using [merge.yaml](merge.yaml) as a source file.

**Utilities**

The code for these are found in the [utils](project_name/utils) folder.

- [ROBOT](https://github.com/ontodev/robot) for transforming ontology.OWL to ontology.JSON
