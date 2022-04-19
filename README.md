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
  ```
  
  If the above setup.py install command doesn't work, try:
  ```
  
  python -m pipinstall.
  
  #Download the data - Use -o to set a path for the download if you don't like the default
  python run.py download
  
  #Transform the data - Use -o and -i to change the default path if you want
  python run.py transform
  
  #Merge the transformed data into a single KG
  python run.py merge
```

Documentation
------------------------------------------------

Starting with Planteome, EOL TraitBank, and the Gene Expression Atlas

Currently only contains data from Arabidopsis, Zea, Oryza, Sorghum, and Populus.

IMPORTANT NOTE ABOUT PO AND TO
Currently unable to transform PO and TO fully within this workflow. I am also using a modified version of TO with added terms. Links to these specific versions are available in download.yaml.

I used ROBOT to transform the PO and TO OWL files to JSON:
* Open a new terminal
* cd to the bin directory in your robot folder
* Run the code below

```
./robot convert --input path/to.owl --format json --output path/to.json
```
Then, using find and replace in PO and TO json files:
* Change plant_anatomy to anatomical_entity. 
* Change plant_structure_development_stage to life_stage. 
* Change plant_trait_ontology to phenotypic_feature.

IMPORTANT NOTE ABOUT EOL TRAITBANK
The only way to get the full data file is to go to the command line and grab the zip file using
```
wget https://editors.eol.org/other_files/SDR/traits_all.zip
```
and then use
```
unzip traits_all.zip
```
Be sure that traits.csv is at least 6GB

Now you can do the transform step.

**Components**

- Download: The [download.yaml](download.yaml) contains all the URLs for the source data.
- Transform: The [transform_utils](project_name/transform_utils) contains the code relevant to trnsformations of the raw data into nodes and edges (tsv format)
- Merge: Implementation of the 'merge' function from [KGX](https://github.com/biolink/kgx) using [merge.yaml](merge.yaml) as a source file.

**Utilities**

The code for these are found in the [utils](project_name/utils) folder.

- [ROBOT](https://github.com/ontodev/robot) for transforming ontology.OWL to ontology.JSON
