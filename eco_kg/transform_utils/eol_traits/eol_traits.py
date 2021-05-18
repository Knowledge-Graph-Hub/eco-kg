import csv
import re
import os
import gzip
import pandas as pd
from typing import Dict, List, Optional
from collections import defaultdict
from zipfile import ZipFile

from eco_kg.transform_utils.transform import Transform
from eco_kg.utils.transform_utils import parse_header, parse_line, write_node_edge_item
from eco_kg.utils import biohub_converter as bc
from eco_kg.utils.nlp_utils import *
from eco_kg.utils.robot_utils import *
from eco_kg.utils.transform_utils import unzip_to_tempdir
from kgx.cli.cli_utils import transform
"""
Ingest traits dataset (EOL TraitBank)

Essentially just ingests and transforms this file:
https://editors.eol.org/other_files/SDR/traits_all.zip

And extracts the following columns:
    - 
"""

class EOLTraitsTransform(Transform):

    def __init__(self, input_dir: str = None, output_dir: str = None) -> None:
        source_name = "EOL_TraitBank"
        super().__init__(source_name, input_dir, output_dir)  # set some variables

    def run(self, EOL_TraitBank: Optional[str] = None) -> None:
        """Method is called and performs needed transformations to process the 
        trait data (EOL TraitBank)
        
        Args:
        	EOL_TraitBank: entire contents of EOL TraitBank [traits_all.zip]
        """
        if not EOL_TraitBank:
            EOL_TraitBank = os.path.join(self.input_base_dir, "traits_all.zip")
        #need to define output files earlier
        node_handle = open(os.path.join(self.output_dir, "eol_traits_nodes.tsv"), 'w')
        edge_handle = open(os.path.join(self.output_dir, "eol_traits_edges.tsv"), 'w')
        self.node_header = ['id', 'name', 'category', 'description', 'provided_by'] #can I add node properties?
        self.edge_header = ['subject', 'edge_label', 'object', 'relation', 'provided_by', 'type']
        node_handle = open(self.output_node_file, 'w')
        edge_handle = open(self.output_edge_file, 'w')
        node_handle.write("\t".join(self.node_header) + "\n")
        edge_handle.write("\t".join(self.edge_header) + "\n")
        self.parse_annotations(node_handle, edge_handle, EOL_TraitBank)
        #self.parse_cooccurrence(node_handle, edge_handle, co_occur_zipfile)

    def parse_annotations(self, node_handle: IO, edge_handle: IO,
                          data_file1: str,
                          ) -> None:
        """Parse annotations from CORD-19_1_5.zip.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            data_file1: Path to traits_all.zip
        Returns:
             None.
        """
        #progress bar showing the unzipping of files on cmd line
        pbar = tqdm(total=1, desc="Unzipping files")

        # unzip to tmpdir, remove after use, to avoid cluttering raw/ with processed
        # data
        with tempfile.TemporaryDirectory(dir=self.input_base_dir) as tmpdir:
            unzip_to_tempdir(data_file1, tmpdir)
            pbar.update(1)
            pbar.close()


            subsets = ['EOL_TraitBank']
            for subset in subsets:
                subset_dir = os.path.join(tmpdir, subset)
                for filename in tqdm(os.listdir(subset_dir)):
                    if filename.startswith('.'):
                        print(f"skipping file {filename}")
                        continue
                    file = os.path.join(subset_dir, filename)
                    doc = json.load(open(file))
                    df = pd.read_csv(file, sep=',', low_memory=False)
                    if filename = 'pages':
                    	self.parse_pages(node_handle, edge_handle, df)
                    if filename = 'traits':
                    	self.parse_traits(node_handle, edge_handle, df)
                    if filename = 'metadata':
                    	self.parse_metadata(node_handle, edge_handle, df)
                    if filename = 'inferred':
                    	self.parse_inferred(node_handle, edge_handle, df)
                    if filename = 'terms':
                    	self.parse_terms(node_handle, edge_handle, df)

    def parse_pages(self, node_handle, edge_handle, df) -> None:
        """Parse a pandas dataframe.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            df: pandas dataframe
            subset: The subset name for this dataset.
        Returns:
            None.
        """
        pages_df = df['page_id','canonical','rank','parent','trait','inferred_trait']

    def parse_traits(self, node_handle, edge_handle, df) -> None:
        """Parse a pandas dataframe.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            df: pandas dataframe
            subset: The subset name for this dataset.
        Returns:
            None.
        """
        traits_df = df['eol_pk','resource_pk','citation','source','predicate','object_term','object_page','sex_term','lifestage_term']
        #find specific predicates
        
    def parse_metadata(self, node_handle, edge_handle, df) -> None:
        """Parse a pandas dataframe.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            df: pandas dataframe
            subset: The subset name for this dataset.
        Returns:
            None.
        """
        metadata_df = df['eol_pk','trait_eol_pk','predicate','value_uri','measurement','units_uri','literal']
        
    def parse_inferred(self, node_handle, edge_handle, df) -> None:
        """Parse a pandas dataframe.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            df: pandas dataframe
            subset: The subset name for this dataset.
        Returns:
            None.
        """
        inferred_df = df['page_id','inferred_trait']
        
    def parse_terms(self, node_handle, edge_handle, df) -> None:
        """Parse a pandas dataframe.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            df: pandas dataframe
            subset: The subset name for this dataset.
        Returns:
            None.
        """
        terms_df = df['uri','name','type','parent_uri']
"""
        # transform data, something like:
        with open(input_file, 'r') as f, \
                open(self.output_node_file, 'w') as node, \
                open(self.output_edge_file, 'w') as edge, \
                open(self.subset_terms_file, 'w') as terms_file:

            # write headers (change default node/edge headers if necessary
            node.write("\t".join(self.node_header) + "\n")
            edge.write("\t".join(self.edge_header) + "\n")
            
            header_items = parse_header(f.readline(), sep=',')
            
            seen_node: dict = defaultdict(int)
            seen_edge: dict = defaultdict(int)
"""
	# transform
	
	header_items = parse_header(f.readline(), sep=',')
	seen_node: dict = defaultdict(int)
	seen_edge: dict = defaultdict(int)
	# Nodes
	org_node_type = "biolink:OrganismTaxon" 
	trait_node_type = 'biolink:PhenotypicFeature'
	env_node_type = "biolink:AbstractEntity" 
	curie = 'NEED_CURIE'

	#Prefixes
	org_prefix = "EOL:"
	trait_prefix = "Carbon:"#not sure what this should be
	env_prefix = "ENVO:"

	# Edges
	org_to_trait_edge_label = "biolink:has_phenotype" 
	org_to_trait_edge_relation = "RO:0002200" 
	org_to_parent_edge_label = 'biolink:OrganismTaxonToOrganismTaxonSpecialization' #need to make a new biolink predicate
	org_to_parent_edge_relation = 'RO:'
	org_to_org_edge_label = "biolink:interacts_with"
	org_to_org_edge_relation = "RO:" 
	org_to_env_edge_label = "biolink:location_of" 
	org_to_env_edge_relation = "RO:" 
	
	#need a function to pare down the terms file and make a dict/json of what I want to use
	
	#write organism nodes
	for index, row in pages_df.iterrows():
		org_id = org_prefix + str(pages_df['page_id'])
		if org_id not in seen_node:
			write_node_edge_item(fh=node,
								 header=self.node_header,
								 data=[org_id,
									   org_name,
									   org_node_type,
									   org_id])
			seen_node[org_id] += 1
		if org_id+parent_id not in seen_edge:
			parent_id = org_prfix + str(pages_df['parent'])
			write_node_edge_item(fh=edge,
									header=self.edge_header,
									data=[org_id,
										org_to_parent_edge_label,
										parent_id,
										org_to_parent_edge_relation])
			seen_edge[org_id+parent_id] += 1

	# Write trait node
	for index, row in traits_df.iterrows():
		trait_id = traits_df['eol_pk']
		trait_uri = traits_df['predicate']
		value_uri = traits_df['value_uri']
		org_id = traits_df['page_id']
		#need to write rest of function when I transform terms.csv to json
		#this will also create env nodes
		#need to figure out if I should have traits and values separate nodes or not. What does kg-microbe do?
		if traits_df['type'] == 'measurement':
			trait = trait_dict[trait_uri]
			value = value_dict[value_uri]
			if trait+value not in seen_node:
				write_node_edge_item(fh=node,
									 header=self.node_header,
									 data=[org_id,
										   org_name,
										   org_node_type,
										   org_id])
				seen_node[trait_id] += 1
			if org_id+trait_id not in seen_edge:
				trait_id

	# Write Edge
		# org-chem edge
		if not chem_id.endswith(':na') and org_id+chem_id not in seen_edge:
			write_node_edge_item(fh=edge,
									header=self.edge_header,
									data=[org_id,
										org_to_chem_edge_label,
										chem_id,
										org_to_chem_edge_relation])
			seen_edge[org_id+chem_id] += 1

		# org-shape edge
		if  not shape_id.endswith(':na') and org_id+shape_id not in seen_edge:
			write_node_edge_item(fh=edge,
									header=self.edge_header,
									data=[org_id,
										org_to_shape_edge_label,
										shape_id,
										org_to_shape_edge_relation])
			seen_edge[org_id+shape_id] += 1
		
		# org-source edge
		if not source_id.endswith(':na') and org_id+source_id not in seen_edge:
			write_node_edge_item(fh=edge,
									header=self.edge_header,
									data=[org_id,
										org_to_source_edge_label,
										source_id,
										org_to_source_edge_relation])
			seen_edge[org_id+source_id] += 1
	# Files write ends

	# Extract the 'cellular organismes' tree from NCBITaxon and convert to JSON
	'''
	NCBITaxon_131567 = cellular organisms 
	(Source = http://www.ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_131567)
	'''
	subset_ontology_needed = 'NCBITaxon'
	extract_convert_to_json(self.input_base_dir, subset_ontology_needed, self.subset_terms_file)
