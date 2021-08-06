import csv
import re
import os
import gzip
import tempfile
import pandas as pd
from typing import Dict, List, Optional, IO
from collections import defaultdict
from zipfile import ZipFile, BadZipFile

from eco_kg.transform_utils.transform import Transform
from eco_kg.utils.transform_utils import parse_header, parse_line, write_node_edge_item
from eco_kg.utils import biohub_converter as bc
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
        source_name = 'EOLTraitbank'
        super().__init__(source_name, input_dir, output_dir)  # set some variables

    def run(self, EOLTraitbank: Optional[str] = None) -> None:
        """Method is called and performs needed transformations to process the 
        trait data (EOL Traitbank)
        
        Args:
            EOL_Traitbank: entire contents of EOL TraitBank [EOLTraitbank.zip]
        """
        if not EOLTraitbank:
            EOLTraitbank = os.path.join(self.input_base_dir, 'EOLTraitbank.zip')
        #need to define output files earlier
        node_handle = open(os.path.join(self.output_dir, "eol_traits_nodes.tsv"), 'w')
        edge_handle = open(os.path.join(self.output_dir, "eol_traits_edges.tsv"), 'w')
        self.node_header = ['id', 'name', 'category', 'description', 'provided_by'] #can I add node properties?
        self.edge_header = ['subject', 'edge_label', 'object', 'relation', 'provided_by', 'type']
        node_handle = open(self.output_node_file, 'w')
        edge_handle = open(self.output_edge_file, 'w')
        node_handle.write("\t".join(self.node_header) + "\n")
        edge_handle.write("\t".join(self.edge_header) + "\n")
        self.parse_annotations(node_handle, edge_handle, EOLTraitbank)
        #self.parse_cooccurrence(node_handle, edge_handle, co_occur_zipfile)
            
        #header_items = parse_header(f.readline(), sep=',')
        seen_node: dict = defaultdict(int)
        seen_edge: dict = defaultdict(int)
        # Nodes
        org_node_type = "biolink:OrganismTaxon" 
        trait_node_type = 'biolink:PhenotypicFeature'
        env_node_type = "biolink:AbstractEntity" 
        curie = 'NEED_CURIE'

        #Prefixes
        org_prefix = "EOL:"
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

    def parse_annotations(self, node_handle: IO, edge_handle: IO,
                          EOLTraitbank: str,
                          ) -> None:
        """Parse annotations from traits_all.zip.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            data_file1: Path to traits_all.zip
        Returns:
             None.
        """

        # unzip files and make dataframes
#        with ZipFile('data/raw/EOLtraitbank.zip', 'r') as zf:
#            listOfFileNames = zf.namelist()
#            print(listOfFileNames)
#            print('Unzipping files now')
#            zf.extractall('data/raw/')
#            zf.extract('trait_bank/metadata.csv')
#            try:
#                zf.extractall('data/raw/')
#            except BadZipFile:
#                subprocess.getoutput('zip -FF '+'data/raw/EOLtraitbank.zip')
#                zf.extractall('data/raw/')
#            print('Done')
        listOfFileNames = ['trait_bank/terms.csv', 'trait_bank/term_parents.csv', 'trait_bank/pages.csv', 'trait_bank/inferred.csv', 'trait_bank/traits.csv', 'trait_bank/metadata.csv']
        for filename in listOfFileNames:
            if filename.startswith('.'):
                print(f"skipping file {filename}")
                continue
            file = os.path.join(filename)
            df = pd.read_csv(file, sep=',', low_memory=False)
            if filename == 'trait_bank/pages.csv':
                df = pd.read_csv(file, sep=',', low_memory=False, header=0)
                pages_df = df[['page_id','parent_id','rank','canonical']]
                print(len(df))
                print(len(pages_df))
            if filename == 'trait_bank/traits.csv':
                df = pd.read_csv(file, sep=',', low_memory=False, header=0)
                traits_df = df[['eol_pk','page_id','resource_pk','resource_id','citation','source','predicate','object_page_id','value_uri','normal_measurement','normal_units_uri','normal_units','measurement','units_uri','units','literal','method','remarks','sample_size','name_en']]
                print(len(df))
                print(len(traits_df))
            if filename == 'trait_bank/metadata.csv':
                df = pd.read_csv(file, sep=',', low_memory=False, header=0)
                metadata_df = df[['eol_pk','trait_eol_pk','predicate','value_uri','measurement','units_uri','literal']]
                print(len(df))
                print(len(metadata_df))
            if filename == 'trait_bank/inferred.csv':
                df = pd.read_csv(file, sep=',', low_memory=False, header=0)
                inferred_df = df[['page_id','inferred_trait']]
                print(len(df))
                print(len(inferred_df))
            if filename == 'trait_bank/terms.csv':
                df = pd.read_csv(file, sep=',', low_memory=False, header=0)
                terms_df = df[['uri','name','type']]
                print(len(df))
                print(len(terms_df))
        
        # transform
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
