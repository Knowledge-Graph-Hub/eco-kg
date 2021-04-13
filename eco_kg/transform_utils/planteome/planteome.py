import re
import os
from typing import Dict, List, Optional
from collections import defaultdict

from eco_kg.transform_utils.transform import Transform
from eco_kg.utils.transform_utils import parse_header, parse_line, write_node_edge_item
from eco_kg.utils import biohub_converter as bc
#from eco_kg.utils.nlp_utils import *
from eco_kg.utils.robot_utils import *
from kgx.cli.cli_utils import transform
"""
Ingest plant annotations from Planteome

Essentially just ingests and transforms GAF files in file_names.tsv

And extracts the following columns:
    - 
"""

class PlanteomeTransform(Transform):

    def __init__(self, input_dir: str = None, output_dir: str = None) -> None:
        source_name = "Planteome"
        super().__init__(source_name, input_dir, output_dir)  # set some variables
        self.node_header = ['id', 'name', 'category', 'description', 'provided_by']
        self.edge_header = ['subject', 'edge_label', 'object', 'relation', 'provided_by', 'type']

    def run(self, data_file: Optional[str] = None) -> None:
        """Method is called and performs needed transformations to process the 
        plant data (Planteome)
        
        Args:
        	Planteome: entire contents of Planteome []
        """
        if data_file is None:
            data_file = self.source_name + '.csv')
        input_file = os.path.join(self.input_base_dir, data_file)
        
        #make directory in data/transformed
        os.makedirs(self.output_dir, data_file)

	def gaf_header(line):
		if line.startswith('!'):
			return True
		return False

    def parse_annotations(self, node_handle: IO, edge_handle: IO,
                          data_file1: str,
                          ) -> None:
        """Parse annotations from Planteome.
        Args:
            node_handle: File handle for nodes.csv.
            edge_handle: File handle for edges.csv.
            data_file1: Path to Planteome GAF
        Returns:
             None.
        """
        #transform data
        with open(input_file, 'r') as f, \
            open(self.output_node_file, 'w') as node, \
            open(self.output_edge_file, 'w') as edge, \
            #open(self.subset_terms_file, 'w') as terms_file:   # If need to capture CURIEs for ROBOT STAR extraction

            # write headers (change default node/edge headers if necessary
            node.write("\t".join(self.node_header) + "\n")
            edge.write("\t".join(self.edge_header) + "\n")
            
            header_items = parse_header(f.readline(), sep=',')
            seen_node: dict = defaultdict(int)
            seen_edge: dict = defaultdict(int)
            
        	df = pd.read_csv(input_file, sep='\t', skiprows= lambda x: gaf_header(x), low_memory=False)
        	gaf_df = df.columns['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','GO_ID','DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']

            # Nodes
            org_node_type = "biolink:OrganismTaxon"
            gene_node_type = "biolink:GenomicEntity"
            cellular_component_node_type = "biolink:CellularComponent"
            process_node_type = 'biolink:BiologicalProcess'
            molecular_function_node_type = 'biolink:MolecularFunction'
            curie = 'NEED_CURIE'
            
            #Prefixes
            org_prefix = "NCBITaxon:"
            gene_prefix = "GO:"
            cellular_component_prefix = "GO:"
            process_prefix = "GO:"
            molecular_function_prefix = "GO:"

            # Edges
            gene_to_org_edge_label = "biolink:in_taxon" #  [org_name -> cell_shape, metabolism]
            gene_to_org_edge_relation = "RO:0002200" #  [org_name -> has phenotype -> cell_shape, metabolism]
            gene_to_cellular_component_edge_label = "biolink:has_gene_product" # [org_name -> carbon_substrate]
            gene_to_cellular_component_edge_relation = "RO:0002438" # [org_name -> 'trophically interacts with' -> carbon_substrate]
            gene_to_process_edge_label = "biolink:location_of" # [org -> isolation_source]
            gene_to_process_edge_relation = "RO:0001015" #[org -> location_of -> source]
            gene_to_molecular_function_label = 
            gene_to_molecular_function_relation = 

            
            
            # transform
            for line in f:
                """
                This dataset is a csv and also has commas 
                present within a column of data. 
                Hence a regex solution
                """
                # transform line into nodes and edges
                # node.write(this_node1)
                # node.write(this_node2)
                # edge.write(this_edge)
                

                line = re.sub(r'(?!(([^"]*"){2})*[^"]*$),', '|', line) # alanine, glucose -> alanine| glucose
                items_dict = parse_line(line, header_items, sep=',')

                org_name = items_dict['org_name']
                tax_id = items_dict['tax_id']
                metabolism = items_dict['metabolism']
                carbon_substrates = set([x.strip() for x in items_dict['carbon_substrates'].split('|')])
                cell_shape = items_dict['cell_shape']
                isolation_source = set([x.strip() for x in items_dict['isolation_source'].split('|')])
                

            # Write Node ['id', 'entity', 'category']
                # Write organism node 
                org_id = org_prefix + str(tax_id)
                if not org_id.endswith(':na') and org_id not in seen_node:
                    write_node_edge_item(fh=node,
                                         header=self.node_header,
                                         data=[org_id,
                                               org_name,
                                               org_node_type,
                                               org_id])
                    seen_node[org_id] += 1
                    if org_id.startswith('NCBITaxon:'):
                        terms_file.write(org_id + "\n")

                # Write chemical node
                for chem_name in carbon_substrates:
                    chem_curie = curie
                    #chem_node_type = chem_name

                    # Get relevant NLP results
                    if chem_name != 'NA':
                        relevant_tax = oger_output.loc[oger_output['TaxId'] == int(tax_id)]
                        relevant_chem = relevant_tax.loc[relevant_tax['TokenizedTerm'] == chem_name]
                        if len(relevant_chem) == 1:
                            chem_curie = relevant_chem.iloc[0]['CURIE']
                            chem_node_type = relevant_chem.iloc[0]['Biolink']
                        

                    if chem_curie == curie:
                        chem_id = chem_prefix + chem_name.lower().replace(' ','_')
                    else:
                        chem_id = chem_curie

                    
                    if  not chem_id.endswith(':na') and  chem_id not in seen_node:
                        write_node_edge_item(fh=node,
                                            header=self.node_header,
                                            data=[chem_id,
                                                chem_name,
                                                chem_node_type,
                                                chem_curie])
                        seen_node[chem_id] += 1

                # Write shape node
                shape_id = shape_prefix + cell_shape.lower()
                if  not shape_id.endswith(':na') and shape_id not in seen_node:
                    write_node_edge_item(fh=node,
                                         header=self.node_header,
                                         data=[shape_id,
                                               cell_shape,
                                               shape_node_type,
                                               curie])
                    seen_node[shape_id] += 1

                # Write source node
                for source_name in isolation_source:
                    #   Collapse the entity
                    #   A_B_C_D => [A, B, C, D]
                    #   D is the entity of interest
                    source_name_split = source_name.split('_')
                    source_name_collapsed = source_name_split[-1]
                    env_curie = curie
                    env_term = source_name_collapsed
                    source_node_type = "" # [isolation_source] left blank intentionally

                    # Get information from the environments.csv (unique_env_df)
                    relevant_env_df = unique_env_df.loc[unique_env_df['Type'] == source_name]

                    if len(relevant_env_df) == 1:
                            '''
                            If multiple ENVOs exist, take the last one since that would be the curie of interest
                            after collapsing the entity.
                            TODO(Maybe): If CURIE is 'nan', it could be sourced from OGER o/p (ENVO backend)
                                  of environments.csv
                            '''
                            env_curie = str(relevant_env_df.iloc[0]['ENVO_ids']).split(',')[-1].strip()
                            env_term = str(relevant_env_df.iloc[0]['ENVO_terms']).split(',')[-1].strip()
                            if env_term == 'nan':
                                env_curie = curie
                                env_term = source_name_collapsed
                            
                                 

                    #source_id = source_prefix + source_name.lower()
                    if env_curie == curie:
                        source_id = source_prefix + source_name_collapsed.lower()
                    else:
                        source_id = env_curie
                        if source_id.startswith('CHEBI:'):
                            source_node_type = chem_node_type

                    if  not source_id.endswith(':na') and source_id not in seen_node:
                        write_node_edge_item(fh=node,
                                            header=self.node_header,
                                            data=[source_id,
                                                env_term,
                                                source_node_type,
                                                env_curie])
                        seen_node[source_id] += 1

                


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
