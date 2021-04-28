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

Essentially just ingests and transforms GAF files

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
        if data_file is None: # = if not data_file
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
        with open(input_file, 'r') as f, \ #is input_file set by download?
            open(self.output_node_file, 'w') as node, \
            open(self.output_edge_file, 'w') as edge, \
            #open(self.subset_terms_file, 'w') as terms_file:   # If need to capture CURIEs for ROBOT STAR extraction

            # write headers (change default node/edge headers if necessary
            node.write("\t".join(self.node_header) + "\n")
            edge.write("\t".join(self.edge_header) + "\n")
            
            header_items = parse_header(f.readline(), sep=',')
            seen_node: dict = defaultdict(int)
            seen_edge: dict = defaultdict(int)
            
            #do I need a dataframe?
        	df = pd.read_csv(input_file, sep='\t', skiprows= lambda x: gaf_header(x), low_memory=False)
        	gaf_df = df.columns['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','Ontology_ID','DB:Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']

            # Nodes
            org_node_type = "biolink:OrganismTaxon"
            gene_node_type = "biolink:GenomicEntity"
            cellular_component_node_type = "biolink:CellularComponent" #Aspect = C
            process_node_type = 'biolink:BiologicalProcess' #Aspect = P
            molecular_function_node_type = 'biolink:MolecularFunction' #Aspect = F
            exposure_node_type = 'biolink:EnvironmentalExposure' #Aspect = E I suspect environment
            anatomy_node_type = 'biolink:AnatomicalEntity' #Aspect = A I suspect anatomy
            trait_node_type = 'biolink:PhenotypicFeature' #Aspect = T I suspect trait
            curie = 'NEED_CURIE'
            
            #Prefixes - may not need this - only if I'm missing the first part of the CURIE
            org_prefix = "NCBITaxon:"
            gene_prefix = "GO:" #This is a crazy mess, won't be GO
            cellular_component_prefix = "GO:"
            process_prefix = "GO:"
            molecular_function_prefix = "GO:"
            exposure_prefix = 'PECO:'
            anatomy_prefix = 'PO:'
            trait_prefix = 'TO:'

            # Edges
            gene_to_org_edge_label = "biolink:in_taxon"
            gene_to_org_edge_relation = "RO:0002162" 
            gene_to_cellular_component_edge_label = "biolink:has_gene_product" # not sure about this one
            gene_to_cellular_component_edge_relation = "RO:" #part of or located in
            gene_to_process_edge_label = "biolink:regulates" # not sure acts upstream of or within - involved in
            gene_to_process_edge_relation = "RO:0011002" #regulates activity of - not sure about this
            gene_to_molecular_function_edge_label = 'biolink:enables'
            gene_to_molecular_function_edge_relation = 'RO:0002327'
            gene_to_exposure_edge_label = 
            gene_to_exposure_edge_relation = 
            gene_to_anatomy_edge_label = 'biolink:expressed_in'
            gene_to_anatomy_edge_relation = 'RO:0002206'
            gene_to_trait_edge_label = 'biolink:has_phenotype'
            gene_to_trait_edge_relation = 'RO:0002200'
            trait_to_org_edge_label = 'biolink:phenotype_of'
            trait_to_org_edge_relation = 'RO:0002201'

            # transform
            for index, row in gaf_df.iterrows():
                tax_id = row['Taxon']
                ontology_id = row['Ontology_ID']
                evidence = row['Evidence']
                gene_id = ['DB_Object_ID']
                #create organism node
                org_id = org_prefix + str(tax_id)
                if org_id not in seen_node:
                    write_node_edge_item(fh=node,
                                         header=self.node_header,
                                         data=[org_id,
                                               org_name,
                                               org_node_type,
                                               org_id])
                    seen_node[org_id] += 1
                #create gene node
				if gene_id not in seen_node:
					write_node_edge_item(fh=node,
										 header=self.node_header,
										 data=[gene_id,
											   gene_name,
											   gene_node_type,
											   gene_id])
					seen_node[gene_id] += 1
				#create org to gene edges
				
                if row['Aspect'] = 'T':
                	trait_id = ['TOid']
					#create trait node
					if trait_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[trait_id,
												   trait_name,
												   trait_node_type,
												   trait_id])
						seen_node[trait_id] += 1
                elif row['Aspect'] = 'A':
                	anatomy_id = ['POid']
					#create anatomy node
					if anatomy_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[anatomy_id,
												   anatomy_name,
												   anatomy_node_type,
												   anatomy_id])
						seen_node[anatomy_id] += 1
                elif row['Aspect'] = 'E':
                	env_id = ['PECOid']
					#create anatomy node
					if env_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[env_id,
												   env_name,
												   env_node_type,
												   env_id])
						seen_node[env_id] += 1
                elif row['Aspect'] = 'C':
                	GO_id = ['GOid']
					#create anatomy node
					if GO_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[go_id,
												   go_name,
												   go_node_type,
												   go_id])
						seen_node[go_id] += 1
                elif row['Aspect'] = 'F':
                	GO_id = ['GOid']
					#create anatomy node
					if GO_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[go_id,
												   go_name,
												   go_node_type,
												   go_id])
						seen_node[go_id] += 1
                elif row['Aspect'] = 'P':
                	GO_id = ['GOid']
					#create anatomy node
					if GO_id not in seen_node:
						write_node_edge_item(fh=node,
											 header=self.node_header,
											 data=[go_id,
												   go_name,
												   go_node_type,
												   go_id])
						seen_node[go_id] += 1
                else:
                	print('Error')
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
