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
            data_file = self.source_name + '.csv'
        input_file = os.path.join(self.input_base_dir, data_file)
        
        #make directory in data/transformed
        os.makedirs(self.output_dir, data_file)

    def gaf_header(line):
        if line.startswith('!'):
            return True
        return False

        #transform data
        with open(input_file, 'r') as f, \
            open(self.output_node_file, 'w') as node, \
            open(self.output_edge_file, 'w') as edge:
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
            growth_stage_node_type = 'biolink:LifeStage' #Aspect = G gene expressed in growth stage form
            
            #Prefixes - may not need this - only if I'm missing the first part of the CURIE
            org_prefix = "NCBITaxon:"
            gene_prefix = "GO:" #This is a crazy mess, won't be GO
            cellular_component_prefix = "GO:"
            process_prefix = "GO:"
            molecular_function_prefix = "GO:"
            exposure_prefix = 'PECO:'
            anatomy_prefix = 'PO:'
            growth_prefix = 'PO:'
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
            gene_to_exposure_edge_label = 'biolink:increases_expression_of'
            gene_to_exposure_edge_relation = 'RO:0003003'
            gene_to_anatomy_edge_label = 'biolink:expressed_in'
            gene_to_anatomy_edge_relation = 'RO:0002206'
            gene_to_growth_stage_edge_label = 'biolink:expressed_in'
            gene_to_growth_stage_edge_relation = 'RO:0002206'
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
                
                if row['Aspect'] == 'T':
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
                elif row['Aspect'] == 'A':
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
                elif row['Aspect'] == 'E':
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
                elif row['Aspect'] == 'C':
                    GO_id = ['GOid']
                    #create cellular component node
                    if GO_id not in seen_node:
                        write_node_edge_item(fh=node,
                                             header=self.node_header,
                                             data=[go_id,
                                                   go_name,
                                                   go_node_type,
                                                   go_id])
                        seen_node[go_id] += 1
                elif row['Aspect'] == 'F':
                    GO_id = ['GOid']
                    #create molecular function node
                    if GO_id not in seen_node:
                        write_node_edge_item(fh=node,
                                             header=self.node_header,
                                             data=[go_id,
                                                   go_name,
                                                   go_node_type,
                                                   go_id])
                        seen_node[go_id] += 1
                elif row['Aspect'] == 'P':
                    GO_id = ['GOid']
                    #create biological process node
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

            # Write Edge
                # gene to org edge
                if gene_id+tax_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_org_edge_label,
                                                tax_id,
                                                gene_to_org_edge_relation])
                    seen_edge[gene_id+tax_id] += 1

                # gene to cellular component edge
                if gene_id+GO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_cellular_component_edge_label,
                                                GO_id,
                                                gene_to_cellular_component_edge_relation])
                    seen_edge[gene_id+GO_id] += 1
                
                # gene to process edge
                if gene_id+GO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_process_edge_label,
                                                GO_id,
                                                gene_to_process_edge_relation])
                    seen_edge[gene_id+GO_id] += 1

                # gene to molecular function edge
                if gene_id+GO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_molecular_function_edge_label,
                                                GO_id,
                                                gene_to_molecular_function_edge_relation])
                    seen_edge[gene_id+GO_id] += 1

                # gene to exposure edge
                if gene_id+PECO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_exposure_edge_label,
                                                PECO_id,
                                                gene_to_exposure_edge_relation])
                    seen_edge[gene_id+PECO_id] += 1

                # gene to anatomy edge
                if gene_id+PO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_anatomy_edge_label,
                                                PO_id,
                                                gene_to_anatomy_edge_relation])
                    seen_edge[gene_id+PO_id] += 1

                # gene to growth stage edge
                if gene_id+PO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_growth_stage_edge_label,
                                                PO_id,
                                                gene_to_growth_stage_edge_relation])
                    seen_edge[gene_id+PO_id] += 1

                # gene to trait edge
                if gene_id+TO_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[gene_id,
                                                gene_to_trait_edge_label,
                                                TO_id,
                                                gene_to_trait_edge_relation])
                    seen_edge[gene_id+TO_id] += 1

                # trait to org edge
                if TO_id+org_id not in seen_edge:
                    write_node_edge_item(fh=edge,
                                            header=self.edge_header,
                                            data=[TO_id,
                                                trait_to_org_edge_label,
                                                org_id,
                                                trait_to_org_edge_relation])
                    seen_edge[TO_id+org_id] += 1
        # Files write ends

        """
        Implement ROBOT 
        """
        # Convert OWL to JSON for Ontologies
        convert_to_json(self.input_base_dir, 'PO')
        convert_to_json(self.input_base_dir, 'TO')
        convert_to_json(self.input_base_dir, 'PECO')
        # Extract the 'plant' tree from NCBITaxon and convert to JSON
        '''
        NCBITaxon_33090 = viridiplantae
        (Source = http://www.ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_33090)
        '''
        subset_ontology_needed = 'NCBITaxon'
        extract_convert_to_json(self.input_base_dir, subset_ontology_needed, self.subset_terms_file, 'BOT')
        #Can I chop out branches of the hierarchy? NCBITaxon_144314
