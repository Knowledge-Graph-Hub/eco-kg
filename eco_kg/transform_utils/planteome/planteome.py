import re
import os
import pandas as pd
from typing import Dict, List, Optional
from collections import defaultdict
import gzip

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
        self.node_header = ['id', 'name', 'category', 'provided_by']
        self.edge_header = ['subject', 'edge_label', 'object', 'predicate', 'provided_by']

    def run(self, data_files: List[str] = None):
        source_name = 'Planteome'#how do I pull this from the top?
        """Method is called and performs needed transformations to process the 
        plant data (Planteome)
        
        Args:
            Planteome: entire contents of Planteome []
        """
        """
        Get rice gene ID map
        """
        rice_gene_ids = {}
        with gzip.open(os.path.join(self.input_base_dir, 'rice_map.gz'),'rt') as rg:
            for line in rg.readlines():
                line = line.strip('\n')
                row = line.split('\t')
                Os = row[0]
                LOC = [row[1]]
                if Os == 'None' or LOC == 'None':
                    continue
                OC = []
                for x in LOC:
                    x = x.split('.')[0]
                    if x not in OC:
                        OC.append(x)
                rice_gene_ids[Os] = OC

        corn_gene_ids = {}
        with gzip.open(os.path.join(self.input_base_dir, 'corn_map.gz'), 'rt') as cg:
            for line in cg.readlines():
                k = []
                v = []
                line = line.strip('\n')
                row = line.split('\t')
                for g in row:
                    if g == '':
                        continue
                    if 'B73v' in g:
                        g = re.sub('B73v\d_','',g)
                    if 'Zm00001eb' in g and g not in v:
                        v.append(g)
                    else:
                        if g not in k:
                            k.append(g)
                for i in k:
                    corn_gene_ids[i] = v

        if not data_files: # = if not data_file
            data_files = []
            for file in os.listdir(self.input_base_dir):
                if file.endswith('.assoc'):
                    data_files.append(file)

        seen_node: dict = defaultdict(int)
        seen_edge: dict = defaultdict(int)
        with open(self.output_node_file, 'w') as node, \
            open(self.output_edge_file, 'w') as edge:
            node.write("\t".join(self.node_header) + "\n")
            edge.write("\t".join(self.edge_header) + "\n")

            for data_file in data_files:
                input_file = os.path.join(self.input_base_dir, data_file)
                print(str(data_file))
                #make directory in data/transformed
                os.makedirs(self.output_dir, exist_ok=True)
            #transform data
                with open(input_file, 'r') as f:
                    #header_items = parse_header(f.readline(), sep=',')
                    header_row = ['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','Ontology_ID','DB_Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']
                    gaf_df = pd.read_csv(input_file, sep='\t', comment='!', names=header_row, low_memory=False)
                    #gaf_df = df.columns['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','Ontology_ID','DB_Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']

                    # Nodes
                    org_node_type = "biolink:OrganismTaxon"
                    gene_node_type = "biolink:GenomicEntity"
                    cellular_component_node_type = "biolink:CellularComponent" #Aspect = C
                    process_node_type = 'biolink:BiologicalProcess' #Aspect = P
                    molecular_function_node_type = 'biolink:MolecularFunction' #Aspect = F
                    anatomy_node_type = 'biolink:AnatomicalEntity' #Aspect = A The root class in PO is plant anatomical entity - Is this a problem?
                    trait_node_type = 'biolink:PhenotypicFeature' #Aspect = T I suspect trait
                    growth_stage_node_type = 'biolink:LifeStage' #Aspect = G (might be S) gene expressed in growth stage form
            
                    #Prefixes - may not need this - only if I'm missing the first part of the CURIE
                    org_prefix = "NCBITaxon:"
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
                    gene_to_cellular_component_edge_relation = "RO:0002205" #part of or located in
                    gene_to_process_edge_label = "biolink:regulates" # not sure acts upstream of or within - involved in
                    gene_to_process_edge_relation = "RO:0011002" #regulates activity of - not sure about this
                    gene_to_molecular_function_edge_label = 'biolink:enables'
                    gene_to_molecular_function_edge_relation = 'RO:0002327'
                    gene_to_anatomy_edge_label = 'biolink:expressed_in'
                    gene_to_anatomy_edge_relation = 'RO:0002206'
                    gene_to_growth_stage_edge_label = 'biolink:expressed_in'
                    gene_to_growth_stage_edge_relation = 'RO:0002206'
                    gene_to_trait_edge_label = 'biolink:has_phenotype'
                    gene_to_trait_edge_relation = 'RO:0002200'
                    gene_to_orth_edge_label = 'biolink:orthologous_to'
                    gene_to_orth_edge_relation = 'RO:HOM0000017'
                    trait_to_org_edge_label = 'biolink:phenotype_of'
                    trait_to_org_edge_relation = 'RO:0002201'

                    # transform
                    for index, row in gaf_df.iterrows():
                        tax_id = row['Taxon'].split(':')[1]
                        ontology_id = row['Ontology_ID']
                        evidence = row['Evidence_Code']
                        type = row['DB_Object_Type']
                        db = row['DB']
                        gene_name = row['DB_Object_Symbol']
                        label = row['DB_Object_Symbol']
                        #source = row['DB']
                        if tax_id == '3702':
                            gene_id = str(row['DB_Object_Name'])
                            org_name = 'Arabidopsis thaliana'
                            if 'AT' in gene_id:
                                if '.' in gene_id:
                                    gene_id = gene_id.split('.')[0]
                                elif len(gene_id) > 9:
                                    gene_id = row['DB_Object_Synonym'].split('|')[0]
                            else:
                                gene_id = str(row['DB_Object_Synonym']).split('|')[0]
                                if 'AT' not in gene_id:
                                    gene_id = row['DB_Object_ID']
                        elif tax_id == '4530':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa'
                        elif tax_id == '39947':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa japonica'
                        elif tax_id == '4558':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Sorghum bicolor'
                        elif tax_id == '3694':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Populus trichocarpa'
                        elif tax_id == '4577':#need to add another conditional to make it look in synonyms
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays'
                        elif tax_id == '381124':#need to add another conditional to make it look in synonyms
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays mays'
                        else:
                            gene_id = row['DB_Object_ID']
                            org_name = 'taxon name'
                        provided_by = db + '-' + source_name
                        #create organism node
                        org_id = org_prefix + str(tax_id)
                        #if the data are from rice...
                        if tax_id == '4530' or tax_id == '39947':#changing Os gene IDs to LOC gene IDs
                            if 'LOC' not in gene_id:
                                gene_id = row['DB_Object_Name']
                            if 'LOC' not in gene_id:
                                gene_id = row['DB_Object_Symbol']
                            if 'LOC' not in gene_id and gene_id in rice_gene_ids:
                                #print(gene_id)
                                gid = rice_gene_ids[gene_id]
                                #print(gid)
                        #if the data are from corn...
                        if tax_id == '381124' or tax_id == '4577':
                            #print(row)
                            if isinstance(gene_id, int):
                                try:
                                    gene_id = row['DB_Object_Synonym'].split('|')[0]
                                    if 'GRMZM' not in gene_id:
                                        gene_id = row['DB_Object_Name']
                                        if 'GRMZM' not in gene_id:
                                            gene_id = row['DB_Object_Symbol']
                                        #print(gene_id)
                                except AttributeError:
                                    #print(gene_id)
                                    gene_id = row['DB_Object_Symbol']
                            #print(gene_id)
                            if '_' in gene_id:
                                gene_id = gene_id.split('_')[0]
                            if 'Zm00001eb' not in gene_id and gene_id in corn_gene_ids:
                                gene_id = corn_gene_ids[gene_id]
                        #create the organism node
                        if org_id not in seen_node:
                            write_node_edge_item(fh=node,
                                                 header=self.node_header,
                                                 data=[org_id,
                                                       org_name,
                                                       org_node_type,
                                                       provided_by])
                            seen_node[org_id] += 1
                        #create gene node
                        if isinstance(gene_id, str):
                            genes = [gene_id]
                        if isinstance(gene_id, int):
                            gi = str(gene_id)
                            genes = [gi]
                        #print(genes)
                        for g in genes:
                            if 'AGI_LocusCode' in g:
                                g = g.split(':')[1]
                            if g not in seen_node:
                                g = str(g)
                                write_node_edge_item(fh=node,
                                                     header=self.node_header,
                                                     data=[g,
                                                           gene_name,
                                                           gene_node_type,
                                                           provided_by])
                                seen_node[g] += 1
                        #create all other node types
                        if ontology_id not in seen_node:
                            if row['Aspect'] == 'T':
                                node_type = trait_node_type
                                label = row['DB_Object_Name']
                            elif row['Aspect'] == 'A':
                                node_type = anatomy_node_type
                                label = 'Need PO label'
                            elif row['Aspect'] == 'G':
                                node_type = growth_stage_node_type
                                label = 'Need PO label'
                            elif row['Aspect'] == 'C':
                                node_type = cellular_component_node_type
                                label = 'Need GO label'
                            elif row['Aspect'] == 'F':
                                node_type = molecular_function_node_type
                                label = 'Need GO label'
                            elif row['Aspect'] == 'P':
                                node_type = process_node_type
                                label = 'Need GO label'
                            else:
                                print('Error. New Aspect.')
                                print(row['Aspect'])
                            write_node_edge_item(fh=node,
                                                 header=self.node_header,
                                                 data=[ontology_id,
                                                       label,
                                                       node_type,
                                                       provided_by])
                            seen_node[g] += 1

                        #create additional nodes for orthologs
                        if 'ortholog' in data_file:
                            orth = row['With_or_From']
                            if '|' in orth:
                                orth = orth.split('|')
                                h = []
                                for t in orth:
                                    a = t.split(':')[1]
                                    if a not in h:
                                        h.append(a)
                                orth = h
                            else:
                                orth = orth.split(':')[1]
                                orth = [orth]
                            #print(orth)
                            for o in orth:
                                if o not in seen_node:
                                    gene_name = 'none'
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[o,
                                                               gene_name,
                                                               gene_node_type,
                                                               provided_by])
                                    seen_node[g] += 1
                    # Write Edge
                        # gene to org edge
                        if isinstance(gene_id, str):
                            genes = [gene_id]
                        if isinstance(gene_id, int):
                            gi = str(gene_id)
                            genes = [gene_id]
                        for g in genes:
                            #print(g)
                            if str(g)+str(tax_id) not in seen_edge:
                                #print(str(g)+str(tax_id))
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[g,
                                                            gene_to_org_edge_label,
                                                            tax_id,
                                                            gene_to_org_edge_relation,
                                                            provided_by])
                                seen_edge[str(g)+str(tax_id)] += 1

                            # gene to cellular component edge
                            if str(g)+ontology_id not in seen_edge:
                                #print(str(g)+ontology_id)
                                if row['Aspect'] == 'T':
                                    edge_label = gene_to_trait_edge_label
                                    edge_relation = gene_to_trait_edge_relation
                                elif row['Aspect'] == 'A':
                                    edge_label = gene_to_anatomy_edge_label
                                    edge_relation = gene_to_anatomy_edge_relation
                                elif row['Aspect'] == 'C':
                                    edge_label = gene_to_cellular_component_edge_label
                                    edge_relation = gene_to_cellular_component_edge_relation
                                elif row['Aspect'] == 'P':
                                    edge_label = gene_to_process_edge_label
                                    edge_relation = gene_to_process_edge_relation
                                elif row['Aspect'] == 'F':
                                    edge_label = gene_to_molecular_function_edge_label
                                    edge_relation = gene_to_molecular_function_edge_relation
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[g,
                                                            edge_label,
                                                            ontology_id,
                                                            edge_relation,
                                                            provided_by])
                                #print(edge_label)
                                seen_edge[str(g)+ontology_id] += 1
                        #ortholog edges
                        if 'ortholog' in data_file:
                            for g in genes:
                                #print(g)
                                for o in orth:
                                    if str(g)+str(o) not in seen_edge:
                                        write_node_edge_item(fh=edge,
                                                                header=self.edge_header,
                                                                data=[g,
                                                                    gene_to_orth_edge_label,
                                                                    o,
                                                                    gene_to_orth_edge_relation,
                                                                    provided_by])
                                        seen_edge[str(g)+str(tax_id)] += 1

                        # trait to org edge
                        if 'TO' in ontology_id:
                            if ontology_id+org_id not in seen_edge:
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[ontology_id,
                                                            trait_to_org_edge_label,
                                                            org_id,
                                                            trait_to_org_edge_relation,
                                                            provided_by])
                                #print(trait_to_org_edge_label)
                                seen_edge[ontology_id+org_id] += 1
                # Files write ends

        """
        Implement ROBOT 
        """
        # Convert OWL to JSON for Ontologies
        convert_to_json(self.input_base_dir, 'PO')
        convert_to_json(self.input_base_dir, 'TO')
        # Extract the 'plant' tree from NCBITaxon and convert to JSON
        '''
        NCBITaxon_33090 = viridiplantae
        (Source = http://www.ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_33090)
        '''
        subset_ontology_needed = 'NCBITaxon'
        extract_convert_to_json(self.input_base_dir, subset_ontology_needed, 'NCBITaxon:33090','TOP')#maybe should be TOP
        #Can I chop out branches of the hierarchy? NCBITaxon_144314
