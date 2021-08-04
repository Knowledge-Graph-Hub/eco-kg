import re
import os
import pandas as pd
from typing import Dict, List, Optional
from collections import defaultdict
import gzip
import json

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
        self.edge_header = ['subject', 'predicate', 'object', 'relation', 'provided_by']

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
                
        #Get corn gene id map
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
                    
        #Get plant germplasm categorical traits mapped_keys
        with open(os.path.join(self.input_base_dir, 'plant_traits.txt'), 'r') as pt:
            plant_trait_ids = json.load(pt)
        #Get plant germplasm numerical traits list
#        plant_traits_num = ['TO:0000269','TO:0000933','TO:0000934','TO:0002757','TO:0002758','TO:0000402','TO:0000734','TO:0000344','TO:0002644','TO:0000602','TO:0001024','TO:0000696','TO:0000919','TO:0000354','TO:0000207','TO:0001002','TO:0000040','TO:0000604','TO:0000797','TO:0000687','TO:0000346','TO:0002746','TO:0000658','TO:0000598']
        num_only = ['TO:0000269','TO:0000933','TO:0000934','TO:0002757','TO:0002758','TO:0000402','TO:0000734','TO:0000344','TO:0002644','TO:0000602','TO:0001024','TO:0000696','TO:0000919','TO:0000354','TO:0000207','TO:0001002','TO:0000040','TO:0000604','TO:0000797','TO:0000687','TO:0000346','TO:0002746','TO:0000658','TO:0000598','TO:0000531','TO:0000072','TO:0001035','TO:0000196','TO:0000734','TO:0000462','TO:0000382','TO:0000138','TO:0000240','TO:0000433','TO:0000431']
        ignore_traits = ['TO:0000439','TO:0020010','TO:0000261','TO:0000664','TO:0000273','TO:0000468','TO:0000315','TO:0000148']
        compound_traits = ['TO:0000969','TO:0002629','TO:0000068']
        with open(os.path.join(self.input_base_dir, 'plant_numerical_traits.txt'), 'r') as pc:
            plant_traits_num = json.load(pc)

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
                file_type = 'genetic'
                if 'germplasm' in data_file:
                    file_type = 'germplasm'
                if 'qtl' in data_file:
                    file_type = "qtl"
                #make directory in data/transformed
                os.makedirs(self.output_dir, exist_ok=True)
            #transform data
                with open(input_file, 'r') as f:
                    #header_items = parse_header(f.readline(), sep=',')
                    header_row = ['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','Ontology_ID','DB_Reference','Evidence_Code','With_or_From','Aspect','DB_Object_Name','DB_Object_Synonym','DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension','Gene_Product_Form_ID']
                    gaf_df = pd.read_csv(input_file, sep='\t', comment='!', names=header_row, low_memory=False)

                    # Nodes
                    org_node_type = "biolink:OrganismTaxon"
                    gene_node_type = "biolink:GenomicEntity"
                    cellular_component_node_type = "biolink:CellularComponent" #Aspect = C
                    process_node_type = 'biolink:BiologicalProcess' #Aspect = P
                    molecular_function_node_type = 'biolink:MolecularFunction' #Aspect = F
                    anatomy_node_type = 'biolink:AnatomicalEntity' #Aspect = A The root class in PO is plant anatomical entity - Is this a problem?
                    trait_node_type = 'biolink:PhenotypicFeature' #Aspect = T I suspect trait
                    growth_stage_node_type = 'biolink:LifeStage' #Aspect = G (might be S) gene expressed in growth stage form
                    qtl_node_type = 'biolink:GenomicEntity'
                    germplasm_node_type = 'biolink:OrganismalEntity'
                    trait_value_node_type = 'biolink:PhenotypicQuality'
                    value_node_type = 'biolink:QuantityValue'
                    #unit_node_type = 'biolink:Unit'
            
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
                    gene_to_cellular_component_edge_label = "biolink:has_gene_product"
                    gene_to_cellular_component_edge_relation = "RO:0002205"
                    gene_to_process_edge_label = "biolink:regulates"
                    gene_to_process_edge_relation = "RO:0011002"
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
                    qtl_to_org_edge_label = "biolink:in_taxon"
                    qtl_to_org_edge_relation = "RO:0002162"
                    germplasm_to_org_edge_label = "biolink:in_taxon"
                    germplasm_to_org_edge_relation = "RO:0002162"
                    germplasm_to_trait_edge_label = 'biolink:has_attribute'
                    germplasm_to_trait_edge_relation = 'RO:0002200'
                    qtl_to_trait_edge_label = 'biolink:has_phenotype'
                    qtl_to_trait_edge_relation = 'RO:0002200'
                    trait_to_value_edge_label = 'biolink:has_quantitative_value'
                    trait_to_value_edge_relation = ''
                    value_to_number_edge_label = 'biolink:has_numeric_value'
                    value_to_number_edge_relation = 'OBI:0001938'
                    value_to_unit_edge_label = 'biolink:has_unit'
                    value_to_unit_edge_relation = 'RO:0002536'

                    # transform
                    for index, row in gaf_df.iterrows():
                        tax_id = row['Taxon'].split(':')[1]
                        ontology_id = row['Ontology_ID']
                        evidence = row['Evidence_Code']
                        object_type = row['DB_Object_Type']
                        db = row['DB']
                        label = row['DB_Object_Symbol']
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
                            if 'At' in gene_id:
                                gene_id = gene_id.upper()
                        elif tax_id == '4530':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa'
                        elif tax_id == '39947':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa japonica'
                        elif tax_id == '39946':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa indica'
                        elif tax_id == '1080340':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Oryza sativa japonica x indica'
                        elif tax_id == '4558':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Sorghum bicolor'
                        elif tax_id == '3694':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Populus trichocarpa'
                        elif tax_id == '4577':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays'
                        elif tax_id == '381124':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays mays'
                        elif tax_id == '112001':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays huehuetenangensis'
                        elif tax_id == '76912':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays parviglumis'
                        elif tax_id == '4579':
                            gene_id = row['DB_Object_ID']
                            org_name = 'Zea mays mexicana'
                        else:
                            print('YOU NEED TO ADD A NEW TAXON ' + str(tax_id))
                        provided_by = db + '-' + source_name
                        org_id = org_prefix + str(tax_id)
                        #if the data are from rice normalize to LOC gene IDs
                        if file_type == 'genetic':
                            if tax_id == '4530' or tax_id == '39947' or tax_id == '39946' or tax_id == '1080340':
                                if 'LOC' not in gene_id:
                                    gene_id = row['DB_Object_Name']
                                if 'LOC' not in gene_id:
                                    gene_id = row['DB_Object_Symbol']
                                if '-' in gene_id:
                                    gene_id = gene_id.split('-')[0]
                                if 'LOC' not in gene_id and gene_id in rice_gene_ids:
                                    gid = rice_gene_ids[gene_id]
                            #if the data are from corn normalize to GRMZM gene IDs or is it Zm IDs?
                            if tax_id == '381124' or tax_id == '4577' or tax_id == '112001' or tax_id == '76912' or tax_id == '4579':
                                if isinstance(gene_id, int):
                                    try:
                                        gene_id = row['DB_Object_Synonym'].split('|')[0]
                                        if 'GRMZM' not in gene_id:
                                            gene_id = row['DB_Object_Name']
                                            if 'GRMZM' not in gene_id:
                                                gene_id = row['DB_Object_Symbol']
                                    except AttributeError:
                                        gene_id = row['DB_Object_Symbol']
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
                        if isinstance(gene_id, str):
                            genes = [gene_id]
                        if isinstance(gene_id, int):
                             gi = str(gene_id)
                             genes = [gi]
                        if isinstance(gene_id, list):
                            genes = gene_id
                        #create germplasm and germplasm trait nodes
                        #print(data_file + 'this is data file')
                        if file_type == 'germplasm':
                            #print('germ')
                            trait_type = None
                            for g in genes:
                                if g not in seen_node:
                                    #print('A')
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[g,
                                                               label,
                                                               germplasm_node_type,
                                                               provided_by])
                                    seen_node[g] += 1
                            #pheno_id = ontology_id
                            #print(row['Annotation_Extension'])
                            #if there is no phenotype annotation, then skip
                            if isinstance(row['Annotation_Extension'], float):
                                pass
                            else:
                                if '=' in row['Annotation_Extension']:
                                    #print('B')
                                    pheno = row['Annotation_Extension'].split('=')[1].strip(')')
                                    if '_' in pheno:
                                        pheno = pheno.split('_')[0]
                                elif '(' in row['Annotation_Extension']:
                                    #print('C')
                                    q = ['Very_weak_all_plants_flat_','Weak_most_plants_nearly_flat_','Strong_no_lodging_','Intermediate_most_plants_moderately_lodged_','Moderately_strong_most_plants_leaning_','Short__<30cm_','Intermediate__30-59_cm_','Tall__>59_cm_']
                                    pheno = row['Annotation_Extension'].split('(')[1].strip(')')
                                    if pheno not in q:
                                        if '_' in pheno:
                                            pheno = pheno.split('_')[0]
                                else:
                                    print('new delimiter in phenotype data')
                            if ontology_id in compound_traits:
                                #print('here')
                                if 'EAR_DIAMETER' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000433'
                                if 'EAR_HEIGHT' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000683'
                                if 'EAR_LENGTH' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000431'
                                if 'EAR_NUMBER' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000443'
                                if 'EAR_SHAPE' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000964'
                                if 'KERNEL_ROW_ARRANGEMENT' in row['Annotation_Extension']:
                                    ontology_id = 'TO:2000109'
                                if 'KERNEL_TYPE' in row['Annotation_Extension']:
                                    ontology_id = 'TO:0000575'
                                if 'ROOT_LODGING' in row['Annotation_Extension']:
                                    ontology_id = 'TO:2000158'
                                if 'STALK_LODGING' in row['Annotation_Extension']:
                                    ontology_id = 'TO:2000159'
                                #print(ontology_id)
                            z = ['(cm)','(count)','(mm)','(gm)']
                            if ontology_id in plant_trait_ids:
                                #print('D')
                                trait_type = 'categorical'
                                if ontology_id == 'TO:0000019':
                                    trait_type = 'categorical'
                                else:
                                    for y in z:
                                        if y in row['Annotation_Extension']:
                                            trait_type = 'numerical'
                            if ontology_id in plant_traits_num:
                                #print('E')
                                trait_type = 'numerical'
                            if ontology_id in ignore_traits:
                                trait_type = 'ignore'
                            if trait_type == None:
                                print('missed a trait type' + ' ' + ontology_id + row['Annotation_Extension'])
                            if trait_type == 'categorical':
                                try:
                                    #print(pheno)
                                    pheno_id = plant_trait_ids[ontology_id]['pheno'][pheno]
                                    pheno_label = ''
                                    if pheno_id not in seen_node:
                                        write_node_edge_item(fh=node,
                                                             header=self.node_header,
                                                             data=[pheno_id,
                                                                   pheno_label,
                                                                   trait_value_node_type,
                                                                   provided_by])
                                        seen_node[pheno_id] += 1
                                except KeyError:
                                    print('Phenotype error')
                                    print(ontology_id)
                                    print(pheno)
                                    print(pheno_label)
                            if trait_type == 'numerical':
                                r = str(row['DB_Object_ID'])+'|'+str(row['Ontology_ID'])+'|'+str(pheno)
                                if ontology_id not in seen_node:
                                    label = plant_traits_num[ontology_id]['label']
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[ontology_id,
                                                               label,
                                                               trait_value_node_type,
                                                               provided_by])
                                    seen_node[ontology_id] += 1
                                if r not in seen_node:
                                    label = plant_traits_num[ontology_id]['label'] + ' measurement'
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[r,
                                                               label,
                                                               value_node_type,
                                                               provided_by])
                                    seen_node[r] += 1
                        #create qtl nodes
                        elif file_type == 'qtl':
                            #print('qtl')
                            for g in genes:
                                if g not in seen_node:
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[g,
                                                               label,
                                                               qtl_node_type,
                                                               provided_by])
                                    seen_node[g] += 1
                        #create gene node
                        else:
                            for g in genes:
                                if 'AGI_LocusCode' in g:
                                    g = g.split(':')[1]
                                if g not in seen_node:
                                    g = str(g)
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[g,
                                                               label,
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
                                label = ''
                            elif row['Aspect'] == 'G':
                                node_type = growth_stage_node_type
                                label = ''
                            elif row['Aspect'] == 'C':
                                node_type = cellular_component_node_type
                                label = ''
                            elif row['Aspect'] == 'F':
                                node_type = molecular_function_node_type
                                label = ''
                            elif row['Aspect'] == 'P':
                                node_type = process_node_type
                                label = ''
                            else:
                                print('Error. New Aspect.')
                                print(row['Aspect'])
                            #print('I')
                            write_node_edge_item(fh=node,
                                                 header=self.node_header,
                                                 data=[ontology_id,
                                                       label,
                                                       node_type,
                                                       provided_by])
                            seen_node[ontology_id] += 1

                        #create additional nodes for orthologs
                        if 'ortholog' in data_file:
                            #print('ortholog')
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
                            for o in orth:
                                if o not in seen_node:
                                    gene_name = 'none'
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[o,
                                                               label,
                                                               gene_node_type,
                                                               provided_by])
                                    seen_node[o] += 1
                    # Write Edge
                        #for the germplasm files
                        if file_type == 'germplasm':
                            if str(g)+str(org_id) not in seen_edge:
                                #print('F')
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[g,
                                                            germplasm_to_org_edge_label,
                                                            org_id,
                                                            germplasm_to_org_edge_relation,
                                                            provided_by])
                                seen_edge[str(g)+str(org_id)] += 1
                            if trait_type == 'categorical':
                                #print(ontology_id + ' ' + pheno)
                                if str(g)+str(pheno_id) not in seen_edge:
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[g,
                                                                germplasm_to_trait_edge_label,
                                                                pheno_id,
                                                                germplasm_to_trait_edge_relation,
                                                                provided_by])
                                    seen_edge[str(g)+str(pheno_id)] += 1
                            if trait_type == 'numerical':
                                #print('H')
                                if r not in seen_edge:
                                    unit = plant_traits_num[ontology_id]['unit']
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[g,
                                                                germplasm_to_trait_edge_label,
                                                                ontology_id,
                                                                germplasm_to_trait_edge_relation,
                                                                provided_by])
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[ontology_id,
                                                                trait_to_value_edge_label,
                                                                r,
                                                                trait_to_value_edge_relation,
                                                                provided_by])
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[r,
                                                                value_to_number_edge_label,
                                                                pheno,
                                                                value_to_number_edge_relation,
                                                                provided_by])
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[r,
                                                                value_to_unit_edge_label,
                                                                unit,
                                                                value_to_unit_edge_relation,
                                                                provided_by])
                                    seen_edge[r] += 1

                        #for the ortholog files
                        elif 'ortholog' in data_file:
                            for g in genes:
                                for o in orth:
                                    if str(g)+str(o) not in seen_edge:
                                        write_node_edge_item(fh=edge,
                                                                header=self.edge_header,
                                                                data=[g,
                                                                    gene_to_orth_edge_label,
                                                                    o,
                                                                    gene_to_orth_edge_relation,
                                                                    provided_by])
                                        seen_edge[str(g)+str(o)] += 1
                        #for the qtl files
                        elif file_type == 'qtl':
                            if g + org_id not in seen_edge:
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[gene_id,
                                                            qtl_to_org_edge_label,
                                                            org_id,
                                                            qtl_to_org_edge_relation,
                                                            provided_by])
                                seen_edge[gene_id + org_id] += 1
                            if g + ontology_id not in seen_edge:
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[gene_id,
                                                            qtl_to_trait_edge_label,
                                                            ontology_id,
                                                            qtl_to_trait_edge_relation,
                                                            provided_by])
                                seen_edge[g + ontology_id] += 1

                        # gene to org edge
                        else:
#                            if isinstance(gene_id, str):
#                                genes = [gene_id]
#                            if isinstance(gene_id, int):
#                                gi = str(gene_id)
#                                genes = [gene_id]
                            for g in genes:
                                if str(g)+str(org_id) not in seen_edge:
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[g,
                                                                gene_to_org_edge_label,
                                                                org_id,
                                                                gene_to_org_edge_relation,
                                                                provided_by])
                                    seen_edge[str(g)+str(org_id)] += 1

                            # gene to all other nodes edges
                            if str(g)+ontology_id not in seen_edge:
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
                                seen_edge[str(g)+ontology_id] += 1
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
                                seen_edge[ontology_id+org_id] += 1
                # Files write ends

        """
        Implement ROBOT 
        """
        # Convert OWL to JSON for Ontologies
        #convert_to_json(self.input_base_dir, 'PO')
        #convert_to_json(self.input_base_dir, 'TO')
        # Extract the 'plant' tree from NCBITaxon and convert to JSON
        '''
        NCBITaxon_33090 = viridiplantae
        (Source = http://www.ontobee.org/ontology/NCBITaxon?iri=http://purl.obolibrary.org/obo/NCBITaxon_33090)
        '''
        subset_ontology_needed = 'NCBITaxon'
        extract_convert_to_json(self.input_base_dir, subset_ontology_needed, 'NCBITaxon:33090','TOP')#maybe should be TOP
        #Can I chop out branches of the hierarchy? NCBITaxon_144314
