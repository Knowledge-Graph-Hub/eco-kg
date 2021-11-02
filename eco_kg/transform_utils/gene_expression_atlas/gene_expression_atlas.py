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
Ingest plant annotations from Gene Expression Atlas



"""

class GeneExpressionAtlasTransform(Transform):

    def __init__(self, input_dir: str = None, output_dir: str = None) -> None:
        source_name = "Gene Expression Atlas"
        super().__init__(source_name, input_dir, output_dir)  # set some variables
        self.node_header = ['id', 'name', 'category', 'provided_by']
        self.edge_header = ['subject','predicate','object','relation','has_attribute','has_attribute_type','has_quantitative_value','has_unit','has_qualitative_value','provided_by']

    def run(self, data_files: List[str] = None):
        source_name = 'Gene Expression Atlas'
        """Method is called and performs needed transformations to process the 
        plant data (Gene Expression Atlas)
        
        Args:
            Gene Expression Atlas: curated contents of Gene Expression Atlas []
        """
        #get rice gene ids
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
        #get corn gene ids
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

        data_files = ['drought_salt_gene_expression.txt']

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
                #open each file and load data into a pandas dataframe
                with open(input_file, 'r') as f:
                    header_row = ['Gene','Species','Experiment accession','Comparison','PECOid','log_2 fold change','Adjusted p-value','t-statistic','ROid','developmental stage','POid LifeStage','anatomical part','POid Anatomy','cultivar']
                    df = pd.read_csv(input_file, sep='\t', comment='!', names=header_row, low_memory=False)

                    # Nodes
                    env_node_type = "biolink:EnvironmentalExposure"
                    gene_node_type = "biolink:GenomicEntity"
                    
                    # Edges
                    increases_edge_label = "biolink:increases_expression_of"
                    increases_edge_relation = "RO:0003003"
                    decreases_edge_label = "biolink:decreases_expression_of"
                    decreases_edge_relation = "RO:0003002"

                    # iterate over each dataframe to do the transform
                    mystery_genes = ['Zm00001d027187','ENSRNA049469775','Zm00001d001633','Zm00001d039240','Zm00001d000237',
                    'Zm00001d013903','Zm00001d000909','ZeamMp108','Zm00001d016401','Zm00001d008017','Zm00001d001694',
                    'Zm00001d022886','Zm00001d000942','Zm00001d000840','Zm00001d001311']
                    for index, row in df.iterrows():
                        gene_id = row['Gene']
                        env_id = row['PECOid']
                        predicate = row['ROid']
                        taxon = row['Species']
                        provided_by = row['Experiment accession'] + '-Gene Expression Atlas'
                        #normalize rice gene ids
                        if taxon == 'oryza sativa':
                            gene_id = rice_gene_ids[gene_id]
                        #normalize corn gene ids
                        if taxon == 'zea mays':
                            if gene_id in mystery_genes:
                                pass
                            else:
                                gene_id = corn_gene_ids[gene_id]
                        #normalize sorghum gene ids
                        if taxon == 'sorghum bicolor':
                            gene_id = re.sub('SORBI_3','Sobic.',gene_id)
                        #create the gene node
                        if type(gene_id) == list:
                            for g in gene_id:
                                if g not in seen_node:
                                    label = ''
                                    write_node_edge_item(fh=node,
                                                         header=self.node_header,
                                                         data=[g,
                                                               label,
                                                               gene_node_type,
                                                               provided_by])
                                    seen_node[g] += 1
                        else:
                            if gene_id not in seen_node:
                                label = ''
                                write_node_edge_item(fh=node,
                                                     header=self.node_header,
                                                     data=[gene_id,
                                                           label,
                                                           gene_node_type,
                                                           provided_by])
                                seen_node[gene_id] += 1
                        #create the exposure node
                        if env_id not in seen_node:
                            label = ''
                            write_node_edge_item(fh=node,
                                                 header=self.node_header,
                                                 data=[env_id,
                                                       label,
                                                       env_node_type,
                                                       provided_by])
                            seen_node[env_id] += 1
                        #create exposure to gene edge
                        if type(gene_id) == list:
                            for g in gene_id: 
                                if env_id + g not in seen_edge:
                                    if predicate == 'biolink:increases_expression_of':
                                        relation = 'RO:0003003'
                                    if predicate == 'biolink:decreases_expression_of':
                                        relation = 'RO:0003002'
                                    has_attribute = ''
                                    has_attribute_type = ''
                                    has_quantitative_value = ''
                                    has_unit = ''
                                    has_qualitative_value = ''
                                    write_node_edge_item(fh=edge,
                                                            header=self.edge_header,
                                                            data=[env_id,
                                                                predicate,
                                                                g,
                                                                relation,
                                                                has_attribute,
                                                                has_attribute_type,
                                                                has_quantitative_value,
                                                                has_unit,
                                                                has_qualitative_value,
                                                                provided_by])
                                    seen_edge[env_id + g] += 1
                        else:
                            if env_id + gene_id not in seen_edge:
                                if predicate == 'biolink:increases_expression_of':
                                    relation = 'RO:0003003'
                                if predicate == 'biolink:decreases_expression_of':
                                    relation = 'RO:0003002'
                                has_attribute = ''
                                has_attribute_type = ''
                                has_quantitative_value = ''
                                has_unit = ''
                                has_qualitative_value = ''
                                write_node_edge_item(fh=edge,
                                                        header=self.edge_header,
                                                        data=[env_id,
                                                            predicate,
                                                            gene_id,
                                                            relation,
                                                            has_attribute,
                                                            has_attribute_type,
                                                            has_quantitative_value,
                                                            has_unit,
                                                            has_qualitative_value,
                                                            provided_by])
                                seen_edge[str(env_id) + str(gene_id)] += 1