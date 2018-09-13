# -*- coding: utf-8 -*-

"""This module has utilities method for parsing and handling KEGG KGML files."""

import os
import pandas as pd
import requests
import tqdm

from bio2bel_kegg import Manager as KeggManager
from pathme.wikipathways.utils import get_files_in_folder
from pathme.kegg.kegg_xml_parser import import_xml_etree, populate_graph, get_node_types, get_edge_types, get_reaction_edge_types, get_xml_types
from pathme.kegg.convert_to_bel import get_bel_types
from ..constants import DATA_DIR, KEGG, KEGG_KGML_URL, KEGG_STATS_COLUMN_NAMES


def get_kegg_pathway_ids(connection=None):
    """Return a list of all pathway identifiers stored in the KEGG database.

    :param Optional[str] connection: connection to the database
    :returns: list of all kegg_pathway_ids
    :rtype: list
    """
    kegg_manager = KeggManager(connection=connection)
    kegg_pathways_ids = [
        pathway.resource_id.replace('path:', '')
        for pathway in kegg_manager.get_all_pathways()
    ]

    if not kegg_pathways_ids:
        raise EnvironmentError('Your database is empty. Please run python3 -m bio2bel_kegg populate')

    return kegg_pathways_ids


def download_kgml_files(kegg_pathway_ids):
    """Downloads KEGG KGML files by querying the KEGG API.

    :param list kegg_pathway_ids: list of kegg ids
    """
    for kegg_id in tqdm.tqdm(kegg_pathway_ids, desc='Downloading KEGG files'):
        request = requests.get(KEGG_KGML_URL.format(kegg_id))
        with open(os.path.join(DATA_DIR, KEGG, '{}.xml'.format(kegg_id)), 'w+') as file:
            file.write(request.text)
            file.close()


def get_kegg_statistics(path):
    """Parse a folder and get KEGG statistics.

    :param graph: path
    :param str path: path to folder containing XML files
    :return: relation edge types in XML
    :rtype: pandas.DataFrame
    """
    df = pd.DataFrame()
    export_file_name = 'KEGG_pathway_stats_unflattened.csv'

    # Get list of all files in folder
    files = get_files_in_folder(path)

    for file_name in tqdm.tqdm(files, desc='Parsing KGML files and BEL graphs for entities and relation stats'):
        pathway_names = []
        file_path = os.path.join(path, file_name)
        tree = import_xml_etree(file_path)
        root = tree.getroot()
        pathway_names.append(root.attrib['title'])

        # Get dictionary of all entity and interaction types in XML
        xml_statistics_dict = get_xml_types(tree)

        # Get dictionary of all node and edge types in BEL Graph
        bel_statistics_dict = get_bel_types(file_path, flatten=False)

        # Get dictionary with all XML and BEL graph stats
        xml_statistics_dict.update(bel_statistics_dict)

        # Update dictionary of all XML and BEL graph stats with corresponding column names
        all_kegg_statistics = {
            KEGG_STATS_COLUMN_NAMES[key]: value
            for key, value in xml_statistics_dict.items()
        }

        # Add pathway stat rows to dataframe
        pathway_data = pd.DataFrame(all_kegg_statistics,
                                    index=pathway_names,
                                    columns=KEGG_STATS_COLUMN_NAMES.values(),
                                    dtype=int)
        df = df.append(pathway_data)

    df.to_csv(export_file_name, sep='\t')
    return df
