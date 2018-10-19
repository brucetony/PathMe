# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling wikipathways RDF and data."""

import logging
import os
import tarfile
from typing import List

from pybel.dsl import protein

from pathme.utils import parse_id_uri
from ..constants import DATA_DIR, WIKIPATHWAYS, HGNC, UNKNOWN, UNIPROT, ENSEMBL

WIKIPATHWAYS_DIR = os.path.join(DATA_DIR, WIKIPATHWAYS)

log = logging.getLogger(__name__)

"""Download utilities"""


def get_hgnc_node_info(gene):
    """Return HGNC identifier, symbol and namespace from HGNC entry.

    :param bio2bel_hgnc.manager.models.HGNC gene:
    :rtype: tuple[str,str,str]
    """
    return gene.identifier, gene.symbol, HGNC

def get_valid_node_parameters(node, hgnc_manager):

    if 'display_name' in node:
        name = node['display_name']
    else:
        if 'name' in node:
            name = node['name']
            if isinstance(name, set):
                name = list(name)[0]
        else:
            name = UNKNOWN

    namespace = node.get('db')
    identifier = node.get('identifier')

    if namespace is None:
        _, _, namespace, identifier = parse_id_uri(node['uri_id'])

    # Look up in HGNC Manager the HGNC Symbol for a given UniProt or ENSEMBL identifier.
    # If not matches anything, leave it as it is and give a warning.
    if namespace == 'uniprot':

        hgnc_entry = hgnc_manager.get_gene_by_uniprot_id(identifier)

        if not hgnc_entry:
            log.debug('UniProt id: %s could not be converted to HGNC', identifier)
            namespace = UNIPROT

        # Multiple HGNC entries match the UniProt ID
        elif 1 < len(hgnc_entry):
            identifier = hgnc_entry
            namespace = 'hgnc_multiple_entry'

        else:
            identifier, name, namespace = get_hgnc_node_info(hgnc_entry[0])

    elif namespace == 'ensembl':
        hgnc_entry = hgnc_manager.get_gene_by_ensembl_id(identifier)

        if not hgnc_entry:
            log.debug('ENSEMBL id: %s could not be converted to HGNC', identifier)
            namespace = ENSEMBL

        else:
            identifier, name, namespace = get_hgnc_node_info(hgnc_entry)

    elif (namespace == 'obo' and 'CHEBI' in identifier) or namespace == 'chebi':
        namespace = 'chebi'

    elif 'uri_reactome_id' in node:
        namespace = 'reactome'
        identifier = node.get('reactome_id')

        if identifier is None:
            _, _, _, identifier = parse_id_uri(node['uri_reactome_id'])
            if '#' in identifier:
                identifier = str(identifier).split('#')[1]
        if 'Complex' not in identifier or 'SmallMolecule' not in identifier:
            log.debug('Adding Reactome identifier for %s ', node['uri_reactome_id'])

    else:
        log.debug('Not found HGNC Symbol neither Reactome id for %s ', node['uri_id'])

    return identifier, name, namespace


def process_multiple_proteins(hgnc_entries: List) -> List:
    """Create multiple nodes when UniProt identifer refers to multiple HGNC symbols.

    :param hgnc_entries: Results from query
    :return: List of Protein BEL nodes
    """
    protein_group = list()

    for hgnc_entry in hgnc_entries:
        protein_group.append(
            protein(namespace='HGNC', name=hgnc_entry.symbol, identifier=hgnc_entry.id)
        )

    return protein_group

def untar_file(file_path, export_folder):
    """Unzip file into a destination folder.

    :param str file_path: name of the file
    :param str export_folder: name of the file
    """
    tar_ref = tarfile.open(file_path, 'r:bz2')
    tar_ref.extractall(export_folder)
    tar_ref.close()
