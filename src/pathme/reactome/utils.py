# -*- coding: utf-8 -*-

"""This module has utilities method for parsing, handling WikiPathways RDF and data."""

import logging
import tarfile
from typing import List, Tuple

from bio2bel_chebi import Manager as ChebiManager
from bio2bel_hgnc import Manager as HgncManager
from bio2bel_hgnc.models import HumanGene
from pybel.dsl import protein
from ..constants import ENSEMBL, HGNC, UNIPROT, UNKNOWN, ZFIN
from ..utils import parse_id_uri, ebel

logger = logging.getLogger(__name__)

"""Download utilities"""


def get_hgnc_node_info(gene: HumanGene) -> Tuple[str, str, str]:
    """Return HGNC identifier, symbol and namespace from HGNC entry.
    :param bio2bel_hgnc.manager.models.HGNC gene:
    """
    return str(gene.identifier), gene.symbol, HGNC


def get_valid_node_parameters(
    node,
    hgnc_manager: HgncManager,
    chebi_manager: ChebiManager,
) -> Tuple[str, str, str]:
    """Get valid node parameters."""
    namespace = None

    if 'uri_id' in node:
        _, _, namespace, identifier = parse_id_uri(node['uri_id'])

    # Look up in HGNC Manager the HGNC Symbol for a given UniProt or ENSEMBL identifier.
    # If not matches anything, leave it as it is and give a warning.
    if namespace == 'uniprot':

        up_sql = "SELECT symbol FROM zfin WHERE uniprot.id = '{}'"
        zfin_results = ebel.client.command(up_sql.format(identifier))
        zfin_symbol = zfin_results[0].oRecordData['symbol']

        if not zfin_symbol:
            logger.debug('UniProt id: %s could not be converted to HGNC', identifier)
            namespace = UNIPROT

        # Multiple HGNC entries match the UniProt ID
        elif 1 < len(zfin_symbol):
            identifier = zfin_symbol
            namespace = 'hgnc_multiple_entry'

        else:
            identifier, name, namespace = identifier, zfin_symbol, ZFIN

    elif namespace == 'ensembl':
        ens_sql = "SELECT symbol FROM zfin WHERE ensembl.gene_id_short = '{}'"
        zfin_results = ebel.client.command(ens_sql.format(ensembl_id))
        zfin_symbol = zfin_results[0].oRecordData['symbol']

        if not zfin_symbol:
            logger.debug('ENSEMBL id: %s could not be converted to ZFIN', identifier)
            namespace = ENSEMBL

        else:
            identifier, name, namespace = identifier, zfin_symbol, ZFIN

    elif namespace == 'obo' and 'CHEBI' in identifier:
        namespace = 'chebi'

    elif 'uri_reactome_id' in node:
        namespace = 'reactome'
        identifier = node.get('reactome_id')

        if identifier is None:
            _, _, _, identifier = parse_id_uri(node['uri_reactome_id'])
            if '#' in identifier:
                identifier = str(identifier).split('#')[1]
        if 'Complex' not in identifier or 'SmallMolecule' not in identifier:
            logger.debug('Adding Reactome identifier for %s ', node['uri_reactome_id'])

    else:
        logger.debug('Not found HGNC Symbol neither Reactome id for %s ', node['uri_id'])

    if 'display_name' in node:
        name = node['display_name']
        if namespace == 'chebi' or namespace == 'CHEBI':
            if not chebi_manager.get_chemical_by_chebi_name(node['display_name']):
                identifier = identifier.replace('CHEBI:', '')
                chem = chebi_manager.get_chemical_by_chebi_id(identifier)

                # In case chebi id is outdated use the identifier as the name
                if chem:
                    name = chem.safe_name
                else:
                    name = identifier

    else:
        if 'name' in node:
            name = node['name']
            if isinstance(name, set):
                name = list(name)[0]
        else:
            name = UNKNOWN

    return identifier, name, namespace


def process_multiple_proteins(hgnc_entries: List) -> List:
    """Create multiple nodes when UniProt identifer refers to multiple HGNC symbols.
    :param hgnc_entries: Results from query
    :return: List of Protein BEL nodes
    """
    return [
        protein(namespace='ZFIN', name=hgnc_entry.symbol, identifier=hgnc_entry.id)
        for hgnc_entry in hgnc_entries
    ]


def untar_file(file_path: str, export_folder: str) -> None:
    """Unzip file into a destination folder.
    :param file_path: name of the file
    :param export_folder: name of the file
    """
    tar_ref = tarfile.open(file_path, 'r:bz2')
    tar_ref.extractall(export_folder)
    tar_ref.close()