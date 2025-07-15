#!/usr/bin/env python3
"""
UniProt Swiss-Prot XML Parser

This script parses UniProt Swiss-Prot XML files and extracts protein information
including accessions, names, sequences, organisms, and other annotations.

Usage:
    python uniprot_parser.py --input uniprot_sprot.xml --output proteins.csv
    python uniprot_parser.py --input uniprot_sprot.xml --output proteins.json --format json
    python uniprot_parser.py --input uniprot_sprot.xml --limit 100 --output sample.csv
"""

import xml.etree.ElementTree as ET
import csv
import json
import argparse
import sys
from typing import Dict, List, Any, Optional
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# UniProt namespace
UNIPROT_NS = "https://uniprot.org/uniprot"


class UniProtParser:
    """Parser for UniProt Swiss-Prot XML files."""
    
    def __init__(self, xml_file: str):
        """Initialize parser with XML file path."""
        self.xml_file = xml_file
        self.namespace = {'uniprot': UNIPROT_NS}
        
    def parse_entry(self, entry_elem) -> Dict[str, Any]:
        """Parse a single UniProt entry."""
        protein_data = {
            'accession': '',
            'name': '',
            'full_name': '',
            'organism': '',
            'organism_common': '',
            'taxonomy_id': '',
            'taxonomy_lineage': [],
            'sequence': '',
            'sequence_length': 0,
            'gene_name': '',
            'function': '',
            'similarity_info': [],
            'protein_families': [],
            'protein_domains': [],
            'protein_domain_descriptions': [],
            'keywords': [],
            'references': [],
            'features': [],
            'db_references': [],
            'protein_existence': '',
            'created': '',
            'modified': '',
            'version': '',
            'go_terms': [],  # List of dicts: {'id': GO_ID, 'aspect': aspect, 'description': description}
            'go_ids': [],
            'go_C_descriptions': [],
            'go_P_descriptions': [],
            'go_F_descriptions': []
        }
        
        # Get entry attributes
        protein_data['created'] = entry_elem.get('created', '')
        protein_data['modified'] = entry_elem.get('modified', '')
        protein_data['version'] = entry_elem.get('version', '')
        
        # Parse accession
        accession_elem = entry_elem.find('uniprot:accession', self.namespace)
        if accession_elem is not None:
            protein_data['accession'] = accession_elem.text or ''
        
        # Parse name
        name_elem = entry_elem.find('uniprot:name', self.namespace)
        if name_elem is not None:
            protein_data['name'] = name_elem.text or ''
        
        # Parse protein information
        protein_elem = entry_elem.find('uniprot:protein', self.namespace)
        if protein_elem is not None:
            recommended_name = protein_elem.find('.//uniprot:recommendedName/uniprot:fullName', self.namespace)
            if recommended_name is not None:
                protein_data['full_name'] = recommended_name.text or ''
        
        # Parse gene information
        gene_elem = entry_elem.find('uniprot:gene', self.namespace)
        if gene_elem is not None:
            gene_name = gene_elem.find('uniprot:name', self.namespace)
            if gene_name is not None:
                protein_data['gene_name'] = gene_name.text or ''
        
        # Parse organism information
        organism_elem = entry_elem.find('uniprot:organism', self.namespace)
        if organism_elem is not None:
            scientific_name = organism_elem.find('uniprot:name[@type="scientific"]', self.namespace)
            if scientific_name is not None:
                protein_data['organism'] = scientific_name.text or ''
            
            common_name = organism_elem.find('uniprot:name[@type="common"]', self.namespace)
            if common_name is not None:
                protein_data['organism_common'] = common_name.text or ''
            
            # Get taxonomy ID
            taxonomy_ref = organism_elem.find('uniprot:dbReference[@type="NCBI Taxonomy"]', self.namespace)
            if taxonomy_ref is not None:
                protein_data['taxonomy_id'] = taxonomy_ref.get('id', '')
            
            # Parse taxonomy lineage
            lineage_elem = organism_elem.find('uniprot:lineage', self.namespace)
            if lineage_elem is not None:
                taxons = lineage_elem.findall('uniprot:taxon', self.namespace)
                protein_data['taxonomy_lineage'] = [taxon.text for taxon in taxons if taxon.text]
        
        # Parse sequence
        sequence_elem = entry_elem.find('uniprot:sequence', self.namespace)
        if sequence_elem is not None:
            protein_data['sequence'] = sequence_elem.text or ''
            protein_data['sequence_length'] = int(sequence_elem.get('length', 0))
        
        # Parse comments (function, similarity, etc.)
        comments = entry_elem.findall('uniprot:comment', self.namespace)
        for comment in comments:
            comment_type = comment.get('type', '')
            if comment_type == 'function':
                text_elem = comment.find('uniprot:text', self.namespace)
                if text_elem is not None:
                    protein_data['function'] = text_elem.text or ''
            elif comment_type == 'similarity':
                text_elem = comment.find('uniprot:text', self.namespace)
                if text_elem is not None:
                    similarity_data = {
                        'text': text_elem.text or '',
                        'evidence': text_elem.get('evidence', '')
                    }
                    protein_data['similarity_info'].append(similarity_data)
                    
                    # Extract protein families from similarity text
                    similarity_text = text_elem.text or ''
                    if 'family' in similarity_text.lower():
                        # Extract family information
                        if 'belongs to the' in similarity_text.lower():
                            family_start = similarity_text.lower().find('belongs to the') + 14
                            family_end = similarity_text.find('.', family_start)
                            if family_end == -1:
                                family_end = len(similarity_text)
                            family_name = similarity_text[family_start:family_end].strip()
                            if family_name:
                                protein_data['protein_families'].append(family_name)
        
        # Parse keywords
        keywords = entry_elem.findall('.//uniprot:keyword', self.namespace)
        protein_data['keywords'] = [kw.text for kw in keywords if kw.text]
        
        # Parse references
        references = entry_elem.findall('uniprot:reference', self.namespace)
        for ref in references:
            ref_data = {}
            citation = ref.find('uniprot:citation', self.namespace)
            if citation is not None:
                title_elem = citation.find('uniprot:title', self.namespace)
                if title_elem is not None:
                    ref_data['title'] = title_elem.text or ''
                
                authors = citation.findall('.//uniprot:person', self.namespace)
                ref_data['authors'] = [author.get('name', '') for author in authors]
            
            scope_elem = ref.find('uniprot:scope', self.namespace)
            if scope_elem is not None:
                ref_data['scope'] = scope_elem.text or ''
            
            protein_data['references'].append(ref_data)
        
        # Parse features
        features = entry_elem.findall('uniprot:feature', self.namespace)
        for feature in features:
            feature_data = {
                'type': feature.get('type', ''),
                'description': feature.get('description', ''),
                'id': feature.get('id', '')
            }
            
            location = feature.find('uniprot:location', self.namespace)
            if location is not None:
                begin = location.find('uniprot:begin', self.namespace)
                end = location.find('uniprot:end', self.namespace)
                if begin is not None:
                    feature_data['begin'] = begin.get('position', '')
                if end is not None:
                    feature_data['end'] = end.get('position', '')
            
            protein_data['features'].append(feature_data)
        
        # Parse database references
        db_refs = entry_elem.findall('uniprot:dbReference', self.namespace)
        for db_ref in db_refs:
            ref_data = {
                'type': db_ref.get('type', ''),
                'id': db_ref.get('id', '')
            }
            
            properties = db_ref.findall('uniprot:property', self.namespace)
            ref_data['properties'] = {prop.get('type', ''): prop.get('value', '') for prop in properties}
            
            protein_data['db_references'].append(ref_data)
            
            # Extract protein domains from domain databases
            db_type = db_ref.get('type', '')
            if db_type in ['Pfam', 'InterPro', 'SMART', 'SUPFAM', 'PROSITE', 'CDD']:
                domain_info = {
                    'database': db_type,
                    'id': db_ref.get('id', ''),
                    'properties': ref_data['properties']
                }
                protein_data['protein_domains'].append(domain_info)
                # Extract entry name property as description if present
                for prop in properties:
                    if prop.get('type', '') == 'entry name':
                        protein_data['protein_domain_descriptions'].append(prop.get('value', ''))
            
            # Extract GO terms and descriptions
            if db_type == 'GO':
                go_id = db_ref.get('id', '')
                go_term = ''
                go_aspect = ''
                go_desc = ''
                for prop in properties:
                    if prop.get('type', '') == 'term':
                        go_term = prop.get('value', '')
                        if ':' in go_term:
                            go_aspect, go_desc = go_term.split(':', 1)
                            go_aspect = go_aspect.strip()
                            go_desc = go_desc.strip()
                        else:
                            go_aspect = ''
                            go_desc = go_term
                        break
                protein_data['go_terms'].append({'id': go_id, 'aspect': go_aspect, 'description': go_desc})
                protein_data['go_ids'].append(go_id)
                if go_aspect == 'C':
                    protein_data['go_C_descriptions'].append(go_desc)
                elif go_aspect == 'P':
                    protein_data['go_P_descriptions'].append(go_desc)
                elif go_aspect == 'F':
                    protein_data['go_F_descriptions'].append(go_desc)
        
        # Parse protein existence
        existence_elem = entry_elem.find('uniprot:proteinExistence', self.namespace)
        if existence_elem is not None:
            protein_data['protein_existence'] = existence_elem.get('type', '')
        
        return protein_data
    
    def parse_file(self, limit: Optional[int] = None) -> List[Dict[str, Any]]:
        """Parse the entire XML file and return list of protein data."""
        logger.info(f"Starting to parse {self.xml_file}")
        
        proteins = []
        entry_count = 0
        
        # Use iterparse for memory-efficient parsing of large files
        context = ET.iterparse(self.xml_file, events=('end',))
        
        for event, elem in context:
            if elem.tag.endswith('entry'):
                try:
                    protein_data = self.parse_entry(elem)
                    proteins.append(protein_data)
                    entry_count += 1
                    
                    if entry_count % 1000 == 0:
                        logger.info(f"Processed {entry_count} entries...")
                    
                    if limit and entry_count >= limit:
                        logger.info(f"Reached limit of {limit} entries")
                        break
                        
                except Exception as e:
                    logger.error(f"Error parsing entry {entry_count + 1}: {e}")
                
                # Clear the element to free memory
                elem.clear()
        
        logger.info(f"Successfully parsed {len(proteins)} protein entries")
        return proteins


def save_to_csv(proteins: List[Dict[str, Any]], output_file: str):
    """Save protein data to CSV file."""
    if not proteins:
        logger.warning("No protein data to save.")
        return
    # Flatten nested fields for CSV
    flattened_proteins = []
    for protein in proteins:
        flat_protein = {
            'accession': protein['accession'],
            'name': protein['name'],
            'full_name': protein['full_name'],
            'organism': protein['organism'],
            'organism_common': protein['organism_common'],
            'taxonomy_id': protein['taxonomy_id'],
            'taxonomy_lineage': '; '.join(protein['taxonomy_lineage']),
            'sequence': protein['sequence'],
            'sequence_length': protein['sequence_length'],
            'gene_name': protein['gene_name'],
            'function': protein['function'],
            'protein_existence': protein['protein_existence'],
            'created': protein['created'],
            'modified': protein['modified'],
            'version': protein['version'],
            'go_ids': '; '.join(protein.get('go_ids', [])),
            'go_C_descriptions': '; '.join(protein.get('go_C_descriptions', [])),
            'go_P_descriptions': '; '.join(protein.get('go_P_descriptions', [])),
            'go_F_descriptions': '; '.join(protein.get('go_F_descriptions', [])),
            'protein_domains': '; '.join([f"{d['database']}:{d['id']}" for d in protein['protein_domains']]),
            'protein_domain_descriptions': '; '.join(protein.get('protein_domain_descriptions', [])),
        }
        flattened_proteins.append(flat_protein)
    # Write to CSV
    fieldnames = list(flattened_proteins[0].keys())
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in flattened_proteins:
            writer.writerow(row)
    logger.info(f"Saved {len(flattened_proteins)} proteins to {output_file}")


def save_to_json(proteins: List[Dict[str, Any]], output_file: str):
    """Save protein data to JSON file."""
    with open(output_file, 'w', encoding='utf-8') as jsonfile:
        json.dump(proteins, jsonfile, indent=2, ensure_ascii=False)
    
    logger.info(f"Saved {len(proteins)} proteins to {output_file}")


def main():
    """Main function to parse command line arguments and run the parser."""
    parser = argparse.ArgumentParser(
        description='Parse UniProt Swiss-Prot XML files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input UniProt XML file path'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output file path'
    )
    
    parser.add_argument(
        '--format', '-f',
        choices=['csv', 'json'],
        default='csv',
        help='Output format (default: csv)'
    )
    
    parser.add_argument(
        '--limit', '-l',
        type=int,
        help='Limit number of entries to parse (for testing)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check if input file exists
    if not Path(args.input).exists():
        logger.error(f"Input file {args.input} does not exist")
        sys.exit(1)
    
    try:
        # Parse the XML file
        parser = UniProtParser(args.input)
        proteins = parser.parse_file(limit=args.limit)
        
        if not proteins:
            logger.warning("No proteins found in the XML file")
            return
        
        # Save to output file
        if args.format == 'csv':
            save_to_csv(proteins, args.output)
        elif args.format == 'json':
            save_to_json(proteins, args.output)
        
        logger.info(f"Parsing completed successfully!")
        logger.info(f"Total proteins processed: {len(proteins)}")
        
    except Exception as e:
        logger.error(f"Error during parsing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 