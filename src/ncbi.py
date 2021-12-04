#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain FASTA with gene sequences from list of NCBI accession ids
Dependencies: ncbi-acc-download
"""

import os
import subprocess
from typing import OrderedDict

from Bio import SeqIO, SeqFeature, SeqRecord


def terminalExecute(command_str: str,
                    suppress_shell_output=False,
                    work_dir: str = None,
                    return_output=False) -> subprocess.STDOUT:
    """
    Execute given command in terminal through Python
    """
    if suppress_shell_output:
        stdout = subprocess.DEVNULL
    else:
        stdout = None
    output = subprocess.run(
        command_str, shell=True,
        cwd=work_dir, capture_output=return_output,
        stdout=stdout
        )
    return output

class GeneBank2Fasta():

    def __init__(self, entry_ids: list, data_dir: str = None) -> None:
        """
        Tools to download selected GenBank feactures from NCBI records.

        @Arguments:
        entry_ids: list containing valid GenBank entry ids
        output_dir: path to directory where genbank files will be downloaded to
        """
        self.entry_ids = entry_ids
        if data_dir is None:
            self._gbk_dir = os.path.join(os.getcwd(), 'gbk_data')
        else:
            self._gbk_dir = os.path.abspath(data_dir)
        if not os.path.isdir(self._gbk_dir):
            os.mkdir(self._gbk_dir)
        if not os.listdir(self._gbk_dir):
            self._downloadGBKfromNCBI()
    
    @classmethod
    def fromEntryIDs(entry_ids: list):
        pass
    @classmethod
    def fromGBKdirectory(gbk_dir: str):
        pass
    
    def _downloadGBKfromNCBI(self) -> None: 
        """
        Download GenBank files from NCBI from given list of entry IDs
        """
        print('Downloading GenBank files')
        for n, entry_id in enumerate(self.entry_ids): 
            print(f'Downloading entry: {entry_id} ({n + 1} / {len(self.entry_ids)})', end='\r')
            outfasta = os.path.join(self._gbk_dir, f'{entry_id}.gbk')
            cmd_str = f'ncbi-acc-download -o {outfasta} {entry_id}'
            terminalExecute(cmd_str)
    
    def _getGBKobject(self, entry_id: str) -> list:
        gbk = GBK(os.path.join(self._gbk_dir, entry_id))
        return gbk

    def _getCDSMatchingKeywords(self, entry_id: str, feature_keywords: dict,
                                case_insensitive: bool = True) -> SeqFeature:
        """
        Extract cds record matching keywords from gbk file.
        """
        gbk = self._getGBKobject(entry_id)
        try:
            return (gbk.gbk_object, gbk.cds.get_by_keywords(feature_keywords, case_insensitive))
        except Exception:
            return None

    def _writeProteinSequencesFromCDSrecords(self, records: dict, output_fasta: str = None) -> None: 
        """
        Write FASTA file from dict of gbk cds SeqFeatures

        records = {str: (SeqRecord, SeqFeature)}
        """
        with open(output_fasta, 'w') as file:
            for record_id, record in records.items():
                if record is not None:
                    record = record[1].qualifiers
                    ref_id = f'{record_id}_{record["protein_id"][0]}_{"_".join(record["product"][0].split())}'
                    file.write(f'>{ref_id}\n')
                    file.write(f'{record["translation"][0]}\n')

    def _writeNucleotideSequencesFromCDSrecords(self, records: dict, output_fasta: str = None) -> None: 
        """
        Write FASTA file from dict of gbk record qualifiers (Ordered dict)
        """
        with open(output_fasta, 'w') as file:
            for record_id, record in records.items():
                if record is not None:
                    gbk_obj, cds = record
                    q = cds.qualifiers
                    seq = str(cds.extract(gbk_obj).seq)
                    ref_id = f'{record_id}_{q["protein_id"][0]}_{"_".join(q["product"][0].split())}'
                    file.write(f'>{ref_id}\n')
                    file.write(f'{seq}\n')

    def writeSequencesInFasta(self, gene_keywords: dict,
                              output_fasta: str = None,
                              case_insensitive: bool = True,
                              sequence: str = 'protein') -> None:
        """
        Write FASTA from list of GenBank files and cds entry field keywords
        @argument:
        sequence: str, either 'protein' or 'nucleotide'
        """
        if output_fasta is None:
            output_fasta = os.path.join(os.getcwd(), 'sequences.fasta')
        records_dict = {
            gbk_file.split('.gbk')[0]: self._getCDSMatchingKeywords(
                os.path.join(self._gbk_dir, gbk_file),
                feature_keywords=gene_keywords,
                case_insensitive=case_insensitive
                ) for gbk_file in os.listdir(self._gbk_dir)
        }
        if 'protein' in sequence.lower():
            self._writeProteinSequencesFromCDSrecords(records_dict, output_fasta)
        else:
            self._writeNucleotideSequencesFromCDSrecords(records_dict, output_fasta)

class GBK():
    """
    Tools to parse GenBank files
    """
    def __init__(self, gbk_file: str) -> None:
        self._gbk = list(SeqIO.parse(gbk_file, 'genbank'))[0]
    
    @property
    def cds(self) -> list:
        return GBKfeatureList([f for f in self._gbk.features if 'cds' in f.type.lower()])

    @property
    def meta(self) -> OrderedDict:
        return self._gbk.features[0].qualifiers

    @property
    def gbk_object(self) -> SeqRecord:
        return self._gbk

class GBKfeatureList(list):
    """
    Allow filtering cds feature objects from gbk by keywords
    """
    def __init__(self, seq=None):
        super(self.__class__, self).__init__(seq)

    def __getslice__(self, start, stop):
        return self.__class__(super(self.__class__, self).__getslice__(start, stop))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.__class__(super(self.__class__, self).__getitem__(key))
        else:
            return super(self.__class__, self).__getitem__(key)
    
    @staticmethod
    def _text_contains_keywords(text: str, keywords: list,
                                case_insensitive: bool = True) -> bool: 
        if case_insensitive:
            return all([key.lower() in text.lower() for key in keywords])
        else:
            return all([key in text for key in keywords])
    
    def _cds_matched(self, cds: SeqFeature, feature_keywords: list,
                     case_insensitive: bool = True) -> bool:
        cds = cds.qualifiers
        return all(
            [(field in cds.keys() and
              self._text_contains_keywords(cds[field][0], keywords, case_insensitive))
            for field, keywords in feature_keywords.items()]
            )

    def get_by_keywords(self, keywords: dict, case_insensitive: bool = True) -> SeqFeature:
        """
        Extract cds record matching keywords.
        @Arguments:
        keywords is a dictionary in which keys correspond to  cds fields 
        and values to keywords to find in each field. For instace,
        keywords = {
            'gene': ['ureC'],
            'product': ['urease', 'alpha'] 
        }
        case_insensitive: whether or not to care for case when matching keywords 
        """
        try:
            return [
                cds for cds in self if self._cds_matched(cds, keywords, case_insensitive)
                ][0]
        except Exception as e:
            raise ValueError(f'Feature not found for given keyword(s). Exception: {e}')

    def get_by_gene_id(self, gene_id: str, case_insensitive: bool = True) -> SeqFeature:
        return self.get_by_keywords({'gene': [gene_id]}, case_insensitive)

    def get_by_ec_number(self, ec: str, case_insensitive: bool = True) -> SeqFeature:
        return self.get_by_keywords({'ec_number': [ec]}, case_insensitive)
