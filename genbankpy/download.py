#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download genome data from NCBI using NCBI's datasets CLI
https://github.com/ncbi/datasets
conda install -c conda-forge ncbi-datasets-cli
"""

from __future__ import annotations
import shutil
import tempfile
from pathlib import Path
from genbankpy.utils import terminalExecute, getJSONlinesObject, unZipDirectory



class NCBIdownloader:

    def __init__(self, data_directory: Path = None) -> None:
        self._data_dir = data_directory

    def updateDataDirectory(self, data_directory: Path) -> None: 
        """
        Update data directory where to download files to
        """
        if not data_directory.exists():
            data_directory.mkdir(parents=True, exist_ok=True)
        self._data_dir = data_directory

    def fromAccessionIDs(self, ids: list[str], data_directory: Path = None) -> dict[Path]:
        """
        Download genomes from NCBI by accession ID
        """
        data_dir = data_directory if data_directory is not None else self._data_dir
        self.updateDataDirectory(data_dir)
        zipfilename = "ncbi_dataset.zip"
        ids = [f"'{aid}'" for aid in ids]
        cmd_str = (
            f"datasets download genome accession {','.join(ids)} "
            "--exclude-gff3 --exclude-protein --exclude-rna --exclude-seq --exclude-genomic-cds --include-gbff "
            f"--filename {zipfilename} --no-progressbar"
        )
        with tempfile.TemporaryDirectory() as tmpdirname:
            terminalExecute(cmd_str, work_dir=tmpdirname)
            self._getDownloadedFiles(Path(tmpdirname), zipfilename)

    def fromSpecies(self, species: list[str], data_directory: Path = None) -> dict[Path]:
        """
        Download genomes from NCBI by species name
        NOTE: seems like only one taxon at a time is possible, thus implement iterative method
        """
        data_dir = data_directory if data_directory is not None else self._data_dir
        self.updateDataDirectory(data_dir)
        zipfilename = "ncbi_dataset.zip"
        species = [f"'{sp}'" for sp in species]
        cmd_str = (
            f"datasets download genome taxon {','.join(species)} "
            f"--exclude-gff3 --exclude-protein --exclude-rna --exclude-seq --exclude-genomic-cds --include-gbff "
            f"--reference --annotated --assembly-level 'complete_genome' --assembly-source 'refseq' "
            f"--filename {zipfilename} --no-progressbar"
        )
        with tempfile.TemporaryDirectory() as tmpdirname:
            terminalExecute(cmd_str, work_dir=tmpdirname)
            self._getDownloadedFiles(Path(tmpdirname), zipfilename)

    def _getDownloadedFiles(self, download_dir: Path, zipfilename: str) -> None:
        unzip_dir = download_dir / zipfilename.split('.')[0]
        unZipDirectory(input_zip= download_dir / zipfilename, output_dir=download_dir)
        self.writeMetadataTable(data_directory=Path(unzip_dir / "data"), output_file=self._data_dir / "metadata.tsv")
        gbk_file_dict = self.getDownloadedGBKfiles(data_directory=Path(unzip_dir / "data"))
        for filename, file in gbk_file_dict.items():
            shutil.move(file, self._data_dir / filename)

    def getDownloadedGBKfiles(self, data_directory: Path = None) -> dict[Path]:
        """
        Retrieve GBK files from downloaded data directory
        """
        data_dir = data_directory if data_directory is not None else self._data_dir
        gbk_dirs = [d for d in data_dir.iterdir() if d.is_dir()]
        gbk_dirs.sort()
        return {
            f"{gbk_dir.stem}{file.suffix}": file
            for gbk_dir in gbk_dirs
            for file in gbk_dir.iterdir()
            if "gb" in file.suffix
        }

    def _getMetadata(self, data_directory: Path) -> list[dict]:
        """
        Get metadata of downloaded file as JSON
        """
        data_dir = data_directory if data_directory is not None else self._data_dir
        json_file = [
            file for file in data_dir.iterdir()
            if "data_report.jsonl" in file.as_posix()
            ].pop()
        result = getJSONlinesObject(json_file)
        return result

    def writeMetadataTable(self, data_directory: Path, output_file: Path = None) -> None:
        """
        Write download metadata to a tsv table
        """
        data_dir = data_directory if data_directory is not None else self._data_dir
        if output_file is None:
            output_file = data_dir / "metadata.tsv"
        meta_dicts = self._getMetadata(data_dir)
        with open(output_file, "w+") as meta:
            lines = ["accession\tname\tlevel\tstatus\torganism\ttaxid\trelease\n"]
            for meta_dict in meta_dicts:
                line = (
                    f"{meta_dict['assemblyInfo']['assemblyAccession']}\t"
                    f"{meta_dict['assemblyInfo']['assemblyName']}\t"
                    f"{meta_dict['assemblyInfo']['assemblyLevel']}\t"
                    f"{meta_dict['assemblyInfo']['assemblyStatus']}\t"
                    f"{meta_dict['organismName']}\t"
                    f"{meta_dict['taxId']}\t"
                    f"{meta_dict['annotationInfo']['name']}, {meta_dict['annotationInfo']['releaseDate']}\t"
                )
                lines.append(line + "\n")
            meta.writelines(lines)