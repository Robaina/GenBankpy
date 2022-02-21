#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain FASTA with gene sequences from list of NCBI accession ids
Dependencies: ncbi-acc-download
"""

import os
import gzip
import shutil
import random
import string
import subprocess


class TemporaryFilePath:
    """
    Custom context manager to create a temporary file
    which is removed when exiting context manager
    """
    def __init__(self,
                 work_dir: str = None,
                 extension: str = None,
                 create_file: bool = False):
        self.work_dir = work_dir or ''
        self.extension = extension or ''
        self.create_file = create_file

    def __enter__(self):
        temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
        self.file_path = os.path.join(
            self.work_dir, f'temp_{temp_id}{self.extension}'
            )
        if self.create_file:
            os.mkdir(self.file_path)
        return self.file_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.file_path):
            os.remove(self.file_path)


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

def contains_substring(substring, strings: list) -> bool:
    return any(substring in string for string in strings)

def fullPathListDir(dir: str) -> list:
    """
    Return full path of files in provided directory
    """
    return [os.path.join(dir, file) for file in os.listdir(dir)]

def setDefaultOutputPath(input_path: str, tag: str = None,
                         extension: str = None,
                         only_filename: bool = False,
                         only_dirname: bool = False) -> str:
    """
    Get default path to output file or directory
    """
    basename = os.path.basename(input_path)
    dirname = os.path.abspath(os.path.dirname(input_path))
    fname, ext = os.path.splitext(basename)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    default_file = f'{fname}{tag}{extension}'
    if only_filename:
        return default_file
    if only_dirname:
        return os.path.abspath(dirname)
    else:
        return os.path.abspath(os.path.join(dirname, default_file))

def unZipFile(input_file: str):
    """
    Unzip gz file and remove .gz file
    """
    cmd_str = f"gzip -d {input_file}"
    terminalExecute(cmd_str)

def mergeMultiRecordGBK(input_file: str, output_file: str) -> None:
    """
    Merge multi-record GBK into a single GBK
    via: https://github.com/kblin/merge-gbk-records
    """
    cmd_str = (
        f"merge-gbk-records {input_file} > {output_file}"
    )
    terminalExecute(cmd_str)