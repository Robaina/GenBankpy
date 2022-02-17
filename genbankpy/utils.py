#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Obtain FASTA with gene sequences from list of NCBI accession ids
Dependencies: ncbi-acc-download
"""

import os
import subprocess

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