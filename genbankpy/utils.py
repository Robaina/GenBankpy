#!/usr/bin/env python
# -*- coding: utf-8 -*-



from __future__ import annotations
import os
import random
import string
from zipfile import ZipFile
from pathlib import Path


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
                    work_dir: Path = None,
                    suppress_shell_output=False) -> None:
    """
    Execute given command in terminal through Python
    """
    current_dir = Path.cwd()
    if work_dir is None:
        work_dir = current_dir
    if suppress_shell_output:
        output_str = " > /dev/null 2>&1"
    else:
        output_str = ""
    os.chdir(work_dir)
    os.system(command_str + output_str)
    os.chdir(current_dir)

def unZipDirectory(input_zip: Path, output_dir: Path  = None) -> None:
    """
    Unzip file to specified directory
    """
    if output_dir is None:
        output_dir = Path(input_zip.parent)
    with ZipFile(input_zip, 'r') as zipObj:
        zipObj.extractall(path=output_dir, members=None)

def contains_substring(substring, strings: list) -> bool:
    return any(substring in string for string in strings)

def unZipFile(input_file: Path):
    """
    Unzip gz file and remove .gz file
    """
    cmd_str = f"gzip -d {input_file}"
    terminalExecute(cmd_str)

def mergeMultiRecordGBK(input_file: Path, output_file: Path) -> None:
    """
    Merge multi-record GBK into a single GBK
    via: https://github.com/kblin/merge-gbk-records
    """
    cmd_str = (
        f"merge-gbk-records {input_file} > {output_file}"
    )
    terminalExecute(cmd_str)