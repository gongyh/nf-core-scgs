###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import subprocess
import sys
import gzip
import shutil

"""Helper functions for running external programs."""


def run(cmd):
    """Run executable program.

    Parameters
    ----------
    cmd : str
        Command to run.

    Returns
    -------
    boolean
        True if executed , else False.
    exception
        Program output if executed successfully, else Exception.
    """

    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
        )
        stdout, stderr = proc.communicate()
        return True, stdout, stderr
    except subprocess.CalledProcessError as e:
        return False, e


def is_executable(fpath):
    """Check if file is executable.

    This is a Python implementation of the linux
    command 'which'.

    Parameters
    ----------
    fpath : str
        Path to file.

    Returns
    -------
    boolean
        True if executable, else False.
    """

    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """Return path to program.

    This is a Python implementation of the linux
    command 'which'.

    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Parameters
    ----------
    program : str
        Name of executable for program.

    Returns
    -------
    str
        Path to executable, or None if it isn't on the path.
    """

    fpath, _fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_executable(exe_file):
                return exe_file

    return None


def check_on_path(program, exit_on_fail=True):
    """Check if program is on the system path.

    Parameters
    ----------
    program : str
        Name of executable for program.
    exit_on_fail : boolean
        Exit program with error code -1 if program in not on path.

    Returns
    -------
    boolean
        True if program is on path, else False.
    """

    if which(program):
        return True

    if exit_on_fail:
        print("%s is not on the system path." % program)
        sys.exit(1)

    return False


def check_dependencies(programs, exit_on_fail=True):
    """Check if all required programs are on the system path.

    Parameters
    ----------
    programs : iterable
        Names of executable programs.
    exit_on_fail : boolean
        Exit program with error code -1 if any program in not on path.

    Returns
    -------
    boolean
        True if all programs are on path, else False.
    """

    for program in programs:
        if not check_on_path(program, exit_on_fail):
            return False
    return True


def check_directory(directory_path):
    """Check if directory existed if not create

    Args:
        directory_path (string): Path of directory to validate
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def delete_files(file_path):
    """Delete files

    Args:
        file_path (str): files path
    """
    if os.path.exists(file_path):
        os.remove(file_path)


def delete_temp_dir(dir_path):
    """Delete temporary dir

    Args:
        dir_path (str): temporary dir path
    """
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def decompress_gz_file(file_path, workdir):
    """decompress fna.gz file

    Args:
        file_path (str): fna.gz file path
    """
    with gzip.open(file_path, "rb") as f:
        content = f.read()
    output_name = file_path.split("/")[-1][:-3]
    output_path = os.path.join(workdir, output_name)
    with open(output_path, "wb") as f:
        f.write(content)
