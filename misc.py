
import numpy as np
from io import StringIO
import sys
import subprocess
import pickle

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    return

def save_obj(name, obj):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def load_array(txt):
    s = StringIO(txt)
    arr = np.loadtxt(s)
    return arr


def save_array(arr):
    s = StringIO()
    np.savetxt(s, arr)
    return s.getvalue()


def fix_dir_name(name):

    if not name.endswith("/"):
        name += "/"

    return name


# search

def readlines_reverse(filename):
    with open(filename) as qfile:
        qfile.seek(0, os.SEEK_END)
        position = qfile.tell()
        line = ''
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ''
            else:
                line += next_char
            position -= 1
        yield line[::-1]


def read_line(filename, pattern):

    for i, line in enumerate(readlines_reverse(filename)):
        if line.find(pattern) != -1:
            return line

    return None


def get_index(lines, pattern, offset=None, n_lines=None):

    if offset is None:
        offset = 0

    if n_lines is None:
        n_lines = len(lines)

    for i in range(offset, n_lines):
        line = lines[i]
        if line.find(pattern) != -1:
            return i

    return None


def reverse_enum(L, max_lines=None, lenl=None):

    if lenl is None:
        lenl = len(L)

    if max_lines is None:
        iterator = reversed(range(lenl))
    else:
        iterator = reversed(range(min(lenl, max_lines)))

    for index in iterator:
        yield index, L[index]


def get_rev_index(lines, pattern, max_lines=None, lenl=None, stoppattern=False):

    for i, line in reverse_enum(lines, max_lines=max_lines):

        if line.find(pattern) != -1:
            return i

        if stoppattern and stoppattern in line:
            return None

    return None


def get_indexes(lines, pattern):

    idxs = []

    for i, line in enumerate(lines):
        if pattern in line:
            idxs.append(i)

    return idxs


def get_indexes_with_stop(lines, pattern, stoppattern):

    idxs = []

    for i, line in enumerate(lines):
        # if line.find(pattern) != -1:
        if pattern in line:
            idxs.append(i)
            continue

        # if line.find(stoppattern) != -1:
        if stoppattern in line:
            break

    return idxs


def get_index(lines, pattern):
    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            return i
    return None


def reverse_enum(L):
    for index in reversed(range(len(L))):
        yield index, L[index]


def get_indexes_patterns(lines, patterns):

    n_patterns = len(patterns)
    i_patterns = list(range(n_patterns))

    idxs = [None]*n_patterns

    for i, line in enumerate(lines):

        for ip in i_patterns:

            pattern = patterns[ip]

            if pattern in line:
                idxs[ip] = i
                i_patterns.remove(ip)

    return idxs


def get_rev_indexes(lines, patterns):

    n_patterns = len(patterns)
    i_patterns = list(range(n_patterns))

    idxs = [None]*n_patterns

    for i, line in reverse_enum(lines):

        for ip in i_patterns:

            pattern = patterns[ip]

            if pattern in line:
                idxs[ip] = i
                i_patterns.remove(ip)

    return idxs


def get_rev_index(lines, pattern):

    for i, line in reverse_enum(lines):
        if line.find(pattern) != -1:
            return i

    return None



def shell(cmd, shell=False):
    """

    Run a sh command.

    return the stdout and stderr

    """

    if shell:
        proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        cmd = cmd.split()
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, err = proc.communicate()
    output = output.decode("utf-8")
    err = err.decode("utf-8")
    return output, err

