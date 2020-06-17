
import os

def get_pinfo():
    """
    get process id of parent and current process
    """
    ppid = os.getppid()
    pid = os.getpid()
    return ppid, ppid

