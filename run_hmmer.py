#!/usr/bin/env pytho
import subprocess
import sys
import os


class RunHmmer(object):
    #based on the biopython application abstract base class

    def __init__(self, hmmer_cmd, hmmfile, seqfile, *args):
        
        """instatiate hmmer command line wrapper object"""
        
        self.hmmer_cmd = hmmer_cmd
        self._options = list(args)
        self._hmmfile = hmmfile
        self._seqfile = seqfile
        #HMMER_DB="/home/drewx/Dropbox/In silico identification of caspase-like homologs/HMMs"    
        #assert os.environ.get('HMMER_DB',False), 'HMMER_DB environmental variable is not set' 

    def __str__(self):
        
        cmd = "{} ".format(self.hmmer_cmd)
        for opt in self.options:
            cmd = " ".join([cmd, opt])
            
        cmd  =  " ".join(cmd, self.hmmer_cmd, self.seqfile)
            
        return cmd.strip()

    
    @property
    def options(self):
          return self._options

    @options.setter
    def options(self, option):
        
        if self._options:
            self._options.append(option)
        else:
            self._options = [option]
            
    def cmd_clear(self):    
        self._options = []


    @property
    def hmmfile(self)

        
    def __call__(self, stdin=None, stdout=True, stderr=True, cwd=None, env=None):
        #this code is modified from biopython Applicaitons module:__init__.py
        
        if not stdout:
            stdout_arg = open(os.devnull, "w")
        elif isinstance(stdout, str):
            stdout_arg = open(stdout, "w")
        else:
            stdout_arg = subprocess.PIPE

        if not stderr:
            stderr_arg = open(os.devnull, "w")
        elif isinstance(stderr, str):
            if stdout == stderr:
                stderr_arg = stdout_arg  # Write both to the same file
            else:
                stderr_arg = open(stderr, "w")
        else:
            stderr_arg = subprocess.PIPE

        if sys.platform != "win32":
            use_shell = True
        else:
            win_ver = platform.win32_ver()[0]
            if win_ver in ["7", "8", "post2012Server", "10"]:
                use_shell = True
            else:
                use_shell = False
        child_process = subprocess.Popen(
            str(self),
            stdin=subprocess.PIPE,
            stdout=stdout_arg,
            stderr=stderr_arg,
            universal_newlines=True,
            cwd=cwd,
            env=env,
            shell=use_shell)
        
        # Use .communicate as can get deadlocks with .wait(), see Bug 2804
        stdout_str, stderr_str = child_process.communicate(stdin)
        if not stdout:
            assert not stdout_str, stdout_str
        if not stderr:
            assert not stderr_str, stderr_str
        return_code = child_process.returncode

        # Particularly important to close handles on Jython and PyPy
        # (where garbage collection is less predictable) and on Windows
        # (where cannot delete files with an open handle):
        if not stdout or isinstance(stdout, str):
            # We opened /dev/null or a file
            stdout_arg.close()
        if not stderr or (isinstance(stderr, str) and stdout != stderr):
            # We opened /dev/null or a file
            stderr_arg.close()

        if return_code:
            raise Exception(return_code, str(self), stdout_str, stderr_str)
        return stdout_str, stderr_str
 

if __name__ == '__main__':
    hmmsearch = RunHmmer('hmmersearch')
    print(hmmsearch)
    #hmmalign [-options] <hmmfile> <seqfile>
    #parser = argparse.ArgumentParser(description="")
    #parser.add_argument('seqfile', action='store', type=str)
    #args = parser.parse_args()
    
