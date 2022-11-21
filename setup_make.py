import sys
import subprocess
subprocess.check_call([sys.executable,'-m','pip','install','numpy'])
subprocess.check_call([sys.executable,'-m','pip','install','astropy'])
subprocess.check_call([sys.executable,'-m','pip','install','ccdproc'])
subprocess.check_call([sys.executable,'-m','pip','install','scipy'])
subprocess.check_call([sys.executable,'-m','pip','install','matplotlib'])
