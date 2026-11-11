#!/usr/bin/env python3
"""Compare complete hybrid fields against serial for several rank/thread settings."""
import os, pathlib, subprocess, tempfile, shlex, csv
ROOT=pathlib.Path(__file__).resolve().parents[1]
LAUNCH=shlex.split(os.environ.get('MPIEXEC','mpirun'))
FLAG=os.environ.get('MPI_NP_FLAG','-np')
with tempfile.TemporaryDirectory() as d:
    d=pathlib.Path(d)
    count=0
    for dim in (2,3):
        for n in (4,7):
            for iterations in (1,2,9):
                args=['--dim',str(dim),'--n',str(n),'--iterations',str(iterations)]
                subprocess.run([str(ROOT/'poisson_serial'),*args,'--field',str(d/'serial'),'--output',str(d/'s.csv')],check=True,capture_output=True)
                a=[float(x) for x in (d/'serial').read_text().split()]
                for ranks,threads in ((1,1),(1,2),(2,2),(3,1)):
                    env=os.environ.copy();env['OMP_NUM_THREADS']=str(threads);env['OMP_DYNAMIC']='FALSE'
                    subprocess.run([*LAUNCH,FLAG,str(ranks),str(ROOT/'poisson_hybrid'),*args,'--field',str(d/'hybrid'),'--output',str(d/'h.csv')],check=True,capture_output=True,env=env,timeout=120)
                    b=[float(x) for x in (d/'hybrid').read_text().split()]
                    assert len(a)==len(b)==n**dim
                    diff=max(abs(x-y) for x,y in zip(a,b))
                    assert diff<1e-12,(dim,n,iterations,ranks,threads,diff)
                    count+=1
    print(f'PASS: {count} hybrid full-field comparisons across ranks/threads')
