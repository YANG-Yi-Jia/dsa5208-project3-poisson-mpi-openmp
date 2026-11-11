"""Run real MPI processes and compare complete fields against Serial."""
import csv
import math
import os
import pathlib
import shlex
import subprocess
import tempfile

ROOT=pathlib.Path(__file__).resolve().parents[1]
LAUNCH=shlex.split(os.environ.get('MPIEXEC','mpirun'))
env=os.environ.copy()
NP_FLAG=os.environ.get('MPI_NP_FLAG','-np')
TIMEOUT=int(os.environ.get('MPI_TEST_TIMEOUT','120'))
records=[]
SINGLE=os.environ.get("MPI_TEST_SINGLE")=="1"
def row_from_stdout(output):
    lines=output.splitlines()
    start=next(i for i,line in enumerate(lines) if line.startswith("version,dim,n,"))
    return list(csv.DictReader(lines[start:start+2]))[0]



def run(command,code=0):
    p=subprocess.run(command,capture_output=True,text=True,env=env,timeout=TIMEOUT)
    assert p.returncode==code,(command,p.returncode,p.stdout,p.stderr)
    return p

for ranks in ((1,) if SINGLE else (1,2,3,4)):
    print(run([*LAUNCH,NP_FLAG,str(ranks),str(ROOT/'test_mpi_stencil')]).stdout.strip())

with tempfile.TemporaryDirectory() as tmp:
    tmp=pathlib.Path(tmp)
    for dim in (2,3):
        for n,ranks in (((1,1),(7,1),(4,1)) if SINGLE else ((1,1),(7,1),(7,2),(7,3),(4,4))):
            for iterations in (1,2,17):
                common=['--dim',str(dim),'--n',str(n),'--iterations',str(iterations)]
                serial_field=tmp/'serial.txt'; mpi_field=tmp/'mpi.txt'
                run([str(ROOT/'poisson_serial'),*common,'--output',str(tmp/'serial.csv'),'--field',str(serial_field)])
                p=run([*LAUNCH,NP_FLAG,str(ranks),str(ROOT/'poisson_mpi'),*common,'--output',str(tmp/'mpi.csv'),'--field',str(mpi_field)])
                a=list(map(float,serial_field.read_text().splitlines())); b=list(map(float,mpi_field.read_text().splitlines()))
                assert len(a)==len(b)==n**dim
                difference=max(abs(x-y) for x,y in zip(a,b))
                assert difference<1e-12,(dim,n,ranks,iterations,difference)
                row=row_from_stdout(p.stdout)
                assert int(row['iterations'])==iterations and float(row['boundary_max'])==0
                rho=math.cos(math.pi/(n+1))
                assert math.isclose(float(row['relative_residual']),abs(rho)**iterations,rel_tol=1e-8,abs_tol=1e-12)
                records.append(dict(dim=dim,n=n,ranks=ranks,iterations=iterations,max_field_difference=difference,test='fixed_iterations'))
        common=['--dim',str(dim),'--n','7','--mode','validation','--tol','1e-10','--check-every','23']
        p=run([*LAUNCH,NP_FLAG,('1' if SINGLE else '3'),str(ROOT/'poisson_mpi'),*common,'--output',str(tmp/'mpi.csv')])
        row=row_from_stdout(p.stdout)
        assert row['converged']=='true' and float(row['relative_residual'])<1e-10
        p=run([*LAUNCH,NP_FLAG,('1' if SINGLE else '3'),str(ROOT/'poisson_mpi'),'--dim',str(dim),'--n','7','--mode','validation','--max-iters','3','--check-every','100','--output',str(tmp/'mpi.csv')],code=2)
        row=row_from_stdout(p.stdout)
        assert row['iterations']=='3' and row['converged']=='false'
    if not SINGLE:
        p=subprocess.run([*LAUNCH,NP_FLAG,('1' if SINGLE else '3'),str(ROOT/'poisson_mpi'),'--n','2'],capture_output=True,text=True,env=env,timeout=TIMEOUT)
        if not SINGLE: assert p.returncode!=0
out=ROOT/'results'/('mpi_single_field_comparison.csv' if SINGLE else 'mpi_field_comparison.csv')
with out.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
print(f'PASS: {len(records)} full-field comparisons, 2 convergence cases, 2 iteration-limit cases; scope={"single process only" if SINGLE else "multiple processes including ranks>N rejection"}')
