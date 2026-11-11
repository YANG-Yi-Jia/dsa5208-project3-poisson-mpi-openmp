"""Independent sine-eigenmode oracle, refinement, CLI and stopping checks."""
import csv
import math
import pathlib
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
rows = []

def run(dim, n, *extra, code=0):
    with tempfile.TemporaryDirectory() as tmp:
        p = subprocess.run([str(ROOT/'poisson_serial'), '--dim', str(dim), '--n', str(n),
                            '--output', str(pathlib.Path(tmp)/'result.csv'), *extra], capture_output=True, text=True)
        assert p.returncode == code, (p.args, p.returncode, p.stderr)
        if code == 1:
            return
        row = list(csv.DictReader(p.stdout.splitlines()))[0]
        rows.append(row)
        assert float(row['boundary_max']) == 0
        return row

for dim in (2, 3):
    for n in (1, 4, 9):
        k=7
        row=run(dim,n,'--iterations',str(k))
        theta=math.pi/(n+1)
        rho=math.cos(theta)
        discrete_factor=(theta/(2*math.sin(theta/2)))**2
        amplitude=discrete_factor*(1-rho**k)
        max_sine=max(math.sin(theta*i) for i in range(1,n+1))
        expected_inf=abs(amplitude-1)*max_sine**dim
        rms_sine=math.sqrt((n+1)/(2*n))
        assert math.isclose(float(row['linf_error']),expected_inf,rel_tol=1e-11,abs_tol=1e-13)
        assert math.isclose(float(row['l2_rms_error']),abs(amplitude-1)*rms_sine**dim,rel_tol=1e-11,abs_tol=1e-13)
        assert math.isclose(float(row['relative_residual']),abs(rho)**k,rel_tol=1e-9,abs_tol=1e-13)
        assert int(row['iterations'])==k
    errors=[]
    for n in (7,15,31):
        row=run(dim,n,'--mode','validation','--tol','1e-11','--check-every','37')
        assert float(row['relative_residual']) < 1e-11
        assert row['converged']=='true'
        errors.append(float(row['linf_error']))
    orders=[math.log2(a/b) for a,b in zip(errors,errors[1:])]
    assert all(1.95<p<2.05 for p in orders), orders
    print(f'{dim}D max-error refinement orders: {orders}')
    row=run(dim,9,'--mode','validation','--max-iters','3','--check-every','100',code=2)
    assert row['iterations']=='3' and row['converged']=='false'
for args in [('--dim','4'),('--n','0'),('--n','2abc'),('--tol','nan'),('--check-every','0'),('--mode','other'),('--unknown','2')]:
    run(2,4,*args,code=1)
with tempfile.TemporaryDirectory() as tmp:
    path=pathlib.Path(tmp)/'append.csv'
    for _ in range(2):
        subprocess.run([str(ROOT/'poisson_serial'),'--n','3','--iterations','1','--output',str(path)],check=True,capture_output=True)
    assert len(list(csv.DictReader(path.open())))==2
out=ROOT/'results'/'local_validation.csv'
with out.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
print(f'PASS: {len(rows)} numerical/stopping cases, 7 invalid arguments, CSV append; {out.name}')
