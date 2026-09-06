"""Run unchanged precomputation sources in a new external directory."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,zipfile
ROOT=Path(__file__).resolve().parent.parent
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--work-dir',type=Path,required=True)
    ap.add_argument('--full',action='store_true')
    args=ap.parse_args();work=args.work_dir.resolve()
    if work==ROOT or ROOT in work.parents:
        raise ValueError('Generation must be outside the frozen package')
    freeze=json.loads((ROOT/'stochastic_diffusion/results/FREEZE.json').read_text(encoding='utf-8'))
    with zipfile.ZipFile(ROOT/'provenance/diffusion_precomputation.zip') as z:
        if set(z.namelist())!={'PROTOCOL.md','run_benchmark.py'}:
            raise ValueError('Unexpected precomputation archive contents')
        for name,key in [('PROTOCOL.md','protocol_sha256'),('run_benchmark.py','generator_sha256')]:
            if hashlib.sha256(z.read(name)).hexdigest()!=freeze[key]:
                raise ValueError('Original source hash mismatch: '+name)
        work.mkdir(parents=True,exist_ok=False)
        for name in ('PROTOCOL.md','run_benchmark.py'):(work/name).write_bytes(z.read(name))
    env=os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([sys.executable,'-B','run_benchmark.py'],cwd=work,env=env,check=True)
    if args.full:subprocess.run([sys.executable,'-B','run_benchmark.py','--full'],cwd=work,env=env,check=True)
    print('Retained generation directory:',work)
if __name__=='__main__':main()
