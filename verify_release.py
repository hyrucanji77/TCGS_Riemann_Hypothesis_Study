#!/usr/bin/env python3
"""Re-run packaged diagnostics in a fresh directory and audit printed values.

Run: python verify_release.py --output-dir verification/fresh_run
The working subdirectory receives only nine diagnostic scripts and requirements.txt.
No expected JSON, previous numerical result, source PDF, or zero table is copied
there. Comparisons with expected files occur only after computations finish.
Reproducibility is not interval certification or independent mathematical review.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
import difflib
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

SCRIPTS = (
    'diagnostics.py', 'morse_diagnostics.py', 'robin_shift_diagnostics.py',
    'boundary_trace_diagnostics.py', 'scale_endpoint_diagnostics.py',
    'shifted_tower_diagnostics.py', 'completed_tower_diagnostics.py',
    'kinetic_rate_diagnostics.py', 'sesquilinear_diagnostics.py',
)
JOBS = (
    ('diagnostics.py', ['--degree','4','--dps','40'], 'diagnostics_40.json'),
    ('diagnostics.py', ['--degree','4','--dps','60'], 'diagnostics_60.json'),
    ('morse_diagnostics.py', ['--dps','60'], 'morse_diagnostics.json'),
    ('robin_shift_diagnostics.py', ['--dps','60'], 'robin_shift_diagnostics_60.json'),
    ('robin_shift_diagnostics.py', ['--dps','80'], 'robin_shift_diagnostics_80.json'),
    ('boundary_trace_diagnostics.py', [], 'boundary_trace_diagnostics.json'),
    ('scale_endpoint_diagnostics.py', ['--dps','60'], 'scale_endpoint_diagnostics.json'),
    ('shifted_tower_diagnostics.py', ['--degree','3','--dps','40'], 'shifted_tower_diagnostics_40.json'),
    ('shifted_tower_diagnostics.py', ['--degree','3','--dps','60'], 'shifted_tower_diagnostics_60.json'),
    ('completed_tower_diagnostics.py', ['--degree','4','--dps','40'], 'completed_tower_diagnostics_40.json'),
    ('completed_tower_diagnostics.py', ['--degree','4','--dps','60'], 'completed_tower_diagnostics_60.json'),
    ('kinetic_rate_diagnostics.py', [], 'kinetic_rate_diagnostics.json'),
    ('sesquilinear_diagnostics.py', [], 'sesquilinear_diagnostics.json'),
)

def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()

def strip_environment(value):
    if isinstance(value, dict):
        return {k: strip_environment(v) for k,v in value.items()
                if k not in {'python', 'python_version'}}
    if isinstance(value, list):
        return [strip_environment(v) for v in value]
    return value

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('verification/fresh_run'))
    parser.add_argument('--jobs', type=int, default=2, help='Parallel processes, 1..4.')
    parser.add_argument('--timeout', type=int, default=600, help='Seconds per diagnostic.')
    args=parser.parse_args()
    if not 1 <= args.jobs <= 4 or args.timeout < 1:
        parser.error('Require jobs 1..4 and a positive timeout.')
    import mpmath
    if mpmath.__version__ != '1.3.0':
        parser.error('Install requirements.txt: mpmath==1.3.0.')
    root=Path(__file__).resolve().parent
    out=args.output_dir.resolve()
    if out.exists():
        parser.error('Output directory exists; choose a fresh path.')
    for name in (*SCRIPTS, 'requirements.txt', 'main.tex', *(j[2] for j in JOBS)):
        if not (root/name).is_file():
            parser.error('Missing release input: '+name)
    work=out/'work'
    work.mkdir(parents=True)
    (out/'logs').mkdir()
    for name in (*SCRIPTS,'requirements.txt'):
        shutil.copy2(root/name, work/name)
    initial_files=sorted(p.name for p in work.iterdir())
    assert not list(work.glob('*.json')) and not list(work.glob('*.pdf'))
    env=dict(os.environ)
    env.update(PYTHONPATH='', PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1')
    started=time.time()
    def run(job):
        script,flags,filename=job
        cmd=[sys.executable,'-B',script,*flags,'--output',filename]
        begin=time.time()
        try:
            process=subprocess.run(cmd,cwd=work,env=env,capture_output=True,
                                   timeout=args.timeout,check=False)
            rc,stdout,stderr=process.returncode,process.stdout,process.stderr
        except subprocess.TimeoutExpired as error:
            rc,stdout,stderr=-1,error.stdout or b'',(error.stderr or b'')+b'\nTimeout.\n'
        (out/'logs'/f'{filename}.stdout.txt').write_bytes(stdout)
        (out/'logs'/f'{filename}.stderr.txt').write_bytes(stderr)
        print(f'{filename}: exit {rc}, {time.time()-begin:.2f} s',flush=True)
        return {'script':script,'arguments':flags+['--output',filename],
                'output':filename,'returncode':rc,'elapsed_seconds':round(time.time()-begin,3)}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs=list(pool.map(run,JOBS))
    # Expected files are compared only after all computations complete.
    for row in jobs:
        expected,actual=root/row['output'],work/row['output']
        row['expected_sha256']=digest(expected)
        if not actual.exists():
            row.update(byte_identical=False,json_equal_except_python_metadata=False)
            continue
        row['actual_sha256']=digest(actual)
        row['byte_identical']=actual.read_bytes()==expected.read_bytes()
        row['json_equal_except_python_metadata']=(
            strip_environment(json.loads(actual.read_text())) ==
            strip_environment(json.loads(expected.read_text())))
        if not row['byte_identical']:
            diff=''.join(difflib.unified_diff(expected.read_text().splitlines(True),
                      actual.read_text().splitlines(True),fromfile='expected/'+row['output'],
                      tofile='fresh/'+row['output']))
            (out/'logs'/f'{row["output"]}.diff').write_text(diff,encoding='utf-8')
    checks=[]
    def printed(name, value, token):
        cleaned=token.replace(' ','').replace('$','')
        match=re.fullmatch(r'([-+]?\d+(?:\.\d+)?)(?:\\times10\^\{([-+]?\d+)\})?',cleaned)
        if not match:
            raise ValueError('Unrecognized printed numeric token: '+token)
        mantissa,exponent=match.group(1),int(match.group(2) or '0')
        places=len(mantissa.split('.')[1]) if '.' in mantissa else 0
        with localcontext() as ctx:
            ctx.prec=100
            expected=Decimal(mantissa)*(Decimal(10)**exponent)
            unit=Decimal(10)**(exponent-places)
            delta=abs(Decimal(value)-expected)
            ok=delta<=unit/2
        checks.append({'location':name,'printed_latex':token,'computed':value,
                       'within_half_last_printed_unit':ok})
    if all(r['returncode']==0 for r in jobs):
        tex=(root/'main.tex').read_text(encoding='utf-8')
        def table_rows(label):
            start=tex.index('\\label{'+label+'}')
            region=tex[start:tex.index('\\end{tabular}',start)]
            return [re.findall(r'\$([^$]+)\$',line) for line in region.splitlines()
                    if re.match(r'^\d+\s*&',line)]
        moment_rows=table_rows('tab:moments')
        tower_rows=table_rows('tab:determinants')
        if len(moment_rows)!=10 or len(tower_rows)!=5:
            raise ValueError('Frozen table shape differs from the article.')
        for dps in (40,60):
            data=json.loads((work/f'completed_tower_diagnostics_{dps}.json').read_text())
            for n,row in enumerate(moment_rows):
                printed(f'Table 1 / m_{n} / {dps} dps',data['moments'][n],row[0])
            for n,row in enumerate(tower_rows):
                printed(f'Table 2 / H_{n} / {dps} dps',data['scaled_unshifted_determinants'][n],row[0])
                printed(f'Table 2 / shifted H_{n} / {dps} dps',data['scaled_shifted_determinants'][n],row[1])
        stop=tex.index('\\label{eq:robinvalues}')
        region=tex[tex.rfind('\\begin{equation}',0,stop):stop]
        tokens=re.findall(r'\\approx\s*([-+]?\d+\.\d+)',region)
        if len(tokens)!=6: raise ValueError('Expected six values in Equation (13.20).')
        begin=tex.index('\\label{eq:robinsecondnumerical}')
        region=tex[begin:tex.index('They approach',begin)]
        residuals=re.findall(r'\$([-+]?\d+\.\d+)\$',region)
        if len(residuals)!=4: raise ValueError('Expected four values after Equation (13.21).')
        for dps in (60,80):
            data=json.loads((work/f'robin_shift_diagnostics_{dps}.json').read_text())
            for key,token in zip(('h','W_at_zero_order','W_x_at_zero_order','Q_at_zero_order','C_R','B'),tokens):
                printed(f'Equation (13.20) / {key} / {dps} dps',data['matched'][key],token)
            for row,token in zip(data['asymptotic_checks'],residuals):
                printed(f'Equation (13.21) / mu={row["mu"]} / {dps} dps',
                        row['mu_squared_log_ratio_over_C'],token)
        counts={
            'boundary_trace_diagnostics.py':json.loads((work/'boundary_trace_diagnostics.json').read_text())['number_of_checks'],
            'scale_endpoint_diagnostics.py':json.loads((work/'scale_endpoint_diagnostics.json').read_text())['total_additional_exact_checks'],
            'robin_shift_diagnostics.py':json.loads((work/'robin_shift_diagnostics_60.json').read_text())['exact_algebra']['count'],
            'kinetic_rate_diagnostics.py':json.loads((work/'kinetic_rate_diagnostics.json').read_text())['exact_checks_passed'],
            'sesquilinear_diagnostics.py':json.loads((work/'sesquilinear_diagnostics.json').read_text())['check_count'],
        }
    else:
        counts={}
    success=(all(r['returncode']==0 and r['json_equal_except_python_metadata'] for r in jobs)
             and len(checks)==60 and all(c['within_half_last_printed_unit'] for c in checks))
    result={'status':'PASS' if success else 'FAIL','python':platform.python_version(),
            'mpmath':mpmath.__version__,'platform':platform.platform(),
            'elapsed_seconds':round(time.time()-started,3),'initial_work_files':initial_files,
            'initial_work_had_no_expected_outputs':True,'jobs':jobs,
            'all_outputs_byte_identical':all(r['byte_identical'] for r in jobs),
            'printed_comparisons':checks,'printed_comparison_count':len(checks),
            'distinct_exact_check_counts':counts,'distinct_exact_checks_total':sum(counts.values()),
            'source_script_hashes':{name:digest(root/name) for name in SCRIPTS},
            'manuscript_sha256':digest(root/'main.tex'),
            'limitations':['Same installed Python/mpmath runtime; fresh working directory, not a separate machine.',
              'No interval certificates for floating-point special-function values.',
              'Reproducibility is not independent peer review or an RH proof.']}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'jobs':len(jobs),
                     'byte_identical':result['all_outputs_byte_identical'],
                     'printed_comparisons':len(checks),'distinct_exact_checks':sum(counts.values())},indent=2))
    if not success: sys.exit(1)

if __name__=='__main__':
    main()
