#!/usr/bin/env python3
"""Run fixed workloads in fresh R processes; retain elapsed time, RSS and outputs."""
import argparse, csv, json, platform, re, subprocess
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--source', required=True); p.add_argument('--output',required=True); p.add_argument('--repeats',type=int,default=3); a=p.parse_args()
out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=True)
script=Path(__file__).with_name('measure.R').resolve()
rows=[]
for case in ['small','sparse','branching']:
 for operation in ['fit','transfer','landmarks','embedding']:
  for repeat in range(1,a.repeats+1):
   prefix=out/f'{case}-{operation}-{repeat}'
   time_args=['/usr/bin/time','-l'] if platform.system()=='Darwin' else ['/usr/bin/time','-v']
   result=subprocess.run(time_args+['Rscript',str(script),a.source,case,operation,str(prefix)],capture_output=True,text=True)
   prefix.with_suffix('.log').write_text(result.stdout+'\n'+result.stderr)
   if result.returncode: raise SystemExit(f'Failed: {prefix}; inspect its log')
   row=next(csv.DictReader(prefix.with_suffix('.csv').open()))
   pattern=r'(\d+)\s+maximum resident set size' if platform.system()=='Darwin' else r'Maximum resident set size \(kbytes\):\s*(\d+)'
   match=re.search(pattern,result.stderr)
   row.update({'repeat':repeat,'peak.process.rss.bytes':int(match[1])*(1 if platform.system()=='Darwin' else 1024)})
   rows.append(row)
   print(case,operation,repeat,row['elapsed.seconds'],flush=True)
with (out/'results.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(out/'environment.json').write_text(json.dumps({'platform':platform.platform(),'machine':platform.machine(),'source':str(Path(a.source).resolve()),'repeats':a.repeats,'R':subprocess.check_output(['Rscript','--version'],stderr=subprocess.STDOUT,text=True)},indent=2)+'\n')
