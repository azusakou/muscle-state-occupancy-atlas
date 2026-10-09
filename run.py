#!/usr/bin/env python3
"""Reproduce from shipped processed data; dependencies must already be installed."""
from pathlib import Path
import os,sys,subprocess,importlib,hashlib,json,time
ROOT=Path(__file__).resolve().parent;os.chdir(ROOT)
started=time.time()
# Fail before analysis if an input or code file differs from the delivered release.
if not (ROOT/'SHA256SUMS').is_file():raise SystemExit('Missing required SHA256SUMS release manifest')
for line in (ROOT/'SHA256SUMS').read_text().splitlines():
 expected,path=line.split('  ',1)
 target=(ROOT/path).resolve()
 if not target.is_relative_to(ROOT) or not target.is_file():raise SystemExit('Invalid or missing release file: '+path)
 if hashlib.sha256(target.read_bytes()).hexdigest()!=expected:raise SystemExit('Checksum mismatch: '+path)
for p in ['numpy','pandas']:importlib.import_module(p)
rs=os.environ.get('RSCRIPT','Rscript')
subprocess.run([rs,'-e','for(p in c("DESeq2","data.table","Matrix","glmnet","nnet","xgboost","ggplot2","patchwork","scales","ragg")) if(!requireNamespace(p,quietly=TRUE)) stop(paste("Missing R package:",p))'],check=True)
(ROOT/'results').mkdir(exist_ok=True);(ROOT/'figures').mkdir(exist_ok=True)
steps=[('test_encoding.R',[]),('refit_lineage.R',[]),('derive_dynamics.R',[]),('derive_validation.R',[]),('derive_denervation.R',[]),
 ('occupancy_sensitivity.py',[]),('fusion_design_sensitivity.R',[]),
 ('benchmark_train_only.R',[]),('controls_train_only.R',[]),
 ('controls_train_only.R',['--exclude-ambiguous']),('uncertainty_train_only.R',[]),
 ('revise_figures.R',[]),('validate_analysis.py',[])]
for name,args in steps:
 command=[sys.executable if name.endswith('.py') else rs,'code/'+name,*args]
 print('RUN: '+' '.join(command),flush=True)
 with (ROOT/'results'/('log_'+name+('_11' if args else '')+'.txt')).open('w') as log:
  subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
print('PASS: analyses, five figures and independent numeric checks; elapsed seconds:',round(time.time()-started,1))
