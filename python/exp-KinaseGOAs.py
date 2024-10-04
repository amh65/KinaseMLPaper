#!/usr/bin/env python3
# Time-stamp: <2024-04-26 11:12:21 smathias>
"""Export TCRD target data to a TSV file.

Usage:
    exp-tsv.py  [--dbhost=<str>] [--dbname=<str>]
    exp-tsv.py -? | --help

Options:
  -h --dbhost DBHOST   : MySQL database host name [default: localhost]
  -n --dbname DBNAME   : MySQL database name [default: tcrd]
  -? --help            : print this message and exit 
"""
__author__ = "Steve Mathias"
__email__ = "smathias@salud.unm.edu"
__org__ = "Translational Informatics Division, UNM School of Medicine"
__copyright__ = "Copyright 2023, Steve Mathias"
__license__ = "Creative Commons Attribution-NonCommercial (CC BY-NC)"
__version__ = "1.0.0"

import os,sys,time
from docopt import docopt
from TCRD.DBAdaptor import DBAdaptor
import slm_util_functions as slmf

PROGRAM = os.path.basename(sys.argv[0])

if __name__ == '__main__':
  print("\n{} (v{}) [{}]:\n".format(PROGRAM, __version__, time.strftime("%c")))
  args = docopt(__doc__, version=__version__)
  
  st = time.time()

  dba = DBAdaptor({'dbname': args['--dbname']})
  dbi = dba.get_dbinfo()
  print("Connected to TCRD database {} (schema ver {}; data ver {})".format(args['--dbname'], dbi['schema_ver'], dbi['data_ver']))

  ofn = "../input_files/TCRDv{}_GOAnnotations.tsv".format(dbi['data_ver'])

  # SELECT distinct go_term FROM goa WHERE evidence IN ('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP') AND go_term REGEXP '^[PF]' AND protein_id IN (SELECT id FROM target WHERE fam = 'Kinase')
  kingoterms = dba.get_all_kinase_goterms()
  header = ['TCRD ID', 'Name', 'Description', 'UniProt', 'HGNC Symbol', 'TDL']
  header.extend(kingoterms)
  tids = dba.find_target_ids_by_fam('Kinase')
  ct = len(tids)
  exp_ct = 0
  print(f"\nExporting GO annotations for {ct} Kinase targets to file {ofn}")
  with open(ofn, 'w') as ofh:
    ofh.write("{}\n".format('\t'.join(header)))
    for tid in tids:
      tpd = dba.get_targetprotein(tid)
      tsvrow = [str(tid), tpd['name'], tpd['description'], tpd['uniprot'], tpd['sym'], tpd['tdl']]
      # SELECT go_term FROM goa WHERE evidence IN ('EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP') AND go_term REGEXP '^[PF]' AND protein_id = %s
      goterms = dba.get_goterms(tid)
      for gt in kingoterms:
        if gt in goterms:
          tsvrow.append('1')
        else:
          tsvrow.append('0')
      assert len(tsvrow) == len(header), f"Bad tsvrow length for target {tid}"
      ofh.write("{}\n".format('\t'.join(tsvrow)))
      exp_ct += 1
      slmf.update_progress(exp_ct/ct)
  print(f"Exported {exp_ct} TSV rows")

  ets = slmf.secs2str(time.time() - st)
  print(f"\n{PROGRAM}: Done. Total time: {ets}\n")


  
  
