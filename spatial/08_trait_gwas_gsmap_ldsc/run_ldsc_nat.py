#!/usr/bin/env python
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re
import glob
import logging
import subprocess


# File:       LDSC.py
# Description: This script processes multiple BED files for cell-type-specific LDSC analysis. It generates annotation files, LD scores, and performs enrichment analysis.
# Author:      jmLiu
# Email:       liujiamao.xmu.edu
# Date:        07/08/2025
# Version:     1.2
# Python:      Python 2.7
# Software:    LDSC: https://github.com/bulik/ldsc/
# Ref:         https://github.com/yal054/snATACutils/tree/master/05.LDSC_analysis; https://github.com/ZunpengLiu/Multi-region_AD/tree/main/01_snATAC/08_GWAS_enrichment


MY_USAGE = '''
    python LDSC.py -i /path/to/input -o /path/to/output -q /path/to/qsub -c /path/to/configure [-r]
    Example:
    python LDSC.py -i ./bed_files -o ./results -q ./qsub_scripts -c ./config.txt -r
'''

parser = OptionParser(MY_USAGE, version='Version 1.0')
parser.add_option('-i', '--in_dir', dest='in_dir', type='string', help='Input directory containing BED files')
parser.add_option('-o', '--out_dir', dest='out_dir', type='string', help='Output directory for results')
parser.add_option('-q', '--qsub_dir', dest='qsub_dir', type='string', help='Directory for qsub scripts')
parser.add_option('-c', '--configure_file', dest='configure_file', type='string', help='Configuration file with software paths')
parser.add_option('-r', '--run', dest='run_scripts', action='store_true', default=False, help='Run generated scripts immediately')
(options, args) = parser.parse_args()

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Initialize directories
ensure_dir(options.out_dir)
ensure_dir(options.qsub_dir)

# Load software paths
sft = {}
with open(options.configure_file, 'r') as config:
    for line in config:
        if len(re.findall('=', line)) == 1:
            key = re.findall(r'(\w+?)=.*', line)[0]
            value = re.findall(r'\w+=(.*)', line)[0]
            sft[key] = value

# Get all BED files
bed_files = glob.glob(os.path.join(options.in_dir, '*.bed'))
if not bed_files:
    sys.exit("Error: No BED files found in input directory!")

# Combine all BED files into hsc.bed as background
hsc_bed = os.path.join(options.in_dir, 'hsc.bed')
with open(hsc_bed, 'w') as outfile:
    for bed_file in bed_files:
        with open(bed_file, 'r') as infile:
            for line in infile:
                outfile.write(line)
bed_files.append(hsc_bed)  # Add hsc.bed to the list for further processing

# Create subdirectories
ensure_dir(os.path.join(options.out_dir, '02_ldsc'))
ensure_dir(os.path.join(options.out_dir, '04_Enrichment_result'))
ensure_dir(os.path.join(options.qsub_dir, '02_ldsc_annot'))
ensure_dir(os.path.join(options.qsub_dir, '03_ldsc_l2'))
ensure_dir(os.path.join(options.qsub_dir, '04_Enrichment'))

# Process each BED file
for bed_file in bed_files:
    cell_type = os.path.basename(bed_file).replace('.bed', '')
    
    # 1. Create annotation files
    annot_script = os.path.join(options.qsub_dir, '02_ldsc_annot', 'ldsc_annot.%s.sh' % cell_type)
    with open(annot_script, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write('echo "Creating annotations for %s"\n' % cell_type)
        # Use absolute path for activation
        f.write('source /cluster2/huanglab/jiamao/conda/bin/activate ldsc\n')
        
        for chr in range(1, 23):
            annot_file = os.path.join(options.out_dir, '02_ldsc', '%s.%s.annot.gz' % (cell_type, chr))
            
            # === Check 1: Skip if annot file exists ===
            if os.path.exists(annot_file) and os.path.getsize(annot_file) > 0:
                continue
            # ==========================================

            cmd = ' '.join([
                sft['python'], sft['make_annot'],
                '--bed-file', bed_file,
                '--bimfile', sft['bimfile'] + str(chr) + '.bim',
                '--annot-file', annot_file
            ]) + ' &\n'
            f.write(cmd)
        f.write('wait\necho "Done creating annotations for %s"\n' % cell_type)
    os.chmod(annot_script, 0755)

    # 2. Compute LD scores
    ldscore_script = os.path.join(options.qsub_dir, '03_ldsc_l2', 'ldsc_l2.%s.sh' % cell_type)
    with open(ldscore_script, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write('echo "Computing LD scores for %s"\n' % cell_type)
        # Use absolute path for activation
        f.write('source /cluster2/huanglab/jiamao/conda/bin/activate ldsc\n')
        
        for chr in range(1, 23):
            annot_file = os.path.join(options.out_dir, '02_ldsc', '%s.%s.annot.gz' % (cell_type, chr))
            out_prefix = os.path.join(options.out_dir, '02_ldsc', '%s.%s' % (cell_type, chr))
            
            # === Check 2: Skip if ldscore file exists ===
            # LDSC outputs [prefix].l2.ldscore.gz
            expected_ld = out_prefix + '.l2.ldscore.gz'
            if os.path.exists(expected_ld) and os.path.getsize(expected_ld) > 0:
                continue
            # ============================================

            cmd = ' '.join([
                sft['python'], sft['ldsc_py'],
                '--l2',
                '--bfile', sft['bfile'] + str(chr),
                '--ld-wind-cm', '1',
                '--annot', annot_file,
                '--out', out_prefix,
                '--thin-annot',
                '--print-snps', sft['snps'] + 'txt'
            ]) + '\n'
            f.write(cmd)
        f.write('echo "Done computing LD scores for %s"\n' % cell_type)
    os.chmod(ldscore_script, 0755)

# 3. Generate LDCTS file
ldcts_file = os.path.join(options.out_dir, 'celltype.ldcts')
hsc_prefix = os.path.join(options.out_dir, '02_ldsc', 'hsc.')
with open(ldcts_file, 'w') as f:
    for bed_file in bed_files:
        cell_type = os.path.basename(bed_file).replace('.bed', '')
        if cell_type != 'hsc':  # Skip hsc.bed itself
            l2_prefix = os.path.join(options.out_dir, '02_ldsc', cell_type + '.')
            f.write("%s\t%s,%s\n" % (cell_type, l2_prefix, hsc_prefix))


# 4. Prepare enrichment analysis
def getbed(bed_file):
    beds = {}
    with open(bed_file, 'r') as b:
        for line in b:
            line = line.strip()
            if line and not line.startswith('#'):
                if len(re.findall('=', line)) == 1:
                    key = re.findall(r'(\w+?)=.*', line)[0]
                    value = re.findall(r'\w+=(.*)', line)[0]
                    beds[key] = value
    return beds

# Load GWAS sumstats from sumstat.txt
beds = getbed("sumstat.txt")
if not beds:
    sys.exit("Error: No GWAS sumstats found in sumstat.txt!")

# Create enrichment scripts for each GWAS trait
for trait, sumstats_path in beds.items():
    
    # === Check 3: Skip if Enrichment result exists ===
    # LDSC outputs [out].cell_type_results.txt
    result_file = os.path.join(options.out_dir, '04_Enrichment_result', trait)
    check_result = result_file + '.cell_type_results.txt'
    if os.path.exists(check_result) and os.path.getsize(check_result) > 0:
        continue
    # =================================================

    enrich_script = os.path.join(options.qsub_dir, '04_Enrichment', 'enrichment_%s.sh' % trait)
    with open(enrich_script, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write('echo "Running enrichment analysis for %s"\n' % trait)
        # Use absolute path for activation
        f.write('source /cluster2/huanglab/jiamao/conda/bin/activate ldsc\n')
        
        cmd = ' '.join([
            sft['python'], sft['ldsc_py'],
            '--h2-cts', sumstats_path,
            '--ref-ld-chr', sft['baselineLD'],
            '--out', result_file,
            '--ref-ld-chr-cts', ldcts_file,
            '--w-ld-chr', sft['weights']
        ]) + '\n'
        f.write(cmd)
        f.write('echo "Enrichment analysis for %s complete"\n' % trait)
    os.chmod(enrich_script, 0755)

# logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

logger.info("Processing complete for the following cell types:")
for bed_file in bed_files:
    logger.info(" - %s", os.path.basename(bed_file).replace('.bed', ''))
logger.info("GWAS traits processed:")
for trait in beds.keys():
    logger.info(" - %s", trait)
logger.info("LDCTS file created at: %s", ldcts_file)
logger.info("Scripts generated in: %s", options.qsub_dir)


# Automatic execution
def run_script(script_path):
    """run single script"""
    try:
        subprocess.check_call(['sh', script_path])
        logging.info("Successfully ran script: %s" % script_path)
        return True
    except subprocess.CalledProcessError as e:
        logging.error("Failed to run script %s: %s" % (script_path, str(e)))
        return False

def run_all_scripts(script_dir):
    """run all script in dir"""
    scripts = glob.glob(os.path.join(script_dir, '*.sh'))
    if not scripts:
        logging.warning("No scripts found in %s" % script_dir)
        return
    
    logging.info("Running all scripts in %s" % script_dir)
    for script in sorted(scripts):
        run_script(script)

if options.run_scripts:
    logging.info("Auto-running generated scripts as requested...")
    
    # order 02_ldsc_annot→ 03_ldsc_l2→ 04_Enrichment
    run_all_scripts(os.path.join(options.qsub_dir, '02_ldsc_annot'))
    run_all_scripts(os.path.join(options.qsub_dir, '03_ldsc_l2'))
    run_all_scripts(os.path.join(options.qsub_dir, '04_Enrichment'))
    
    logging.info("All scripts have been submitted for execution.")
else:
    logging.info("Scripts generated but not executed (use -r to run them automatically)")