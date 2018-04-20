#! /usr/bin/env python3

import logging
import json
import os
import subprocess
import sys
import time

ENA_URL = ('https://www.ebi.ac.uk/ena/data/warehouse/search?result=read_run&'
           'display=report')
FIELDS = [
    'study_accession', 'secondary_study_accession', 'sample_accession',
    'secondary_sample_accession', 'experiment_accession', 'run_accession',
    'submission_accession', 'tax_id', 'scientific_name',
    'instrument_platform', 'instrument_model', 'library_name',
    'library_layout', 'nominal_length', 'library_strategy',
    'library_source', 'library_selection', 'read_count',
    'base_count', 'center_name', 'first_public', 'last_updated',
    'experiment_title', 'study_title', 'study_alias', 'experiment_alias',
    'run_alias', 'fastq_bytes', 'fastq_md5', 'fastq_ftp', 'fastq_aspera',
    'fastq_galaxy', 'submitted_bytes', 'submitted_md5', 'submitted_ftp',
    'submitted_aspera', 'submitted_galaxy', 'submitted_format',
    'sra_bytes', 'sra_md5', 'sra_ftp', 'sra_aspera', 'sra_galaxy',
    'cram_index_ftp', 'cram_index_aspera', 'cram_index_galaxy',
    'sample_alias', 'broker_name', 'sample_title', 'nominal_sdev',
    'first_created'
]


def output_handler(output, redirect='>'):
    if output:
        return [open(output, 'w'), '{0} {1}'.format(redirect, output)]
    else:
        return [subprocess.PIPE, '']


def onfinish_handler(cmd, out, err, returncode):
    out = '\n{0}'.format(out) if out else ''
    err = '\n{0}'.format(err) if err else ''
    if returncode != 0:
        logging.error('COMMAND: {0}'.format(cmd))
        logging.error('STDOUT: {0}'.format(out))
        logging.error('STDERR: {0}'.format(err))
        logging.error('END\n'.format(err))
        raise RuntimeError(err)
    else:
        logging.info('COMMAND: {0}'.format(cmd))
        logging.info('STDOUT: {0}'.format(out))
        logging.info('STDERR: {0}'.format(err))
        logging.info('END\n'.format(err))
        return [out, err]


def byte_to_string(b):
    if b:
        return b.decode("utf-8")
    else:
        return ''


def run_command(cmd, cwd=os.getcwd(), stdout=False, stderr=False):
    """Execute a single command and return STDOUT and STDERR."""
    stdout, stdout_str = output_handler(stdout)
    stderr, stderr_str = output_handler(stderr, redirect='2>')
    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd)

    out, err = p.communicate()
    return onfinish_handler(
        '{0} {1} {2}'.format(' '.join(cmd), stdout_str, stderr_str),
        byte_to_string(out), byte_to_string(err), p.returncode
    )


def log_stdout(message, quiet=False):
    logging.info(message)


def md5sum(file):
    """Return the MD5SUM of an input file."""
    if os.path.exists(file):
        stdout, stderr = run_command(['md5sum', file])
        if stdout:
            md5sum, filename = stdout.split()
            return md5sum
        else:
            return None
    else:
        return None


def download_fastq(fasp, ftp, outdir, md5, max_retry=10, use_ftp=False):
    """Download FASTQ from ENA using Apera Connect."""
    success = False
    retries = 0
    fastq = '{0}/{1}'.format(
        outdir, format(os.path.basename(fasp))
    )

    if not os.path.exists(fastq):
        if not os.path.isdir(outdir):
            run_command(['mkdir', '-p', outdir])

        while not success:
            if use_ftp:
                run_command(['wget', '-O', fastq, ftp])
            else:
                run_command([os.environ['ASCP'], '-QT', '-l', '300m',
                            '-P33001', '-i', os.environ['ASCP_KEY'],
                            'era-fasp@{0}'.format(fasp), outdir])

            if md5sum(fastq) != md5:
                retries += 1
                if os.path.exists(fastq):
                    os.remove(fastq)
                if retries > max_retry:
                    if not use_ftp:
                        use_ftp = True
                        retries = 0
                    else:
                        break
                time.sleep(10)
            else:
                success = True
    else:
        success = True

    return [success, fastq]


def merge_runs(runs, output):
    """Merge runs from an experiment."""
    if len(runs) > 1:
        cat_cmd = ['cat']
        rm_cmd = ['rm']
        for run in runs:
            cat_cmd.append(run)
            rm_cmd.append(run)
        run_command(cat_cmd, stdout=output)
        run_command(rm_cmd)
    else:
        run_command(['mv', runs[0], output])


def get_run_info(experiment):
    """Retreive a list of unprocessed samples avalible from ENA."""
    import requests
    url = '{0}&query="{1}"&fields={2}'.format(ENA_URL, query, ",".join(FIELDS))
    r = requests.get(url)
    if r.status_code == requests.codes.ok:
        data = []
        col_names = None
        for line in r.text.split('\n'):
            cols = line.rstrip().split('\t')
            if line:
                if col_names:
                    data.append(dict(zip(col_names, cols)))
                else:
                    col_names = cols
        return data
    else:
        return False


def write_json(data, output):
    """Write input data structure to a json file."""
    with open(output, 'w') as fh:
        json.dump(data, fh, indent=4, sort_keys=True)


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='ena-dl',
        conflict_handler='resolve',
        description=(''))
    group1 = parser.add_argument_group('Options', '')
    group1.add_argument('query', metavar="STRING", type=str,
                        help=('ENA accession to query. (Study, Experiment, or '
                              'Run accession)'))
    group1.add_argument('output', metavar="OUTPUT_DIR", type=str,
                        help=('Directory to output downloads to.'))
    group1.add_argument('--quiet', action='store_true', default=False,
                        help='Do not print current status.')
    group1.add_argument('--group_by_experiment', action='store_true',
                        default=False,
                        help='Group runs by experiment accession.')
    group1.add_argument('--group_by_sample', action='store_true',
                        default=False, help='Group runs by sample accession.')
    group1.add_argument('--is_study', action='store_true', default=False,
                        help='Query is a study accession.')
    group1.add_argument('--is_experiment', action='store_true', default=False,
                        help='Query is an experiment accession.')
    group1.add_argument('--is_run', action='store_true', default=False,
                        help='Query is a run accession.')
    group1.add_argument('--ftp', action='store_true', default=False,
                        help='Skip Aspera Connect and use FTP downloads.')
    group1.add_argument('--nextflow', action='store_true', default=False,
                        help='Output instrument model and paired status.')
    group1.add_argument('--debug', action='store_true', default=False,
                        help='Skip downloads, print what will be downloaded.')

    args = parser.parse_args()
    if args.quiet:
        logging.basicConfig(stream=sys.stdout, level=logging.ERROR)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    query = None
    if args.is_study:
        query = 'study_accession={0}'.format(args.query)
    elif args.is_experiment:
        query = 'experiment_accession={0}'.format(args.query)
    elif args.is_run:
        query = 'run_accession={0}'.format(args.query)
    else:
        # Try to guess...
        if args.query[1:3] == 'RR':
            query = 'run_accession={0}'.format(args.query)
        elif args.query[1:3] == 'RX':
            query = 'experiment_accession={0}'.format(args.query)
        else:
            query = 'study_accession={0}'.format(args.query)

    ena_data = get_run_info(query)

    outdir = os.getcwd() if args.output == './' else '{0}'.format(args.output)
    log_stdout('Query: {0}'.format(args.query), quiet=args.quiet)
    log_stdout('Total Runs To Download: {0}'.format(len(ena_data)))

    # FASTQ file names
    runs = None
    is_miseq = False
    is_paired = False
    r1 = None
    r2 = None
    if args.group_by_experiment or args.group_by_sample:
        runs = {}
    for run in ena_data:
        log_stdout('\tWorking on run {0}...'.format(run['run_accession']))

        aspera = run['fastq_aspera'].split(';')
        ftp = run['fastq_ftp'].split(';')
        md5 = run['fastq_md5'].split(';')
        is_paired = True if run['library_layout'] == 'PAIRED' else False
        for i in range(len(aspera)):
            is_r1 = False
            is_r2 = False
            # If run is paired only include *_1.fastq and *_2.fastq, rarely a
            # run can have 3 files.
            # Example:ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1143237
            if is_paired:
                if aspera[i].endswith('_2.fastq.gz'):
                    # Example: ERR1143237_2.fastq.gz
                    is_r2 = True
                elif aspera[i].endswith('_1.fastq.gz'):
                    # Example: ERR1143237_1.fastq.gz
                    pass
                else:
                    # Example: ERR1143237.fastq.gz
                    # Not apart of the paired end read, so skip this file. Or,
                    # its the only fastq file, and its not a paired
                    obs_fq = os.path.basename(aspera[i])
                    exp_fq = '{0}.fastq.gz'.format(run['run_accession'])
                    if (len(aspera) == 1 and obs_fq == exp_fq):
                        is_paired = False
                    else:
                        continue

            # Download Run
            if md5[i] and not args.debug:
                success, fastq = download_fastq(aspera[i], ftp[i], outdir,
                                                md5[i], use_ftp=args.ftp)
                if success:
                    if args.group_by_experiment or args.group_by_sample:
                        name = run["sample_accession"]
                        if args.group_by_experiment:
                            name = run["experiment_accession"]

                        if name not in runs:
                            runs[name] = {'r1': [], 'r2': []}

                            if 'miseq' in run['instrument_model'].lower():
                                is_miseq = True

                        if is_r2:
                            runs[name]['r2'].append(fastq)
                        else:
                            runs[name]['r1'].append(fastq)
                else:
                    print(
                        "Failed to download matching files after 20 attempts "
                        "(10 via Aspera Connect, and 10 via FTP). Please try "
                        "again later or manually download from ENA."
                    )
                    sys.exit()

    # If applicable, merge runs
    if runs and not args.debug:
        for name, vals in runs.items():
            if len(vals['r1']) and len(vals['r2']):
                # Not all runs labled as paired are actually paired...
                if len(vals['r1']) == len(vals['r2']):
                    log_stdout(
                        "\tMerging paired end runs to {0}...".format(name)
                    )
                    r1 = '{0}/{1}_R1.fastq.gz'.format(outdir, name)
                    r2 = '{0}/{1}_R2.fastq.gz'.format(outdir, name)
                    merge_runs(vals['r1'], r1)
                    merge_runs(vals['r2'], r2)
                else:
                    log_stdout("\tMerging single end runs to experiment...")
                    r1 = '{0}/{1}.fastq.gz'.format(outdir, name)
                    merge_runs(vals['r1'], r1)
            else:
                log_stdout("\tMerging single end runs to experiment...")
                r1 = '{0}/{1}.fastq.gz'.format(outdir, name)
                merge_runs(vals['r1'], r1)
        write_json(runs, "{0}/ena-run-mergers.json".format(outdir))
    write_json(ena_data, "{0}/ena-run-info.json".format(outdir))

    if args.nextflow:
        # Assumes grouped by single experiment/sample was downloaded! Mainly
        # used for Staphopia Nextflow pipeline.
        fq2 = "--fq2 {0}".format(r2) if r2 else ""
        is_miseq = "--is_miseq" if is_miseq else ""
        print("--fq1 {0} {1} {2}".format(r1, fq2, is_miseq))
