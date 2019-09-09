#! /usr/bin/env python3
"""
usage: ena-dl [-h] [--aspera STRING] [--aspera_key STRING]
              [--aspera_speed STRING] [--is_study] [--is_experiment]
              [--is_run] [--group_by_experiment] [--group_by_sample]
              [--outdir OUTPUT_DIR] [--max_retry INT] [--ftp_only] [--silent]
              [--debug] [--version]
              ACCESSION

ena-dl - Download FASTQs from ENA

optional arguments:
  -h, --help            show this help message and exit

Required Options:

  ACCESSION             ENA accession to query. (Study, Experiment, or Run
                        accession)

Aspera Connect Options:
  --aspera STRING       Path to the Aspera Connect tool "ascp" (Default:
                        "which ascp")
  --aspera_key STRING   Path to Aspera Connect private key, if not given,
                        guess based on ascp path
  --aspera_speed STRING
                        Speed at which Aspera Connect will download. (Default:
                        100M)

Query Related Options:
  --is_study            Query is a Study.
  --is_experiment       Query is an Experiment.
  --is_run              Query is a Run.
  --group_by_experiment
                        Group Runs by experiment accession.
  --group_by_sample     Group Runs by sample accession.

Helpful Options:
  --outdir OUTPUT_DIR   Directory to output downloads to. (Default: ./)
  --max_retry INT       Maximum times to retry downloads (Default: 10)
  --ftp_only            FTP only downloads.
  --silent              Only critical errors will be printed.
  --debug               Skip downloads, print what will be downloaded.
  --version             show program's version number and exit
"""
PROGRAM = "ena-dl"
VERSION = "1.0.0"
import logging
import os
import subprocess

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
        return [open(output, 'w'), f'{redirect} {output}']
    else:
        return [subprocess.PIPE, '']


def onfinish_handler(cmd, out, err, returncode):
    out = f'\n{out}' if out else ''
    err = f'\n{err}' if err else ''
    if returncode != 0:
        logging.error(f'COMMAND: {cmd}')
        logging.error(f'STDOUT: {out}')
        logging.error(f'STDERR: {err}')
        logging.error('END\n')
        raise RuntimeError(err)
    else:
        logging.info(f'COMMAND: {cmd}')
        logging.info(f'STDOUT: {out}')
        logging.info(f'STDERR: {err}')
        logging.info(f'END\n')
        return [out, err]


def byte_to_string(b):
    return b.decode("utf-8") if b else ''


def run_command(cmd, cwd=os.getcwd(), stdout=False, stderr=False):
    """Execute a single command and return STDOUT and STDERR."""
    stdout, stdout_str = output_handler(stdout)
    stderr, stderr_str = output_handler(stderr, redirect='2>')
    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, cwd=cwd)
    out, err = p.communicate()
    return onfinish_handler(
        f'{" ".join(cmd)} {stdout_str} {stderr_str}',
        byte_to_string(out), byte_to_string(err), p.returncode
    )


def log_stdout(message):
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


def download_fastq(fasp, ftp, outdir, md5, aspera, max_retry=10, ftp_only=False):
    """Download FASTQ from ENA using Apera Connect."""
    import time
    success = False
    retries = 0
    fastq = f'{outdir}/{os.path.basename(fasp)}'

    if not os.path.exists(fastq):
        if not os.path.isdir(outdir):
            run_command(['mkdir', '-p', outdir])

        while not success:
            if ftp_only:
                log_stdout(f'\t\tFTP download attempt {retries + 1}')
                run_command(['wget', '--quiet', '-O', fastq, ftp])
            else:
                log_stdout(f'\t\tAspera Connect download attempt {retries + 1}')
                run_command([
                    aspera['ascp'], '-QT', '-l', aspera['speed'], '-P33001',
                    '-i',aspera['private_key'], f'era-fasp@{fasp}', outdir
                ])

            if md5sum(fastq) != md5:
                retries += 1
                if os.path.exists(fastq):
                    os.remove(fastq)
                if retries > max_retry:
                    if not ftp_only:
                        ftp_only = True
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
    url = f'{ENA_URL}&query="{query}"&fields={",".join(FIELDS)}'
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
    import json
    with open(output, 'w') as fh:
        json.dump(data, fh, indent=4, sort_keys=True)


def parse_query(query, is_study, is_experiment, is_run):
    "Parse user query, to determine search field value."
    if is_study:
        return f'study_accession={query}'
    elif is_experiment:
        return f'experiment_accession={query}'
    elif args.is_run:
        return f'run_accession={query}'
    else:
        # Try to guess...
        if query[1:3] == 'RR':
            return f'run_accession={query}'
        elif query[1:3] == 'RX':
            return f'experiment_accession={query}'
        else:
            return f'study_accession={query}'


def check_aspera(ascp, private_key, speed):
    "Verify Aspera Connect is available, not if it works."
    error_message = None
    if not os.path.exists(ascp):
        error_message = f'cannot access "{ascp}": No such file or directory'
    else:
        if private_key:
            # User provided path to private key
            if not os.path.exists(private_key):
                error_message = f'cannot access "{private_key}": No such file or directory'
        else:
            # Try to guess private key path, based on ascp path
            key_path = os.path.dirname(ascp).replace('/bin', '/etc')
            private_key = f'{key_path}/asperaweb_id_dsa.openssh'
            if not os.path.exists(private_key):
                error_message = f'cannot access "{private_key}": No such file or directory'

    if error_message:
        logging.ERROR(f'Aspera Related Error: {error_message}')
        sys.exit(1)
    else:
        return {'ascp': ascp, 'private_key': private_key, 'speed': speed}



if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog='ena-dl',
        conflict_handler='resolve',
        description=(f'{PROGRAM} (v{VERSION}) - Download FASTQs from ENA')
    )
    group1 = parser.add_argument_group('Required Options', '')
    group1.add_argument('query', metavar="ACCESSION", type=str,
                        help=('ENA accession to query. (Study, Experiment, or '
                              'Run accession)'))

    group2 = parser.add_argument_group('Aspera Connect Options')
    group2.add_argument('--aspera', metavar="STRING", type=str,
                        help='Path to the Aspera Connect tool "ascp" (Default: "which ascp")')
    group2.add_argument(
        '--aspera_key', metavar="STRING", type=str,
        help='Path to Aspera Connect private key, if not given, guess based on ascp path'
    )
    group2.add_argument('--aspera_speed', metavar="STRING", type=str, default="100M",
                        help='Speed at which Aspera Connect will download. (Default: 100M)')

    group3 = parser.add_argument_group('Query Related Options')
    group3.add_argument('--is_study', action='store_true', help='Query is a Study.')
    group3.add_argument('--is_experiment', action='store_true', help='Query is an Experiment.')
    group3.add_argument('--is_run', action='store_true', help='Query is a Run.')
    group3.add_argument('--group_by_experiment', action='store_true',
                        help='Group Runs by experiment accession.')
    group3.add_argument('--group_by_sample', action='store_true',
                        help='Group Runs by sample accession.')

    group4 = parser.add_argument_group('Helpful Options')
    group4.add_argument('--outdir', metavar="OUTPUT_DIR", type=str, default='./',
                        help=('Directory to output downloads to. (Default: ./)'))
    group4.add_argument('--max_retry', metavar="INT", type=int, default=10,
                        help='Maximum times to retry downloads (Default: 10)')
    group4.add_argument('--ftp_only', action='store_true', help='FTP only downloads.')
    group4.add_argument('--silent', action='store_true',
                        help='Only critical errors will be printed.')
    group4.add_argument('--debug', action='store_true',
                        help='Skip downloads, print what will be downloaded.')
    group4.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(stream=sys.stdout,
                        level=logging.ERROR if args.silent else logging.INFO)

    aspera = check_aspera(args.aspera, args.aspera_key, args.aspera_speed) if args.aspera else None
    if not aspera:
        log_stdout("Aspera Connect not available, using FTP")
        args.ftp_only = True

    outdir = os.getcwd() if args.outdir == './' else f'{args.outdir}'
    query = parse_query(args.query, args.is_study, args.is_experiment, args.is_run)

    # Start Download Process
    ena_data = get_run_info(query)
    log_stdout(f'Query: {args.query}')
    log_stdout(f'Total Runs To Download: {len(ena_data)}')
    runs = {} if args.group_by_experiment or args.group_by_sample else None
    for run in ena_data:
        log_stdout(f'\tWorking on run {run["run_accession"]}...')
        fasp = run['fastq_aspera'].split(';')
        ftp = run['fastq_ftp'].split(';')
        md5 = run['fastq_md5'].split(';')
        for i in range(len(fasp)):
            is_r2 = False
            # If run is paired only include *_1.fastq and *_2.fastq, rarely a
            # run can have 3 files.
            # Example:ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/007/ERR1143237
            if run['library_layout'] == 'PAIRED':
                if fasp[i].endswith('_2.fastq.gz'):
                    # Example: ERR1143237_2.fastq.gz
                    is_r2 = True
                elif fasp[i].endswith('_1.fastq.gz'):
                    # Example: ERR1143237_1.fastq.gz
                    pass
                else:
                    # Example: ERR1143237.fastq.gz
                    # Not apart of the paired end read, so skip this file. Or,
                    # its the only fastq file, and its not a paired
                    obs_fq = os.path.basename(fasp[i])
                    exp_fq = f'{run["run_accession"]}.fastq.gz'
                    if (len(fasp) != 1 and obs_fq != exp_fq):
                        continue

            # Download Run
            if md5[i] and not args.debug:
                success, fastq = download_fastq( fasp[i], ftp[i], outdir, md5[i], aspera,
                                                max_retry=args.max_retry, ftp_only=args.ftp_only)
                if success:
                    if args.group_by_experiment or args.group_by_sample:
                        name = run["sample_accession"]
                        if args.group_by_experiment:
                            name = run["experiment_accession"]

                        if name not in runs:
                            runs[name] = {'r1': [], 'r2': []}

                        if is_r2:
                            runs[name]['r2'].append(fastq)
                        else:
                            runs[name]['r1'].append(fastq)
                else:
                    logging.error(
                        f'Download files after {args.max_retry} attempts (Aspera and/or '
                        'FTP). Please try again later or manually from ENA.'
                    )
                    sys.exit(1)

    # If applicable, merge runs
    if runs and not args.debug:
        for name, vals in runs.items():
            if len(vals['r1']) and len(vals['r2']):
                # Not all runs labled as paired are actually paired.
                if len(vals['r1']) == len(vals['r2']):
                    log_stdout(f'\tMerging paired end runs to {name}...')
                    merge_runs(vals['r1'], f'{outdir}/{name}_R1.fastq.gz')
                    merge_runs(vals['r2'], f'{outdir}/{name}_R2.fastq.gz')
                else:
                    log_stdout('\tMerging single end runs to experiment...')
                    merge_runs(vals['r1'], f'{outdir}/{name}.fastq.gz')
            else:
                log_stdout('\tMerging single end runs to experiment...')
                merge_runs(vals['r1'], f'{outdir}/{name}.fastq.gz')
        write_json(runs, f'{outdir}/ena-run-mergers.json')
    write_json(ena_data, f'{outdir}/ena-run-info.json')
