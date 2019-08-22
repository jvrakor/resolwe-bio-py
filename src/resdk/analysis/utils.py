"""Utility functions for ReSDK."""

import getpass
import os

from pexpect import pxssh
from six.moves import input

SSH_HOSTNAME = 'mhgcp-h00.grid.bcm.edu'
DATA_FOLDER_PATH = '/storage/genialis/bcm.genialis.com/data_by_id/'

__all__ = ('link_project')

def get_resolwe(*resources):
    """Return resolwe object used in given resources.
    Raise an error if there is more than one.
    """
    resolwes = {res_obj.resolwe for res_obj in resources}
    if len(resolwes) != 1:
        raise TypeError('All input objects must be from the same `Resolwe` connection.')

    return list(resolwes)[0]

def link_project(resource, genome_name, path, output_table=''):
    links = [
                {'type': 'data:alignment:bam:bowtie2:', 'field': 'bam', 'subfolder': 'bams'},
                {'type': 'data:alignment:bam:bowtie2:', 'field': 'bai', 'subfolder': 'bams'},
                {'type': 'data:alignment:bam:bwaaln:', 'field': 'bam', 'subfolder': 'bams'},
                {'type': 'data:alignment:bam:bwaaln:', 'field': 'bai', 'subfolder': 'bams'},
                {'type': 'data:alignment:bam:hisat2:', 'field': 'bam', 'subfolder': 'bams'},
                {'type': 'data:alignment:bam:hisat2:', 'field': 'bai', 'subfolder': 'bams'},
                {'type': 'data:chipseq:callpeak:macs14:','field':'peaks_bed','subfolder': 'macs14'},
                {'type': 'data:chipseq:callpeak:macs2:','field':'narrow_peaks','subfolder': 'macs2'},
                {'type': 'data:chipseq:rose2:','field':'ALL','subfolder': 'rose2'},
                {'type': 'data:cufflinks:cuffquant:','field':'ALL','subfolder': 'cufflinks/cuffquant'},
                {'type': 'data:expressionset:cuffnorm:','field':'ALL','subfolder': 'cufflinks/cuffnorm'},
            ]

    # This must be a dict, so the reference doesn't breake even if it
    # is assigned in sub-function.
    ssh = {'connection': None}
    def _create_local_link(src, dest):
        dest_dir = os.path.dirname(dest)
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        if os.path.isfile(dest):
            os.remove(dest)
        os.symlink(src, dest)
    def _create_ssh_link(src, dest):
        if ssh['connection'] is None:
            print('Credentials for connection to {}:'.format(SSH_HOSTNAME))
            username = input('username: ')
            password = getpass.getpass('password: ')
            ssh['connection'] = pxssh.pxssh()
            ssh['connection'].login(SSH_HOSTNAME, username, password)
        dest_dir = os.path.dirname(dest)
        ssh['connection'].sendline('mkdir -p "{}"'.format(dest_dir))
        ssh['connection'].sendline('ln -sf "{}" "{}"'.format(src, dest))

    data_table = [['FILE_PATH','UNIQUE_ID','GENOME','NAME','BACKGROUND','ENRICHED_REGION','ENRICHED_MACS','COLOR','FASTQ_FILE']]

    resolwe = get_resolwe(resource)

    print('Linking results...')
    for link in links:
        for data in resource.data.filter(status='OK', type=link['type']):

            if link['field'] == 'ALL':
                files_filter = {}
            else:
                files_filter = {'field_name': link['field']}

            for file_name in data.files(**files_filter):

                if link['field'] == 'bam':
                    bam_path = os.path.join(path, 'bams/')
                    unique_ID = str(data.id)
                    genome = str(genome_name).upper()
                    name = str(data.name)
                    enriched_region = 'NONE'
                    color = '0,0,0'

                    reads = data.sample.get_reads()
                    fastq = os.path.join(DATA_FOLDER_PATH, str(reads.id), reads.output['fastq'][0]['file'])

                    macs = data.sample.data.filter(type = 'data:chipseq:callpeak:macs14')
                    macs2 = data.sample.data.filter(type = 'data:chipseq:callpeak:macs2')
                    background = 'NONE'
                    if macs:
                        enriched_macs ='{:05}_{}'.format(
                           macs[0].id,
                           macs[0].output['peaks_bed']['file']
                        )
                        if 'control' in macs[0].input:
                            background = resolwe.data.get(macs[0].input['control']).name
                    elif macs2:
                        enriched_macs ='{:05}_{}'.format(
                           macs2[0].id,
                           macs2[0].output['narrow_peaks']['file'].replace('.narrowPeak.gz', '.bed')
                        )
                        if 'control' in macs2[0].input:
                            background = resolwe.data.get(macs2[0].input['control']).name
                    else:
                        enriched_macs = 'NONE'

                    new_line = [
                       bam_path, 
                       unique_ID, 
                       genome, 
                       name, 
                       background,
                       enriched_region, 
                       enriched_macs, 
                       color,
                       fastq,
                    ]
                    
                    data_table.append(new_line)

                file_path = os.path.join(DATA_FOLDER_PATH, str(data.id), file_name)                         
                link_name = '{:05}_{}'.format(
                    data.id,
                    file_name,
                )
                link_path = os.path.join(path, link['subfolder'], link_name)
                if os.path.isfile(file_path):
                    _create_local_link(file_path, link_path)
                else:
                    _create_ssh_link(file_path, link_path)


    sep = '\t'
    if output_table == '':
        output_table = 'data_table.txt'


    fh_out = open(output_table, 'w')
    if len(sep) == 0:
        for i in data_table:
            fh_out.write(str(i) + '\n')
    else:
        for line in data_table:
            fh_out.write(sep.join([str(x) for x in line]) + '\n')

    fh_out.close()


    if ssh['connection'] is not None:
        ssh['connection'].logout()
