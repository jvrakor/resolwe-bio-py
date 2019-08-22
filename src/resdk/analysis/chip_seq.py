"""Chip Seq analysis."""
import os
import pandas as pd

from resdk.resources.utils import is_sample, is_data

from collections import defaultdict
from datetime import datetime
from shutil import copyfile

__all__ = ('create_trackhub', 'download_qc', 'get_bamsplit_qc')

def get_macs2(resource):
    """Return list of ``bed`` objects on the sample."""
    return list(resource.data.filter(type='data:chipseq:callpeak:macs2'))

def get_samples(resource):
    """Get the list of samples from given resources.
    Get the list of samples with:
        * use recursion if given resource is a list
        * return the resource if it is already the sample
        * call ResolweQuery object named `samples` (if exists) and return
          the result
    """
    error_msg = ("Resource should be sample, have `samples` query, be list of multiple "
                 "resources or be data object with not empty `sample` property.")
    if isinstance(resource, list):
        samples = []
        for res in resource:
            samples.extend(get_samples(res))
        return samples

    elif is_data(resource):
        if not resource.sample:
            raise TypeError(error_msg)

        return [resource.sample]

    elif is_sample(resource):
        return [resource]

    elif hasattr(resource, 'samples'):
        return resource.samples

    else:
        raise TypeError(error_msg)

def get_resource_collection(resource, fail_silently=True):
    """Get id of the collection to which resource belongs.
    If resource does not belong to any collection or collection cannot
    be determined uniquely, ``None`` is returned of ``LookupError`` is
    raised (if ``fail_silently`` is set to ``True``).
    """
    if is_collection(resource):
        return resource.id

    elif hasattr(resource, 'collection'):
        return resource.collection.id

    if isinstance(resource, list):
        collections_ids = [get_resource_collection(item) for item in resource]
        if len(set(collections_ids)) == 1:
            return collections_ids[0]

    if fail_silently:
        return None

    raise LookupError('Collection id cannot be determined uniquely.')

def get_resolwe(resources):
    """Return resolwe object used in given resources.
    Raise an error if there is more than one.
    """
    resolwes = {res_obj.resolwe for res_obj in resources}
    if len(resolwes) != 1:
        raise TypeError('All input objects must be from the same `Resolwe` connection.')

    return list(resolwes)[0]

def merge_qc_reports(qc_data):
    """Merge ChIP-seq QC reports."""
    data = defaultdict(dict)
    for report_file, report_type, sample_name in qc_data:
        data[sample_name][report_type] = report_file

    prepeak_list = []
    postpeak_list = []
    for sample in data:
        if 'prepeak' in data[sample]:
            prepeak = pd.read_csv(data[sample]['prepeak'], sep='\t')
            prepeak.index = [sample]
            prepeak_list.append(prepeak)
        if 'postpeak' in data[sample]:
            postpeak = pd.read_csv(data[sample]['postpeak'], sep='\t')
            postpeak.index = [sample]
            postpeak_list.append(postpeak)

        if prepeak_list and postpeak_list:
            prepeaks = pd.concat(prepeak_list)
            postpeaks = pd.concat(postpeak_list)
            report = pd.merge(prepeaks, postpeaks, left_index=True, right_index=True, how='outer')
        elif prepeak_list:
            report = pd.concat(prepeak_list)
        else:
            report = pd.concat(postpeak_list)

    report.to_csv('QC_report.txt', sep='\t', na_rep='N/A', index_label='SAMPLE_NAME', float_format='%.3f')

def download_qc(resource, download_dir=None):
    """Download and merge ChIP(ATAC)-seq QC reports.

    This method downloads and merges all QC reports found on MACS2 data
    objects that are connected to samples in the resource.

    """
    if not isinstance(resource, list):
        resource = [resource]

    resolwe = get_resolwe(resource)

    qc_data = []
    for single_resource in resource:
        for sample in get_samples(single_resource):
            macs2 = get_macs2(sample)
            if not macs2:
                continue
            if len(macs2) > 1:
                raise ValueError("There can be only one macs2 object on a sample.")
            macs2 = macs2[0]

            case_prepeak_field = 'case_prepeak_qc'
            case_postpeak_field = 'chip_qc'
            background_prepeak_field = 'control_prepeak_qc'
            case_prepeak_file = None
            case_postpeak_file = None
            background_prepeak_file = None
            if case_prepeak_field in macs2.output:
                macs2.download(field_name=case_prepeak_field, download_dir=download_dir)
                case_prepeak_file = macs2.output[case_prepeak_field]['file']
                qc_data.append((case_prepeak_file, 'prepeak', macs2.sample.name))
            if case_postpeak_field in macs2.output:
                macs2.download(field_name=case_postpeak_field, download_dir=download_dir)
                case_postpeak_file = macs2.output[case_postpeak_field]['file']
                qc_data.append((case_postpeak_file, 'postpeak', macs2.sample.name))
            if background_prepeak_field in macs2.output:
                macs2.download(field_name=background_prepeak_field, download_dir=download_dir)
                background_prepeak_file = macs2.output[background_prepeak_field]['file']
                background_name = resolwe.data.get(macs2.input['control']).sample.name
                qc_data.append((background_prepeak_file, 'prepeak', background_name))

    merge_qc_reports(qc_data)

    for report_file, _, _ in set(qc_data):
        os.remove(report_file)

def get_bamsplit_qc(resource, download_dir=None):
    """Get ChIP-Rx/hybrid mapping metrics.

    Returns overall mapping percantage and the percentages of reads mapped to each genome.
    """
    if not isinstance(resource, list):
        resource = [resource]

    qc_data = []
    for single_resource in resource:
        secondary_bams = single_resource.data.filter(type='data:alignment:bam:secondary')
        for secondary_bam in secondary_bams:
            secondary_file = secondary_bam.output['stats']['file']
            secondary_bam.download(field_name='stats', download_dir=download_dir)
            with open(secondary_file) as infile:
                secondary_count = float(infile.readline().split(' ')[0])
            original_bam = secondary_bam.sample.data.filter(type='data:alignment:bam:bwaaln')[0]
            original_file = original_bam.output['stats']['file']
            original_bam.download(field_name='stats', download_dir=download_dir)
            with open(original_file) as infile:
                original_count = float(infile.readline().split(' ')[0])
                infile.readline()
                infile.readline()
                infile.readline()
                mapped_count = float(infile.readline().split(' ')[0])
            qc_data.append((
                secondary_bam.sample.name,
                str(round(100 * mapped_count / original_count, 2)),
                str(round(100 * (1 - secondary_count / mapped_count), 2)),
            ))
            os.remove(original_file)
            os.remove(secondary_file)

    with open('bamsplit_qc.txt', 'w') as outfile:
        for line in qc_data:
            outfile.write('\t'.join(line) + '\n')

def make_public_trackhub(file_list, name_list, genome, analysis_name, project_folder, hub_name='',
                         hub_short_lab='', hub_long_lab='', email='', file_type='bigWig',
                         col=['0,0,0'],
                         cyverse_url="https://de.cyverse.org/anon-files/iplant/home/linlabbcm"):
    """Create trackhub and upload it to cyverse.

    Write trackhub files with cyverse paths, upload all data, and return url for public trackhub.
    If trackhub already exists can supply path and data will be uploaded to cyverse and adjusted
    and url returned.

    :param list file_list: list of files to be uploaded (bigWig, bigBed...)
    :param list name_list: list of sample names that correspond to the file_list (same order)
    :param str genome: genome build (hg18, hg19, mm9, mm10, rn4, rn6...)
    :param str analysis_name: name of the created trackhub
    :param str project_folder: folder in which trackhub folder will get created
    :param str hub_name: trackhub name in .hub.txt file
    :param str hub_short_lab: trackhub short label in .hub.txt file
    :param str hub_long_lab: trackhub long label in .hub.txt file
    :param str email: email in .hub.txt file
    :param str file_type: file type (bigWig, bigBed, bigNarrowPeak, BAM...)
    :param list col: list of colors that correspond to file_list (same order); if list has only
        one element then all the tracks will be of that color

    """
    if not isinstance(col, list):
        raise TypeError("Parameter 'col' must be a list.")

    if not col:
        raise ValueError("Parameter 'col' must have at least one element.")

    local_track_dir = '{}_publicTrackFiles'.format(os.path.join(project_folder, analysis_name))
    os.mkdir(local_track_dir)

    track_db_file = os.path.join(local_track_dir, 'trackDb.txt')
    hub_file = "{}.hub.txt".format(os.path.join(local_track_dir, analysis_name))
    genomes_file = "{}.genomes.txt".format(os.path.join(local_track_dir, analysis_name))
    cyverse_path = "{}/{}/".format(cyverse_url, analysis_name)  # This is a url join so os.join() is no good

    os.system("imkdir {}".format(analysis_name))

    # Write out all the fields for hub file
    with open(hub_file, 'w') as hub:
        hub.write("hub {}\n".format(hub_name))
        hub.write("shortLabel {}\n".format(hub_short_lab))
        hub.write("longLabel {}\n".format(hub_long_lab))
        hub.write("genomesFile {}{}.genomes.txt\n".format(cyverse_path, analysis_name))
        hub.write("email {}\n".format(email))

    # Write out genome file
    with open(genomes_file, 'w') as genomes:
        genomes.write("genome {}\n".format(genome.lower()))
        genomes.write("trackDb {}trackDb.txt\n".format(cyverse_path))

    if len(col) == 1:
        col = len(file_list) * col

    # Write out track_db_file
    with open(track_db_file, 'w') as track_db:
        for track_file, name, color in zip(file_list, name_list, col):
            print(name)
            copyfile(track_file, os.path.join(local_track_dir, os.path.basename(track_file)))
            track_db.write("track {}\n".format(name))
            track_db.write("bigDataUrl {}{}\n".format(cyverse_path, os.path.basename(track_file)))
            track_db.write("shortLabel {}\n".format(name))
            track_db.write("longLabel {}\n".format(name))
            track_db.write("type {}\n".format(file_type))
            track_db.write("color {}\n".format(color))
            track_db.write("\n")

    # Upload everything to the correct directory on cyverse
    cmd = "iput -K -f {}/* {}".format(local_track_dir, analysis_name)
    os.system(cmd)

    # Finally generate the url to load trackhub and give to collabs
    hubURL = "{}{}".format(cyverse_path, hub_file.split('/')[-1])
    ucscURL = "http://genome.ucsc.edu/cgi-bin/hgTracks?db={}&hubUrl={}".format(
        genome.lower(),
        hubURL
    )
    print(hubURL)
    print(ucscURL)

def create_trackhub(resource, file_dir=None, name=None, background=False):
    """."""
    genialis_path = '/storage/genialis/bcm.genialis.com/data_by_id'
    file_list = []
    name_list = []
    genomes = []
    project_folder = file_dir or os.getcwd()
    download = False

    if not isinstance(resource, list):
        resource = [resource]

    resolwe = get_resolwe(resource)

    case_bigwig_field = 'treat_pileup_bigwig'
    control_bigwig_field = 'control_lambda_bigwig'
    for single_resource in resource:
        for sample in get_samples(single_resource):
            macs2 = get_macs2(sample)
            if not macs2:
                continue
            if len(macs2) > 1:
                raise ValueError("There can be only one macs2 object on a sample.")
            macs2 = macs2[0]
            genomes.append(macs2.output['build'])

            name_list.append(sample.name)
            if background and 'control' in macs2.input:
                name_list.append(resolwe.data.get(macs2.input['control']).sample.name)

            if 'bcm' in str(resolwe):
                case_bigwig_file = os.path.join(
                    genialis_path,
                    str(macs2.id),
                    macs2.output[case_bigwig_field]['file'],
                )
                file_list.append(case_bigwig_file)
                if background and 'control' in macs2.input:
                    control_bigwig_file = os.path.join(
                        genialis_path,
                        str(macs2.id),
                        macs2.output[control_bigwig_field]['file'],
                    )
                    file_list.append(control_bigwig_file)
            else:
                download = True
                macs2.download(field_name=case_bigwig_field, download_dir=project_folder)
                case_bigwig_file = os.path.join(
                    project_folder,
                    macs2.output[case_bigwig_field]['file'],
                )
                file_list.append(case_bigwig_file)
                if background and 'control' in macs2.input:
                    macs2.download(field_name=control_bigwig_field, download_dir=project_folder)
                    control_bigwig_file = os.path.join(
                        project_folder,
                        macs2.output[control_bigwig_field]['file'],
                    )
                    file_list.append(control_bigwig_file)


    if len(set(genomes)) > 1:
        raise TypeError("Samples must have the same genome build.")

    genome = genomes[0]
    analysis_name = name if name else 'Track_{}'.format(datetime.now().isoformat())

    print(file_list)
    print(name_list)
    print(analysis_name)

    make_public_trackhub(
        file_list,
        name_list,
        genome,
        analysis_name,
        project_folder,
        hub_name=analysis_name,
        hub_short_lab='',
        hub_long_lab='',
        email='',
    )

    if download:
        for f in file_list:
            os.remove(f)
