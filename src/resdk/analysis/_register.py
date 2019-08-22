"""Patch ReSDK resources with analysis methods."""
from resdk.analysis.chip_seq import create_trackhub, download_qc, get_bamsplit_qc
from resdk.analysis.utils import link_project

from resdk.resources import Collection, Sample

Collection.create_trackhub = create_trackhub
Collection.download_qc = download_qc
Collection.get_bamsplit_qc = get_bamsplit_qc
Collection.link_project = link_project

Sample.create_trackhub = create_trackhub
Sample.download_qc = download_qc
Sample.get_bamsplit_qc = get_bamsplit_qc
Sample.link_project = link_project
