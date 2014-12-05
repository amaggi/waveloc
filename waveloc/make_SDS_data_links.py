#!/usr/bin/env python
# encoding: utf-8

import os
import optparse
import glob
import logging
from obspy.core import read


def make_SDS_data_links(datadir, dataglob, outdir):
    """
    Creates an SDS directory structure from a flat directory of data files.
    The SDS structure will contain symlinks to the original data, to avoid
    copying data. If existing links exist, they will be overwritten.

    :param datadir: Directory with original data
    :param dataglob: File glob to select data of interest, e.g. *.BHZ
    :param outdir: Directory in which to put the root of the SDS structure
    """

    data_dir = os.path.abspath(datadir)
    out_dir = os.path.abspath(outdir)
    logging.debug('Data from directory %s' % data_dir)
    logging.debug('SDS data to be put in directory %s' % out_dir)

    all_files = glob.glob(os.path.join(data_dir, dataglob))
    all_files.sort()

    filedict = {}
    for filename in all_files:
        st = read(filename)
        net = st.traces[0].stats.network
        sta = st.traces[0].stats.station
        cha = st.traces[0].stats.channel
        dirid = "%s.%s.%s" % (net, sta, cha)
        if dirid in filedict:
            filedict[dirid].append(filename)
        else:
            filedict[dirid] = [filename]

    for dirid, filelist in filedict.iteritems():
        net = dirid.split('.')[0]
        sta = dirid.split('.')[1]
        cha = dirid.split('.')[2]
        dirname = os.path.join(out_dir, net, sta, "%s.D" % cha)
        try:
            os.makedirs(dirname)
            logging.info("Made directories : %s" % dirname)
        except OSError:
            logging.debug("Directories already exist : %s" % dirname)
            pass

        for my_file in filelist:
            dest_file = os.path.join(dirname, os.path.basename(my_file))
            try:
                os.symlink(my_file, dest_file)
                logging.info("Linked %s" % dest_file)
            except OSError:
                logging.debug("Removing old %s" % dest_file)
                os.remove(dest_file)
                os.symlink(my_file, dest_file)
                logging.info("Linked %s" % dest_file)


if __name__ == '__main__':

    p = optparse.OptionParser()

    p.add_option('--datadir', action='store', help="data directory")
    p.add_option('--dataglob', action='store', help="data glob")
    p.add_option('--outdir', action='store', help="output directory")
    p.add_option('--debug', action='store_true',
                 help="turn on debugging output")

    (options, arguments) = p.parse_args()

    if options.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(levelname)s : %(asctime)s : %(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(levelname)s : %(asctime)s : %(message)s')

    make_SDS_data_links(options.datadir, options.dataglob, options.outdir)
