Command Line Interface
======================

.. code-block:: console

    $ dimspy --help

    Executing dimspy version 2.0.0.
    usage: __main__.py [-h]
                       {process-scans,replicate-filter,align-samples,blank-filter,sample-filter,remove-samples,mv-sample-filter,merge-peaklists,get-peaklists,get-average-peaklist,hdf5-pm-to-txt,hdf5-pls-to-txt,create-sample-list,unzip,licenses}
                       ...

    Python package to process DIMS data

    positional arguments:
      {process-scans,replicate-filter,align-samples,blank-filter,sample-filter,remove-samples,mv-sample-filter,merge-peaklists,get-peaklists,get-average-peaklist,hdf5-pm-to-txt,hdf5-pls-to-txt,create-sample-list,unzip,licenses}
        process-scans       Process scans and/or stitch SIM windows.
        replicate-filter    Filter irreproducible peaks from technical replicate
                            peaklists.
        align-samples       Align peaklists across samples.
        blank-filter        Filter peaks across samples that are present in the
                            blank samples.
        sample-filter       Filter peaks based on certain reproducibility and
                            sample class criteria.
        remove-samples      Remove sample(s) from a peak matrix object or list of
                            peaklist objects.
        mv-sample-filter    Filter samples based on the percentage of missing
                            values.
        merge-peaklists     Merge peaklists from multiple lists of peaklist or
                            peak matrix objects.
        get-peaklists       Get peaklists from a peak matrix object.
        get-average-peaklist
                            Get an average peaklist from a peak matrix object.
        hdf5-pm-to-txt      Write HDF5 output (peak matrix) to text format.
        hdf5-pls-to-txt     Write HDF5 output (peak lists) to text format.
        create-sample-list  Create a sample list from a peak matrix object or list
                            of peaklist objects.
        unzip               Extract files from zip file
        licenses            Show licenses DIMSpy and RawFileReader

    optional arguments:
      -h, --help            show this help message and exit


.. code-block:: console

    $ dimspy process-scans --help

    Executing dimspy version 2.0.0b1.
    usage: __main__.py process-scans [-h] -i source -o OUTPUT [-l FILELIST] -m
                                     {median,mean,mad,noise_packets} -s
                                     SNR_THRESHOLD [-p PPM] [-n MIN_SCANS]
                                     [-a MIN_FRACTION] [-d RSD_THRESHOLD] [-k]
                                     [-r RINGING_THRESHOLD]
                                     [-e start end scan_type]
                                     [-x start end scan_type] [-z start end]
                                     [-u REPORT] [-b BLOCK_SIZE] [-c NCPUS]

    optional arguments:
      -h, --help            show this help message and exit
      -i source, --input source
                            Directory (*.raw, *.mzml or tab-delimited peaklist
                            files), single *.mzml/*.raw file or zip archive
                            (*.mzml only)
      -o OUTPUT, --output OUTPUT
                            HDF5 file to save the peaklist objects to.
      -l FILELIST, --filelist FILELIST
                            Tab-delimited file that include the name of the data
                            files (*.raw or *.mzml) and meta data. Column names:
                            filename, replicate, batch, injectionOrder,
                            classLabel.
      -m {median,mean,mad,noise_packets}, --function-noise {median,mean,mad,noise_packets}
                            Select function to calculate noise.
      -s SNR_THRESHOLD, --snr-threshold SNR_THRESHOLD
                            Signal-to-noise threshold
      -p PPM, --ppm PPM     Mass tolerance in Parts per million to group peaks
                            across scans / mass spectra.
      -n MIN_SCANS, --min_scans MIN_SCANS
                            Minimum number of scans required for each m/z range or
                            event.
      -a MIN_FRACTION, --min-fraction MIN_FRACTION
                            Minimum fraction a peak has to be present. Use 0.0 to
                            not apply this filter.
      -d RSD_THRESHOLD, --rsd-threshold RSD_THRESHOLD
                            Maximum threshold - relative standard deviation
                            (Calculated for peaks that have been measured across a
                            minimum of two scans).
      -k, --skip-stitching  Skip the step where (SIM) windows are 'stitched' or
                            'joined' together. Individual peaklists are generated
                            for each window.
      -r RINGING_THRESHOLD, --ringing-threshold RINGING_THRESHOLD
                            Ringing
      -e start end scan_type, --include-scan-events start end scan_type
                            Scan events to select. E.g. 100.0 200.0 sim or 50.0
                            1000.0 full
      -x start end scan_type, --exclude-scan-events start end scan_type
                            Scan events to select. E.g. 100.0 200.0 sim or 50.0
                            1000.0 full
      -z start end, --remove-mz-range start end
                            M/z range(s) to remove. E.g. 100.0 102.0 or 140.0
                            145.0.
      -u REPORT, --report REPORT
                            Summary/Report of processed mass spectra
      -b BLOCK_SIZE, --block-size BLOCK_SIZE
                            The size of each block of peaks to perform clustering
                            on.
      -c NCPUS, --ncpus NCPUS
                            Number of central processing units (CPUs).
