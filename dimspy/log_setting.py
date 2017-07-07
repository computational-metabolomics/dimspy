#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
log_setting: provide functions to config logging system.

author(s): Albert
origin: Dec. 11, 2014

"""

import logging, sys


class LevelFormatter(logging.Formatter):
    def __init__(self, fmt_dict=None, *args, **kwargs):
        self._fmtdict = fmt_dict if fmt_dict is not None else {
            logging.DEBUG: ' D> %(message)s',
            logging.INFO: '    %(message)s',
            logging.WARN: '\n ?? WARNING: %(message)s \n',
            logging.ERROR: '\n !! ERROR: %(message)s \n',
        }
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record):
        assert self._fmtdict.has_key(record.levelno), 'unknown logging level [%s]' % str(record.levelno)
        self._fmt = self._fmtdict[record.levelno]
        return logging.Formatter.format(self, record)


def set_logging(level=logging.INFO, print_formatter=None, log_file=None, log_mode='w', log_formatter=None,
                exception_format=None):
    logger = logging.getLogger()
    logger.handlers = []  # remove existing handlers
    logger.setLevel(level)

    # log to terminal
    assert print_formatter is None or isinstance(print_formatter,
                                                 logging.Formatter), 'print formatter is not a Formatter object'
    pfmt = LevelFormatter() if print_formatter is None else print_formatter
    chdl = logging.StreamHandler()
    chdl.setFormatter(pfmt)
    logger.addHandler(chdl)

    # log to file
    if log_file is not None:
        assert log_formatter is None or isinstance(log_formatter,
                                                   logging.Formatter), 'log file formatter is not a Formatter object'
        lfmt = logging.Formatter(
            '%(asctime)-15s : [%(levelname)s] >> %(message)s') if log_formatter is None else log_formatter
        fhdl = logging.FileHandler(log_file, mode=log_mode)
        fhdl.setFormatter(lfmt)
        logger.addHandler(fhdl)

    # handle exception / assertion
    assert exception_format is None or isinstance(exception_format, str), 'exception format is not a string object'

    def _excepthook(*args):
        logging.getLogger().error(
            ('%(class)s: %(message)s' if exception_format is None else exception_format) %
            {'class': args[1].__class__.__name__, 'message': args[1].message},
            exc_info=False if level > logging.DEBUG else args
        )

    sys.excepthook = _excepthook

    return logger  # to add extra handlers
