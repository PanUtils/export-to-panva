#!/usr/bin/python

# Author: Astrid van den Brandt
# Code adapted from: https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
# Start date v1: 13-8-2024
# WUR - Bioinformatics
# Functionality:
# Script is used to pre-process PanTools pangenomes for variant analysis in PanVA.

# Imports:
import logging

# Create logger for console messages
class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"
    datefmt="%Y-%m-%d %H:%M:%S"

    FORMATS = {
        logging.DEBUG: grey + "🛠️  " + format + reset,
        logging.INFO: grey + "ℹ️  " + format + reset,
        logging.WARNING: yellow + "⚠️  " + format + reset,
        logging.ERROR: red + "❌ " + format + reset,
        logging.CRITICAL: bold_red + "🚨 " + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

logger = logging.getLogger("export-to-panva")
logger.setLevel(logging.INFO)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

ch.setFormatter(CustomFormatter())

logger.addHandler(ch)