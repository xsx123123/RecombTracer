#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loguru-based logger initialization.

Loguru is imported at module level with a fallback so that ``import``
succeeds even when loguru is not installed.  An :class:`ImportError` is
raised only when logger functions are actually called.
"""

import datetime
import sys
from pathlib import Path

try:
    from loguru import logger
    _LOGURU_AVAILABLE = True
except ImportError:
    _LOGURU_AVAILABLE = False
    logger = None


def _ensure_loguru() -> None:
    if not _LOGURU_AVAILABLE:
        raise ImportError("loguru is required for logging. Install it with: pip install loguru")


def logger_init(
    logger_name: str = None,
    log_level: str = "INFO",
    more_info: bool = False,
) -> "logger":
    """
    Configure Loguru logger with simple, clean output.
    """
    _ensure_loguru()

    from .configuration import load_default_config
    default = load_default_config()
    cfg_logs = default.get("logs", {})

    log_level = log_level or cfg_logs.get("log_level", "INFO")
    more_info = more_info if more_info is not None else cfg_logs.get("more_info", False)

    logger.remove()

    # Simple console format
    if more_info:
        console_fmt = (
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<cyan>{name}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
        )
    else:
        console_fmt = (
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level: <8}</level> | "
            "<level>{message}</level>"
        )

    # Add Console Handler
    logger.add(
        sys.stderr,
        format=console_fmt,
        level=log_level,
        colorize=True,
    )

    # Add File Handler if name provided
    if logger_name:
        logger.add(
            logger_name,
            format=console_fmt,
            level=log_level,
            colorize=False,
        )
    
    return logger


def logger_generator(
    output_dir: str,
    log_level: str = "INFO",
    more_info: bool = False,
) -> tuple:
    """
    Create a logger, log software metadata, and return ``(logger, output_dir)``.

    Software metadata (name, version, author, email) is read from
    ``recombtracer/config/software.yaml``.

    Parameters
    ----------
    output_dir : str
        Directory where the log file will be written.
    log_level : str
        Minimum logging level.  Default: INFO.
    more_info : bool
        If True, use the verbose format with file/line info.  Default: False.

    Returns
    -------
    tuple
        ``(logger, output_dir)``
    """
    _ensure_loguru()

    # Lazy import to avoid circular dependencies.
    from .configuration import load_software_config, load_default_config

    software = load_software_config()
    default = load_default_config()

    sw = software.get("software", {})
    label = default.get("logs", {}).get("Label", "RecombTracer")

    times = datetime.datetime.now().strftime("%Y-%m-%d:%X:%p")
    logger_name = f"{output_dir}/{label}_{times}.log"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger = logger_init(logger_name, log_level=log_level, more_info=more_info)

    logger.info(f"RecombTracer Author : {sw.get('author', 'unknown')}")
    logger.info(f"RecombTracer Version : {sw.get('version', 'unknown')}")
    logger.info(f"RecombTracer Email : {sw.get('email', '')}")
    logger.info(f"Logger initialized, log file : {logger_name}")
    logger.info(f"RecombTracer Analysis Result : {output_dir}")
    logger.debug(f"Software Full config : {software}")
    logger.debug(f"Default Full config : {default}")
    return logger, output_dir
