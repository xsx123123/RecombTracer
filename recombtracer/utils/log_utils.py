#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loguru-based logger initialization with Hydra/OmegaConf integration.

Heavy dependencies (loguru, omegaconf, hydra) are imported lazily inside
functions so that importing this module does not fail when they are absent.
"""

import datetime


def logger_init(logger_name: str, cfg: "DictConfig") -> "Logger":
    """
    Configure and return a Loguru logger instance.

    Parameters
    ----------
    logger_name : str
        Path to the log file.
    cfg : DictConfig
        Configuration object with a ``logs`` section.

    Returns
    -------
    Logger
        Configured Loguru logger.
    """
    from loguru import logger
    from omegaconf import DictConfig

    logger.remove()
    if cfg.logs.more_info:
        fmt = (
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> |"
            "{level.icon}|"
            "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
            "<level>{message}</level>"
        )
        logger.add(logger_name, format=fmt, level=cfg.logs.log_level, colorize=True, serialize=False)
        logger.add(logger_name, format=fmt, level=cfg.logs.log_level, colorize=True, serialize=False)
    else:
        fmt = (
            "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
            "<level>{level}</level> | <level>{message}</level>"
        )
        logger.add(logger_name, format=fmt, level=cfg.logs.log_level, colorize=True, serialize=False)
        logger.add(logger_name, format=fmt, level=cfg.logs.log_level, colorize=True, serialize=False)
    return logger


def logger_generator(cfg: "DictConfig") -> tuple:
    """
    Generate a logger and output directory based on Hydra runtime config.

    Parameters
    ----------
    cfg : DictConfig
        Full OmegaConf configuration.

    Returns
    -------
    tuple
        (logger, output_dir)
    """
    from hydra.core.hydra_config import HydraConfig
    from loguru import logger

    times = datetime.datetime.now().strftime("%Y-%m-%d:%X:%p")
    output_dir = HydraConfig.get().runtime.output_dir
    logger_name = f"{output_dir}/{cfg.logs.project_id}_{times}.log"

    # Remove the default logger to avoid duplicate logs
    logger = logger_init(logger_name, cfg)
    logger.info(f"RecombTracer Author:{cfg.RecombTracer.Author}")
    logger.info(f"RecombTracer Version:{cfg.RecombTracer.Version}")
    logger.info(f"RecombTracer Email:{cfg.RecombTracer.Email}")
    logger.info(f"Logger initialized, log file: {logger_name}")
    logger.info(f"RecombTracer Analysis workspace: {HydraConfig.get().runtime.cwd}")
    logger.info(f"RecombTracer Analysis Result: {output_dir}")
    logger.debug(f"Full config: {cfg}")
    return logger, output_dir
