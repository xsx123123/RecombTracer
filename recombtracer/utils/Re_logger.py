#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import datetime
from loguru import logger
from omegaconf import DictConfig
from hydra.core.hydra_config import HydraConfig
# ------ logger initialization ------ #
def logger_init(logger_name:str,cfg:DictConfig):
    """
    配置并返回一个 Loguru 日志记录器实例。
    这个函数设置了日志的格式、颜色和输出位置。
    """
    logger.remove()
    if cfg.logs.more_info:
        logger.add(sys.stdout,format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> |"
                "{level.icon}|"
                "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
                "<level>{message}</level>",
                level=cfg.logs.log_level,colorize=True,serialize=False)
        logger.add(logger_name,
                format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> |"
                "{level.icon}|"
                "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
                "<level>{message}</level>",
                level=cfg.logs.log_level,colorize=True,serialize=False)
    else:
        logger.add(sys.stdout,format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                "<level>{level}</level> | <level>{message}</level>",
                level=cfg.logs.log_level,colorize=True,serialize=False)
        logger.add(logger_name,
                format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
                "<level>{level}</level> | <level>{message}</level>",
                level=cfg.logs.log_level,colorize=True,serialize=False)
    return logger
 
def logger_generator(cfg:DictConfig): 
    """
    配置并返回一个 Loguru 日志记录器实例。
    这个函数设置了日志的格式、颜色和输出位置。
    """ 
    times = datetime.datetime.now().strftime("%Y-%m-%d:%X:%p")
    logger_name = f"{HydraConfig.get().runtime.output_dir}/{cfg.logs.project_id}_{times}.log"
    # Remove the default logger to avoid duplicate logs
    logger = logger_init(logger_name,cfg)
    logger.info(f"RecombTracer Author:{cfg.RecombTracer.Author}")
    logger.info(f"RecombTracer Version:{cfg.RecombTracer.Version}")
    logger.info(f"RecombTracer Email:{cfg.RecombTracer.Email}")
    logger.info(f"Logger initialized, log file: {logger_name}")
    logger.info(f"RecombTracer Analysis workspace: {HydraConfig.get().runtime.cwd}")
    logger.info(f"RecombTracer Analysis Result: {HydraConfig.get().runtime.output_dir}")
    logger.debug(f"Full config: {cfg}")
    output_dir = HydraConfig.get().runtime.output_dir
    return logger,output_dir
# ------ logger initialization ------ #