#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration loading utilities for RecombTracer.

All YAML configs are stored inside the installed package under
recombtracer/config/ and accessed via importlib.resources so that
paths remain valid regardless of the current working directory.
"""

from importlib import resources
from typing import Any

import yaml


def _load_yaml(config_name: str) -> dict[str, Any]:
    """Load a YAML file from recombtracer.config by basename (no extension)."""
    config_path = resources.files("recombtracer.config") / f"{config_name}.yaml"
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_software_config() -> dict[str, Any]:
    """Load software metadata from software.yaml."""
    return _load_yaml("software")


def load_default_config() -> dict[str, Any]:
    """Load default analysis parameters from default.yaml."""
    return _load_yaml("default")
