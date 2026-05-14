#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import platform
import importlib.metadata as importlib_metadata

def _print_info_dict(info):
    for k, v in info.items():
        print(f"{k:>12}: {v}")

def _get_sys_info():
    return {
        "python": sys.version.replace("\n", " "),
        "executable": sys.executable,
        "machine": platform.platform(),
    }

def _get_deps_info(deps):
    def get_ver(pkg):
        try:
            return importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            return None
    return {pkg: get_ver(pkg) for pkg in deps}

def show_versions(project_name=None, deps=None, extras=None):
    """
    打印调试用的版本与系统信息，风格类似 sklearn.show_versions()。
    - project_name: 你的项目名（用于打印标题），如 "myproj"
    - deps: 关注的依赖包列表，如 ["numpy","scipy","pandas"]
    - extras: 额外的键值信息 dict，如 {"CUDA_VISIBLE_DEVICES": os.getenv("CUDA_VISIBLE_DEVICES")}
    """
    deps = deps or []
    title = f"{project_name} deps:" if project_name else "Dependencies:"
    print(title)
    _print_info_dict(_get_deps_info(deps))

    if extras:
        print("\nExtras:")
        _print_info_dict(extras)

    print("\nSystem:")
    _print_info_dict(_get_sys_info())
