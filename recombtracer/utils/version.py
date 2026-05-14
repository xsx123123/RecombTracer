#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import platform
import importlib.metadata as importlib_metadata

# Load custom functions from the configuration module
# Compatible with both package import and direct execution
if __name__ == "__main__":
    from pathlib import Path
    _project_root = Path(__file__).resolve().parent.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))
    from recombtracer.utils.configuration import load_software_config
    from recombtracer.utils.print_utils import print_info_table, print_section_panel, print_in_columns
else:
    from .configuration import load_software_config
    from .print_utils import print_info_table, print_section_panel, print_in_columns

def _get_sys_info():
    return {
        "Python": sys.version.split("|")[0].strip(),
        "Executable": sys.executable,
        "Machine": platform.platform(),
    }

def _get_deps_info(deps):
    def get_ver(pkg):
        try:
            return importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            return "[red]Not Found[/red]"
    return {pkg: get_ver(pkg) for pkg in deps}

def show_versions(project_name=None, deps=None, extras=None):
    """
    Print debug version and system information with rich beauty.
    """
    deps = deps or []
    tables = []
    
    # 1. Dependencies Table
    deps_info = _get_deps_info(deps)
    tables.append(print_info_table(
        title=f"[bold blue]{project_name}[/bold blue] Deps",
        data=deps_info,
        column_names=["Package", "Version"],
        header_style="bold magenta",
        return_table=True
    ))

    # 2. Software Metadata Table
    if extras:
        tables.append(print_info_table(
            title="[bold green]Software Info[/bold green]",
            data=extras,
            column_names=["Attribute", "Value"],
            header_style="bold yellow",
            border_style="green",
            return_table=True
        ))
    
    # Print tables in columns (horizontal layout)
    print_in_columns(tables)

    # 3. System Info (keep it below as it's usually wider)
    print_info_table(
        title="[bold cyan]System Environment[/bold cyan]",
        data=_get_sys_info(),
        column_names=["Component", "Details"],
        header_style="bold white",
        border_style="cyan"
    )

if __name__ == "__main__":
    # get software configuration from the config file
    _software_conf = load_software_config()
    sw = _software_conf.get("software", {})
    
    # Define interesting dependencies
    deps_to_check = ["numpy", "pandas", "scipy", "cyvcf2", "pyyaml", "rich"]
    
    # Show versions with full metadata
    show_versions(
        project_name=sw.get("app_name", "RecombTracer"),
        deps=deps_to_check,
        extras={
            "Version": sw.get("version", "unknown"),
            "Author": sw.get("author", "unknown"),
            "Email": sw.get("email", "unknown"),
            "Description": sw.get("description", "")
        }
    )
