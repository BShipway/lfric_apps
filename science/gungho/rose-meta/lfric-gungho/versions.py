import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn31_pmsl_solver(MacroUpgrade):
    # Upgrade macro for PMSL alternative solver by B. Shipway

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_pmsl_solver"

    def upgrade(self, config, meta_config=None):
        # Add PMSL solver type (default: jacobi for backwards compatibility)
        self.add_setting(config, ["namelist:physics", "pmsl_solver"], "'jacobi'")
        # Add SOR omega parameter (default 0.0 = adaptive/auto-computed at runtime)
        self.add_setting(config, ["namelist:physics", "pmsl_omega"], "0.0")
        return config, self.reports


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""
