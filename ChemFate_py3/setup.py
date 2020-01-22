import sys
from cx_Freeze import setup, Executable

setup(
    name = "ChemFate",
    version = "1.0",
    description = "ChemFate model",
    executables = [Executable("run_GUI.py")])