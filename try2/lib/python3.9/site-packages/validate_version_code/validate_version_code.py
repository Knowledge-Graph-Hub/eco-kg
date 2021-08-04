import re
import os
from glob import glob


def validate_version_code(version_code: str) -> bool:
    """Return a boolean representing if given version code is valid.
        version_code:str, the version code to validate.
    """
    return bool(re.compile(r"\d+\.\d+\.\d+").match(version_code))


def extract_version_code(path: str = None) -> str:
    """Return the version code from given project."""
    if path is None:
        return extract_version_code(os.path.dirname(glob("*/__version__.py")[0]))
    with open("{path}/__version__.py".format(path=path), "r") as f:
        return re.compile(r"(\d+\.\d+\.\d+)").findall(f.read())[0]
