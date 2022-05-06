# -*- coding: utf-8 -*-
"""This module defines functions and constants related to MOTIF."""

__all__ = ["prositeToRegEx"]


def prositeToRegEx(pattern: str) -> str:
    """Change PROSITE Motif to python regular expression.

    Args:
        pattern (str): PROSITE Motif

    Returns:
        str: regular expression
    """
    pattern = (
        pattern.replace("{", "[^")
        .replace("}", "]")
        .replace("(", "{")
        .replace(")", "}")
        .replace("-", "")
        .replace("x", ".")
        .replace(">", "$")
        .replace("<", "*")
    )
    return pattern
