"""
IS Element Extractor Module

Extract IS element DNA sequences and flanking regions from MGEfinder output.
Annotate ORFs and noncoding regions using Prodigal.
"""

from .extractor import ISExtractor, extract_is_elements, extract_batch
from .orf_finder import annotate_is_element, annotate_is_elements

__version__ = "0.3.0"
__all__ = [
    "ISExtractor",
    "extract_is_elements",
    "extract_batch",
    "annotate_is_element",
    "annotate_is_elements"
]
