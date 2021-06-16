"""Conversion between rdflib and PyLD data formats for compatibility.

Imports two functions from rdflib_pyld_compat.convert :

pyld_jsonld_from_rdflib_graph(graph)
    Get PyLD JSON-LD object from and rdflib input graph.

rdflib_graph_from_pyld_jsonld(json)
    Get rdflib graph from PyLD JSON-LD object.
"""

from .convert import pyld_jsonld_from_rdflib_graph, rdflib_graph_from_pyld_jsonld

__version__ = '0.1.1'
