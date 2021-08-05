"""Conversion between rdflib and PyLD data formats for compatibility.

The rdflib and PyLD libraries use different in-memory data structures
so it is not trivial to use the powerful RDF support of rdflib in conjuction
with the JSON-LD processor implementation of PyLD. This code is designed
to provide a bridge between rdflib and PyLD in-memory formats, avoiding
the need to serialize and then re-parse the data.

rdflib: see <http://rdflib.readthedocs.io/en/stable/apidocs/>

PyLD: see <https://github.com/digitalbazaar/pyld>

FIXME - this is limited in that it assumes the URIRefs in the rdflib graph
do not need expansion with a namespaceManager.
"""

from rdflib import Graph, URIRef, Literal, BNode
from pyld import jsonld


__version__ = '0.0.1'


def _rdflib_term_to_pyld_term(term):
    """Convert rdflib term to a PyLD term.

    See:
    <http://rdflib.readthedocs.io/en/stable/apidocs/rdflib.html#module-rdflib.term>
    and
    <https://github.com/digitalbazaar/pyld/blob/master/lib/pyld/jsonld.py#L1645>a

    Return value is a PyLD term which is represented as a hash where:
      `type` in {'IRI', 'blank node', 'literal'}
      `value` is string of IRI, BNoode or literal
      `language` optional language code for literal (if present `datatype` will
          default to RDF.langString)
      `datatype` optional data type URI for literal, else None
    """
    t = {'type': 'IRI',
         'value': term.toPython()}
    if (isinstance(term, BNode)):
        t['type'] = 'blank node'
    elif (isinstance(term, Literal)):
        t['type'] = 'literal'
        if (term.language):
            t['datatype'] = str(term.datatype) if (term.datatype) else \
                'http://www.w3.org/1999/02/22-rdf-syntax-ns#langString'
            t['language'] = str(term.language)
        else:
            t['datatype'] = str(term.datatype) if (term.datatype) else None
    return(t)


def _pyld_term_to_rdflib_term(term):
    """Convert a PyLD term to a rdflib term.

    Inverse of rdflib_term_to_pyld_term()
    """
    term_type = term.get('type', None)
    if (term_type == 'IRI'):
        return URIRef(term['value'])
    elif (term_type == 'blank node'):
        return BNode(term['value'])
    elif (term_type == 'literal'):
        return Literal(term['value'])
    else:
        raise ValueError("Bad pyld term type '%s'" % (term_type))


def _pyld_dataset_from_rdflib_graph(graph):
    """Get a PyLD dataset from an rdflib graph.

    Returns the contents of the input graph as a list if triples where
    each triple is represented as a dict with term entries for 'subject',
    'predicate', and 'object'.
    """
    g = []
    for s, p, o in graph:
        triple = {}
        triple['subject'] = _rdflib_term_to_pyld_term(s)
        triple['predicate'] = _rdflib_term_to_pyld_term(p)
        triple['object'] = _rdflib_term_to_pyld_term(o)
        g.append(triple)
    return g


def _fix_type_null(j):
    """In-place recursive fix of problem with type: null in JSON-LD.

    "pred": [
      {
        "@type": null,
        "@value": "value"
      }

    to

    "pred": "value"

    FIXME -- how can this be avoided in the first place?
    """
    if (isinstance(j, list)):
        for j2 in j:
            _fix_type_null(j2)
    elif (isinstance(j, dict)):
        if (len(j) == 2 and '@type' in j and '@value' in j):
            if (j['@type'] is None):
                # case to fix if @type is None, simply remove '@type'
                j.pop('@type')
        else:
            for k in j:
                _fix_type_null(j[k])


def _rdflib_graph_from_pyld_dataset(dataset):
    """Get an rdflib graph from a PyLD dataset."""
    g = Graph()
    for triple in dataset:
        s = _pyld_term_to_rdflib_term(triple['subject'])
        p = _pyld_term_to_rdflib_term(triple['predicate'])
        o = _pyld_term_to_rdflib_term(triple['object'])
        g.add([s, p, o])
    return g


def pyld_jsonld_from_rdflib_graph(graph):
    """Get PyLD JSON-LD object from and rdflib input graph."""
    default_graph = _pyld_dataset_from_rdflib_graph(graph)
    json = jsonld.from_rdf({'@default': default_graph})
    _fix_type_null(json)
    return(json)


def rdflib_graph_from_pyld_jsonld(json):
    """Get rdflib graph from PyLD JSON-LD object."""
    default_graph = jsonld.to_rdf(json)['@default']
    return(_rdflib_graph_from_pyld_dataset(default_graph))
