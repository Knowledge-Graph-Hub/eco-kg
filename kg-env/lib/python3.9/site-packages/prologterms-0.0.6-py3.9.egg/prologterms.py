"""
Prolog Terms

All Prolog data structures are called terms. A term is either:

 - A constant, which can be either an atom or a number (represented using python types)
 - A variable (represented using a Var object)
 - A compound term (represented using a Term object)

In prolog, structures such as lists and tuples are also terms. Here these can be
directly represented by the analogous python structure

"""

import re

__author__ = "Chris Mungall <cjmungall@lbl.gov>"
__version__ = "0.0.6"

class Program(object):
    """
    A program is a collection of terms

    note: in future this may also allow for creation of modules
    """
    def __init__(self, *terms):
        self.terms = list(terms)

class Term(object):
    """
    Compound terms, which consist of a predicate plus one or more arguments, each of which is a term
    """
    def __init__(self, pred, *args):
        self.pred = pred
        self.args = list(args)
        self.comments = []

    def add_comment(self, c):
        """
        Terms can be annotated with comments
        """
        self.comments.append(c)
        
    # override "<=" to mean same as prolog ":-"
    def __le__(self, body):
        return Rule(self, body)

    def __lt__(self, body):
        return Term('<', self, body)
    def __gt__(self, body):
        return Term('>', self, body)
    def __ge__(self, body):
        return Term('>=', self, body)
    
    # override "==" to mean same as prolog "="
    def __eq__(self, body):
        return Term('=', self, body)
    
    # override "!=" to mean same as prolog "\="
    def __ne__(self, body):
        return Term('\=', self, body)

    def __sub__(self, body):
        return Term('-', self, body)
    def __add__(self, body):
        return Term('+', self, body)
    def __mul__(self, body):
        return Term('*', self, body)
    def __truediv__(self, body):
        return Term('/', self, body)
    def __pow__(self, body):
        return Term('**', self, body)
    
    # override % for comments
    def __mod__(self, comment):
        self.comments = [comment]
        return self

    # unary
    def __neg__(self):
        return Term('-', self)

    # use ~ for NOT
    def __invert__(self):
        return Term('\+', self)
    

class Rule(Term):
    """
    A Rule is a compound term with 2 args, where the predicate is ":-"
    the first arg is a compound term representing the head (consequent)
    the second arg is a tuple of compound terms representing the body (antecedents)
    """
    def __init__(self, head, body):
        """
        Arguments
        ---------
        head: Term
        body: tuple of Terms
        """
        self.pred = ':-'
        self.args = [head, body]
        self.comments = []

class Var(Term):
    """
    Represents a prolog variable.
    These should use a leading uppercase.
    TODO: allow other naming conventions
    """
    def __init__(self, name):
        self.name = name
        self.comments = []

class TermGenerator(object):
    """
    An object for conveniently generating compound terms.
    After creating a TermGenerator G, call G.PRED(ARG1, ..., )
    to generate a compound term
    """
    def __getattr__(self, pred):
        def method(*args):
            return Term(pred, *args)
        return method

class Renderer(object):
    pass

class PrologRenderer(Renderer):
    """
    Renders internal prolog term or program representations as strings that can be
    fed directly to a prolog engine
    """
    def render(self, t):
        """
        Renders a prolog term, program or python structure holding these.

        Note that single terms are rendered without a closing ".". Include
        terms in a Program object to render with "."s
        """
        s=""
        if isinstance(t, Program):
            s = "{}".format("".join([self.render(x) + ".\n" for x in t.terms]))
        elif isinstance(t, Rule):
            head = self.render(t.args[0])
            bodyt = t.args[1]
            body = ""
            if isinstance(bodyt, tuple):
                body = "{}".format(", ".join([self.render(e) for e in bodyt]))
            else:
                body = self.render(bodyt)
            s = "{} :-\n    {}".format(head, body)
        elif isinstance(t, Var):
            s = "{}".format(t.name)
        elif isinstance(t, Term):
            if len(t.args) == 0:
                s = t.pred
            else:
                s = "{}({})".format(t.pred, ", ".join([self.render(a) for a in t.args]))
        elif isinstance(t, list):
            s = "[{}]".format(", ".join([self.render(e) for e in t]))
        elif isinstance(t, tuple):
            s = "({})".format(", ".join([self.render(e) for e in t]))
        elif isinstance(t, str):
            s = "{}".format(t)
            if not re.match(r"^[a-z]\w*$", s):
                s = s.replace("'","\\'").replace("\n","\\n")
                s = "'{}'".format(s)
        else:
            s = "{}".format(t)
        # add comments
        if isinstance(t, Term) and len(t.comments) > 0:
            cmt = "".join(['% {}\n'.format(c) for c in t.comments])
            s = cmt + s
        return s

class SExpressionRenderer(Renderer):
    """
    Renders internal prolog term or program representations as S-Expression strings that can be
    fed directly to a prolog engine such as kanren
    """
    def render(self, t):
        """
        Renders a prolog term, program or python structure holding these.
        """
        s=""
        if isinstance(t, Program):
            s = "{}".format(" ".join([self.render(x) for x in t.terms]))
        elif isinstance(t, Rule):
            head = self.render(t.args[0])
            bodyt = t.args[1]
            body = ""
            if isinstance(bodyt, tuple):
                body = "{}".format(" ".join([self.render(e) for e in bodyt]))
            else:
                body = self.render(bodyt)
            s = "(<= {} {})".format(head, body)
        elif isinstance(t, Var):
            s = "?{}".format(t.name)
        elif isinstance(t, Term):
            s = "({} {})".format(t.pred, " ".join([self.render(a) for a in t.args]))
        elif isinstance(t, list):
            s = "(list {})".format(" ".join([self.render(e) for e in t]))
        elif isinstance(t, tuple):
            s = "({})".format(" ".join([self.render(e) for e in t]))
        elif isinstance(t, str):
            s = "{}".format(t)
            if not re.match(r"^[a-z]\w*$", s):
                s = s.replace("'","\\'").replace("\n","\\n")
                s = "'{}'".format(s)
        else:
            s = "{}".format(t)
        # add comments
        if isinstance(t, Term) and len(t.comments) > 0:
            cmt = "".join(['; {}\n'.format(c) for c in t.comments])
            s = cmt + s
        return s

    
