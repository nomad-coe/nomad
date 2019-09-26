import json
from optimade.filterparser import LarkParser
from lark import Transformer

from nomad.processing import Upload
from nomad.search import SearchRequest
from nomad.metainfo.optimade import OptimadeStructureEntry


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id

    with published.upload_files.archive_file(calc_id) as f:
        data = json.load(f)

    assert 'OptimadeStructureEntry' in data
    search_result = SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'optimade' in search_result


class ESTransformer(Transformer):

    def __default__(self, tree, children, *args, **kwargs):
        return children[0]

    def and_expr(self, args):
        if len(args) == 1:
            return args[0]

        return dict(op='AND', ops=list(args))

    def or_expr(self, args):
        if len(args) == 1:
            return args[0]

        return dict(op='OR', ops=list(args))

    def not_expr(self, args):
        if len(args) == 1:
            return args[0]

        return dict(op='NOT', ops=list(args))

    def cmp_op(self, args):
        return dict(op=args[1], ops=[args[0], args[2]])

    def list_op(self, args):
        if len(args) == 3:
            return dict(op='HAS', qualifier=args[1], ops=[args[0], args[2]])
        else:
            return dict(op='HAS', ops=[args[0], args[1]])

    def known_op(self, args):
        return dict(op='KNOWN', qualifier=args[1], ops=[args[0]])

    def list(self, args):
        return list(args)

    def tuple(self, args):
        return list(args)

    def predicate(self, args):
        if len(args) == 1:
            return args[0]
        return dict(pred=args[0], op=args[1])

    def quantity(self, args):
        quantity_name = args[0]
        quantity_def = OptimadeStructureEntry.m_section.quantities.get(quantity_name, None)

        if quantity_def is None:
            raise Exception('%s is not a known quantity' % quantity_name)

        return quantity_def.name

    def literal(self, args):
        literal = args[0]

        try:
            int(literal)
        except Exception:
            pass

        try:
            float(literal)
        except Exception:
            pass

        return literal.strip('"')


def test_optimade_parser():
    p = LarkParser(version=(0, 10, 0))
    tree = p.parse('nelements < 3.4e-10  OR elements:elements_ratios HAS ALL "H":>1, "O":>2 AND (elements CONTAINS "H")')
    transformer = ESTransformer()
    result = transformer.transform(tree)
    print(json.dumps(result, indent=2))
