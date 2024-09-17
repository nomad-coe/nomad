from nomad.metainfo import Quantity, MSection, MEnum


def simple_quantity():
    return Quantity(
        name='test',
        type=str,
        shape=[],
        description="""Sample description""",
        aliases=['alias1', 'alias2'],
    )


def test_quantity():
    """
    The following properties affect the hash of a Quantity:
        name
        aliases
        type
        shape
        unit
        default
        virtual

    """
    q1 = simple_quantity()
    q2 = simple_quantity()

    ref_hash = q1.definition_id

    assert ref_hash == q2.definition_id

    # order of aliases does not matter
    q2.aliases = ['alias2', 'alias1']
    q2.hash(regenerate=True)
    assert ref_hash == q2.definition_id

    # description does not matter
    q2.description = 'Some other text'
    q2.hash(regenerate=True)
    assert ref_hash == q2.definition_id

    # different aliases matter
    q2.aliases = ['alias2', 'alias1', 'alias3']
    # cached so does not change
    assert ref_hash == q2.definition_id
    # regenerate to have different hash
    q2.hash(regenerate=True)
    assert ref_hash != q2.definition_id

    # type matters
    q2 = simple_quantity()
    q2.type = float
    assert ref_hash != q2.definition_id

    # default value matters
    q2 = simple_quantity()
    q2.default = 'default'
    assert ref_hash != q2.definition_id

    # shape matters
    q2 = simple_quantity()
    q2.shape = [1]
    assert ref_hash != q2.definition_id

    # default value matters
    q2 = simple_quantity()
    q2.default = 'default'
    assert ref_hash != q2.definition_id

    # virtual matters
    q2 = simple_quantity()
    q2.virtual = True
    assert ref_hash != q2.definition_id

    q2 = simple_quantity()
    q2.type = MEnum('aad', 'wwa', 'qe')
    q2.hash(regenerate=True)
    assert ref_hash != q2.definition_id

    # order of enum values do not matter
    ref_hash = q2.definition_id
    q2.type = MEnum('wwa', 'qe', 'aad')
    q2.hash(regenerate=True)
    assert ref_hash == q2.definition_id


def test_section():
    class Sample(MSection):
        q1 = simple_quantity()
        q2 = simple_quantity()

    ref_hash = Sample.m_def.definition_id

    del Sample

    class Sample(MSection):  # pylint: disable=function-redefined
        q1 = simple_quantity()
        q2 = simple_quantity()

    # assert equality
    assert ref_hash == Sample.m_def.definition_id

    class Sample(MSection):  # pylint: disable=function-redefined
        q2 = simple_quantity()
        q1 = simple_quantity()

    # order of quantities matters
    assert ref_hash != Sample.m_def.definition_id
