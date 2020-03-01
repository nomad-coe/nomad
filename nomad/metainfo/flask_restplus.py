from flask_restplus import fields

from nomad.app.common import RFC3339DateTime

from .metainfo import Section, Quantity, Datetime


def field(quantity: Quantity):
    ''' Returns a flask restplus field with quantity type and shape. '''
    field = None
    if quantity.type == int:
        field = fields.Integer
    elif quantity.type == float:
        field = fields.Float
    elif quantity.type == str:
        field = fields.String
    elif quantity.type == bool:
        field = fields.Boolean
    elif quantity.type == Datetime:
        field = RFC3339DateTime
    else:
        raise NotImplementedError

    result = field(description=quantity.description)

    if len(quantity.shape) == 0:
        return result
    elif len(quantity.shape) == 1:
        return fields.List(result)
    else:
        raise NotImplementedError


def generate_flask_restplus_model(api, section_def: Section):
    return api.model(section_def.name, {
        name: field(quantity)
        for name, quantity in section_def.all_quantities.items()
    })
