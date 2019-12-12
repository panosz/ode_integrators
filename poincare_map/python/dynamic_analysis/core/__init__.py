from dynamic_analysis.core._core import *


def _integration_options_details(options):
    """
    Get options details in a nice format.

    Returns
        class_name: str
        details: dict of attribute - value pairs
    """
    class_name = type(options).__name__
    details = dict([('abs_err', options.abs_err),
                    ('rel_err', options.rel_err),
                    ('dt', options.dt),
                    ('coefficient', options.coefficient)])
    return class_name, details


def _integration_options_str(self):
    class_name, details = _integration_options_details(self)
    details_strings = [f'{name}: {value}' for name, value in details.items()]
    self_strings = [f'{class_name} object:'] + details_strings
    return '\n\t'.join(self_strings)


def _integration_options_repr(self):
    class_name, details = _integration_options_details(self)
    arg_string = ', '.join(f'{name}={value}'
                           for name, value in details.items())
    return class_name + '(' + arg_string + ')'


IntegrationOptions.__str__ = _integration_options_str
IntegrationOptions.__repr__ = _integration_options_repr
