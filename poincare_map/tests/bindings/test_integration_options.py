from action_integration_ext import IntegrationOptions

defaults = {'abs_err': 1e-10,
            'rel_err': 1e-10,
            'dt': 1e-2,
            'coefficient': 1e3}


def test_IntegrationOptions_default_constructor0():
    options = IntegrationOptions()
    assert options.abs_err == defaults['abs_err']
    assert options.rel_err == defaults['rel_err']
    assert options.dt == defaults['dt']
    assert options.coefficient == defaults['coefficient']


def test_IntegrationOptions_default_constructor2():
    options = IntegrationOptions(abs_err=1e-7,
                                 rel_err=1e-9)
    assert options.abs_err == 1e-7
    assert options.rel_err == 1e-9
    assert options.dt == defaults['dt']
    assert options.coefficient == defaults['coefficient']


def test_IntegrationOptions_default_constructor3():
    options = IntegrationOptions(abs_err=1e-7,
                                 rel_err=1e-9,
                                 dt=4e-3)
    assert options.abs_err == 1e-7
    assert options.rel_err == 1e-9
    assert options.dt == 4e-3
    assert options.coefficient == defaults['coefficient']


def test_IntegrationOptions_default_constructor4():
    options = IntegrationOptions(abs_err=1e-7,
                                 rel_err=1e-9,
                                 dt=4e-3,
                                 coefficient=1e5)
    assert options.abs_err == 1e-7
    assert options.rel_err == 1e-9
    assert options.dt == 4e-3
    assert options.coefficient == 1e5


if __name__ == "__main__":

    options = IntegrationOptions(abs_err=1e-12,
                                 rel_err=1e-10)
    print(dir(options))
