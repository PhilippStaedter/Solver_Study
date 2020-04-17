# python script to transfer the observables from the new sbml models to sbml2amici via a dictionary

from petab import sbml


def get_observables(sbml_model, remove=False):
    observables = sbml.assignment_rules_to_dict(
        sbml_model,
        filter_function=sbml_parameter_is_observable,
        remove=remove)
    return observables

def sbml_parameter_is_observable(sbml_parameter):
    return sbml_parameter.getId().startswith('observable_')
