__author__ = "Jwely"


def merge_dicts(*args):

    result = {}
    for arg in args:
        result.update(arg)
    return result