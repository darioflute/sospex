import json
from os.path import dirname

with open(dirname(__file__) + '/version.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']