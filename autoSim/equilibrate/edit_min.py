from .make_metaGroups import make_metaGroups
from .rewrite_settings import rewrite_settings

import json
import os

#####
def edit_min(root):

    with open(f'{root}/parm/meta_atomTypes.json') as f:
        types = json.load(f)
    print(types)

    make_metaGroups(types)

    #rewrite_settings('system.in.settings')

#####
if __name__ == '__main__':
    edit_min()
    os.system('lmp -in in.min')
