#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import os
import re
import sys

script_dir = os.path.abspath( os.path.dirname( __file__ ) )

sys.path.append(script_dir + '/monovar_src/')

from monovar import main, parse_args

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main(parse_args()))
