#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import re
import sys

sys.path.append('monovar_src/')

from monovar import main, parse_args

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main(parse_args()))
