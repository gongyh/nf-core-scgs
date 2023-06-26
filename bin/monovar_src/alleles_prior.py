"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""


def get_prior_matrix(p):
    prior_matrix = {
        ('AA', 'T'): p, ('AA', 'G'): p, ('AA', 'C'): p, ('AA', 'A'): 1 - 3 * p,
        ('TT', 'T'): 1 - 3 * p, ('TT', 'G'): p, ('TT', 'C'): p, ('TT', 'A'): p,
        ('GG', 'T'): p, ('GG', 'A'): p, ('GG', 'C'): p, ('GG', 'G'): 1 - 3 * p,
        ('CC', 'C'): 1 - 3 * p, ('CC', 'G'): p, ('CC', 'T'): p, ('CC', 'A'): p,
        ('AC', 'T'): p, ('AC', 'G'): p,
        ('AC', 'C'): (1 - 2 * p) / 2, ('AC', 'A'): (1 - 2 * p) / 2,
        ('AG', 'T'): p, ('AG', 'C'): p,
        ('AG', 'G'): (1 - 2 * p) / 2, ('AG', 'A'): (1 - 2 * p) / 2,
        ('AT', 'C'): p, ('AT', 'G'): p,
        ('AT', 'T'): (1 - 2 * p) / 2, ('AT', 'A'): (1 - 2 * p) / 2,
        ('CG', 'T'): p, ('CG', 'A'): p,
        ('CG', 'C'): (1 - 2 * p) / 2, ('CG', 'G'): (1 - 2 * p) / 2,
        ('CT', 'A'): p, ('CT', 'G'): p,
        ('CT', 'T'): (1 - 2 * p) / 2, ('CT', 'C'): (1 - 2 * p) / 2,
        ('GT', 'C'): p, ('GT', 'A'): p,
        ('GT', 'T'): (1 - 2 * p) / 2, ('GT', 'G'): (1 - 2 * p) / 2,
    }
    return prior_matrix
