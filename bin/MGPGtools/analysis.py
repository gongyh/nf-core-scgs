import logging
from MGPGtools.utils.taxonomy import Taxonomy
from MGPGtools.static import Stat
from MGPGtools.extern.odgi import Odgi
from MGPGtools.search import Search
from MGPGtools.describe import Describe
from MGPGtools.tree import Tree
from MGPGtools.core import Core
from MGPGtools.exceptions import stat_args_judgment, viz_args_judgment


class Analysis(object):
    def __init__(self):
        self.logger = logging.getLogger("timestamp")
        self.warnings = logging.getLogger("warnings")

    def statistic(self, options):
        stat_args_judgment(options)
        rank_index = Taxonomy.rank_index[options.rank]
        stat = Stat(options.db, rank_index, options.name, options.outdir)
        stat.write()

    def viz(self, options):
        viz_args_judgment(options)
        odgi = Odgi(
            options.db, options.name, options.outdir, options.outName, options.threads
        )
        odgi.viz(options)

    def search(self, options):
        search = Search(options)
        search.ref_info()

    def describe(self, options):
        describe = Describe(options)
        describe.describe_gfa()

    def core(self, options):
        core = Core(options)
        if options.coreGenes and options.fasta:
            core.contigCompleteness()
        else:
            core.staticCoreGene()

    def tree(self, options):
        tree = Tree(options)
        if options.gene:
            tree.drawGeneTree()
        if options.genesFile:
            if options.fasta:
                tree.drawTreeWithContig()
            else:
                tree.drawSpeciesTree()
