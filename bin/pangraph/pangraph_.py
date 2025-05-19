import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np
import networkx as nx
import math
from networkx import NetworkXNoPath
from .utils import getContigsAdjacency
from .utils import getContigsAdjacency_v2
from .utils import generate_fasta_from_dict
from .utils import buildOverlapEdge
from .utils import read_contigs2dict
from .utils import write_fasta, max_common_subsequence, similarity_sequence
from .utils import append_strand, append_strand_undirected, append_strand_reverse
from .utils import top_prev, recent_large_prev, recent_large_after, get_node_length, vote_sign


class PanGraph:
    def __init__(self, sample_info, gene_info, gene_position, grades=None):
        if sample_info is None:
            return None
        else:
            # input paramters
            self.sample_info = sample_info
            self.gene_info = gene_info
            self.gene_position = gene_position
            self.strand = {}
            # self.grades = grades or {}
            self.multiplicity_bool = False

            # computed parameters
            self.gene2cluster_dict = {}
            self.n_clusters = len(self.gene_info.index)
            self.n_samples = len(self.sample_info.index)
            for i in range(self.n_clusters):
                self.gene2cluster_dict[self.gene_info.iloc[i, 0]] = self.gene_info.iloc[i, 2]
            self.head_contig = {}  # the first gene in the contig
            self.tail_contig = {}  # the last gene in the contig
            self.longhead_contig = {}  # the first gene in the contig
            self.longtail_contig = {}  # the last gene in the contig

    def map_edge_fn(self, i, j, N=10000431):
        # map an edge (i, j) to a number.
        return i * N + j

    def compute_number_nucleotides(self, gene_contigs=None):
        # compute number of nucleotides in the contig
        n_nucleo = 0
        for gene in gene_contigs:
            gc = gene.split("@")
            n_nucleo += int(gc[-2])
        return n_nucleo

    def flatten(self, l):
        if len(l) == 1:
            return l[0]
        else:
            return [item for sublist in l for item in sublist]

    def remove_duplicate(self, your_list):
        return [v for i, v in enumerate(your_list) if i == 0 or v != your_list[i - 1]]

    def get_node_coverage(self, node_id):
        return float(node_id.split("_")[5])

    def get_multiplicity(self, node_id):
        return round(float(node_id.split("_")[5][:-1]) / self.basecoverage)

    def compute_multiplicity(self, sample_df):
        gene_position_sub = sample_df.copy()
        nodes_list = list(gene_position_sub.iloc[:, 1].values)
        nodes_len = [int(node.split("_")[3]) for node in nodes_list]
        nodes_coverage = [float(node.split("_")[5]) for node in nodes_list]
        gene_position_sub["length"] = nodes_len
        gene_position_sub["coverage"] = nodes_coverage
        gene_position_sub = gene_position_sub.sort_values(by="length", ascending=False)
        basecoverage = np.median(gene_position_sub["coverage"][:5])
        expected_node_coverage = [
            min(9, round(gene_position_sub.iloc[i, 4] / basecoverage)) for i in range(len(nodes_coverage))
        ]
        self.node_multiplicity = {}
        for i in range(len(nodes_coverage)):
            self.node_multiplicity[gene_position_sub.iloc[i, 1]] = expected_node_coverage[i]
        # return(self.node_multiplicity)

    def get_pangraph_cost(self, source_id, target_id):
        ### Compute the cost between two contigs in the pangraph
        return_val = -1.0
        # if self.multiplicity_bool:
        #     source_id = source_id[:-3]
        #     target_id = target_id[:-3]
        path_weight_vec = []
        for tail_gene in self.longtail_contig[source_id]:
            for head_gene in self.longhead_contig[target_id]:
                source_id_0 = "C-" + str(tail_gene)
                target_id_0 = "C-" + str(head_gene)
                weight_p = 0.0
                try:
                    p = nx.shortest_path(self.H, source=source_id_0, target=target_id_0)
                except NetworkXNoPath:
                    p = None
                if p is not None:
                    for node_p_idx in range(len(p) - 1):
                        weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                    if len(p) >= 2:
                        weight_p = weight_p / float((len(p) - 1) * (len(p) - 1))
                        path_weight_vec.append(weight_p)
        if len(path_weight_vec) > 0:
            return_val = np.mean(path_weight_vec)
        return return_val

    def remove_unsupported_pangraph(self, adj_list_assembly):
        ### Disconnect the path if the two nodes does not support by a pangraph
        print("Refine the scaffolds")
        adj_list_assembly_new = {}
        contig_id = 0
        for key in adj_list_assembly:
            path = adj_list_assembly[key]
            L = []
            for i in range(len(path)):
                if (
                    self.get_multiplicity(path[i]) <= 1
                    and (path[i][:-1] in self.head_contig)
                    and (get_node_length(path[i]) > 1000)
                ):
                    L.append(i)
            break_point = []
            for i in range(len(L) - 1):
                if self.get_pangraph_cost(path[L[i]][:-1], path[L[i + 1]][:-1]) <= 3.0:
                    break_point.append(L[i])
            start = 0
            for bp in break_point:
                new_path = []
                for j in range(start, bp + 1):
                    new_path.append(path[j])
                if (j + 1 < len(path)) and self.get_multiplicity(path[j + 1]) >= 6:
                    new_path.append(path[j + 1])
                adj_list_assembly_new[key[:-6] + "_c" + str(contig_id)] = new_path
                contig_id = contig_id + 1
                start = bp + 1
            new_path = []
            if start > 0:
                for j in range(break_point[-1] + 1, len(path)):
                    new_path.append(path[j])
            else:
                for j in range(len(path)):
                    new_path.append(path[j])
            if len(new_path) > 0:
                adj_list_assembly_new[key[:-6] + "_c" + str(contig_id)] = new_path
                contig_id = contig_id + 1
        return adj_list_assembly_new

    def get_value_long_contigs(self, edge_df0, source_id, target_id):
        return_val = -1.0
        df_res = edge_df0[((edge_df0["source"] == source_id) & (edge_df0["target"] == target_id))]
        if len(df_res.index) >= 1:
            return_val = float(df_res.iloc[0, 2])
        if return_val < 0.2 * self.n_samples / 4.0:
            print("compute long range dependancy")
            if self.multiplicity_bool:
                source_id = source_id[:-3]
                target_id = target_id[:-3]
            path_weight_vec = []
            for tail_gene in self.longtail_contig[source_id]:
                for head_gene in self.longhead_contig[target_id]:
                    # print(tail_gene, head_gene)
                    source_id_0 = "C-" + str(tail_gene)
                    target_id_0 = "C-" + str(head_gene)
                    weight_p = 0.0
                    try:
                        p = nx.shortest_path(self.H, source=source_id_0, target=target_id_0)
                    except NetworkXNoPath:
                        p = None
                    if p is not None:
                        for node_p_idx in range(len(p) - 1):
                            weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                        if len(p) >= 2:
                            weight_p = 0.05 + weight_p / float((len(p) - 1) * (len(p) - 1))
                            path_weight_vec.append(weight_p)
            if len(path_weight_vec) > 0:
                return_val = np.mean(path_weight_vec)
            print(return_val)
        return return_val

    def construct_graph(
        self,
        method="graph_alignment",
        sample_id_ref=None,
        min_nucleotides=200,
        min_genes=1,
        edge_weight="unit",
        target_genome_id=-1,
    ):
        """Construct pangenome graph.
        Parameters
        ----------
        method : string
            Method to construct the graph. Valid methods include:
            * graph_alignment
            * graph_free: use gene direction, the input for panta must be .gff (because we need to
            know the gene direction).
        sample_id_ref : integer
            The reference sample, None if none.
        edge_weight: string
            Scheme for edge weights. Valid schemes include:
            * unit: each edge = 1
            * contig_id: weight = contig_ID
            * sample_id: weight = sample_ID
            * adjusted: compute the weight based on the similarity with the target genome
        Returns
        -------
        H : networkx graph
            The pangenome graph
        """
        self.strand = {}
        n_contigs = len(self.gene_position.index)
        self.gene2contigs_dict = {}
        reverse_bool = {}  # 0: no, 1: yes
        contig_dic = {}
        contigName = self.gene_position["ContigName"]
        rows = []
        cols = []
        weight_contig = []
        n_computed_contig = 0
        edge_id_ref = set()
        ref_id = None
        # genome_edges = {} # store list of edges for each genome: Key is genome, values are edges (in map_edge_fn value)
        # for i in range(self.n_samples):
        #     genome_edges[i] = set()
        if sample_id_ref == None:
            # ref_id = [0] #take the first contig as reference
            ref_id = [1]  # take the first contig as reference
            print("Should we use target sequence as reference, hence, the strand will be ok")
            gene_contigs_ref = self.gene_position.iloc[ref_id[0], 2].split(";")
            edge_id_ref = [
                self.map_edge_fn(
                    self.gene2cluster_dict[gene_contigs_ref[i]], self.gene2cluster_dict[gene_contigs_ref[i + 1]]
                )
                for i in range(len(gene_contigs_ref) - 1)
            ]
            edge_id_ref = set(edge_id_ref)
        else:
            ref_id = [i for i in range(n_contigs) if self.gene_position.iloc[i, 0] == sample_id_ref]
            for ref_ in ref_id:
                gene_contigs_ref = self.gene_position.iloc[ref_, 2].split(";")
                edge_id_ref_temp = [
                    self.map_edge_fn(
                        self.gene2cluster_dict[gene_contigs_ref[i]], self.gene2cluster_dict[gene_contigs_ref[i + 1]]
                    )
                    for i in range(len(gene_contigs_ref) - 1)
                ]
                edge_id_ref_temp = set(edge_id_ref_temp)
                edge_id_ref.update(edge_id_ref_temp)

        ### Compute the target mini-graphgenome
        current_genome_id = 0
        target_genome_edge = []
        if edge_weight == "adjusted":
            print("Use target adjusted weight scheme!!!")
            current_genome_id = min(10, n_contigs)
            for i in range(current_genome_id):
                gene_contigs = self.gene_position.iloc[i, 2].split(";")
                edge_id1 = set(
                    [
                        self.map_edge_fn(
                            self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                        )
                        for i in range(len(gene_contigs) - 1)
                    ]
                )
                n1_value = len(edge_id_ref.intersection(edge_id1))
                gene_contigs.reverse()
                edge_id2 = set(
                    [
                        self.map_edge_fn(
                            self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                        )
                        for i in range(len(gene_contigs) - 1)
                    ]
                )
                n2_value = len(edge_id_ref.intersection(edge_id2))
                if n2_value < n1_value:
                    gene_contigs.reverse()
                    edge_id_ref.update(edge_id1)
                    self.strand[self.gene_position.iloc[i, 1]] = "+"
                else:
                    edge_id_ref.update(edge_id2)
                    self.strand[self.gene_position.iloc[i, 1]] = "-"
            target_id = [i for i in range(n_contigs) if self.gene_position.iloc[i, 0] == target_genome_id]
            for target_ in target_id:
                gene_contigs = self.gene_position.iloc[target_, 2].split(";")
                # edge_id1 = set([self.map_edge_fn(self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i+1]]) for i in range(len(gene_contigs)-1)])
                # n1_value = len(edge_id_ref.intersection(edge_id1))
                edge_id1 = [
                    self.map_edge_fn(
                        self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                    )
                    for i in range(len(gene_contigs) - 1)
                ]
                n1_value = len(edge_id_ref.intersection(set(edge_id1)))
                gene_contigs.reverse()
                # edge_id2 = set([self.map_edge_fn(self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i+1]]) for i in range(len(gene_contigs)-1)])
                # n2_value = len(edge_id_ref.intersection(edge_id2))
                edge_id2 = [
                    self.map_edge_fn(
                        self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                    )
                    for i in range(len(gene_contigs) - 1)
                ]
                n2_value = len(edge_id_ref.intersection(set(edge_id2)))
                if n2_value < n1_value:
                    gene_contigs.reverse()
                    target_genome_edge += edge_id1
                else:
                    target_genome_edge += edge_id2
        # string2 = ''
        # for elem in target_genome_edge:
        #     string2 += '_'+ str(elem)
        # print(string2)
        for i in range(current_genome_id, n_contigs):
            gene_contigs = self.gene_position.iloc[i, 2].split(";")
            adjusted_min_genes = (
                min_genes if self.gene_position.iloc[i, 0] == target_genome_id else 10
            )  # use adjusted gene_contigs
            # if self.compute_number_nucleotides(gene_contigs) >= min_nucleotides and len(gene_contigs) >= adjusted_min_genes:
            if len(gene_contigs) >= adjusted_min_genes:
                current_sequence_edge = None
                n_computed_contig = n_computed_contig + 1
                ### align to reference
                if method == "graph_alignment":
                    # if gene_position.iloc[i,0] != sampleID_ref:
                    if i not in ref_id:
                        edge_id1 = set(
                            [
                                self.map_edge_fn(
                                    self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                                )
                                for i in range(len(gene_contigs) - 1)
                            ]
                        )
                        n1_value = len(edge_id_ref.intersection(edge_id1))
                        gene_contigs.reverse()
                        edge_id2 = set(
                            [
                                self.map_edge_fn(
                                    self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i + 1]]
                                )
                                for i in range(len(gene_contigs) - 1)
                            ]
                        )
                        n2_value = len(edge_id_ref.intersection(edge_id2))
                        # if min(n1_value, n2_value) > 3:
                        #     print("ContigID: ", i, ", Contig Length: ", len(gene_contigs),", sample:",self.gene_position.iloc[i,0], ", # of shared edges: ", n1_value, n2_value)
                        if n2_value < n1_value:
                            gene_contigs.reverse()
                            edge_id_ref = edge_id_ref.union(edge_id1)
                            self.strand[self.gene_position.iloc[i, 1]] = "+"
                            current_sequence_edge = edge_id1
                        else:
                            # print("Reverse the sequence: ", i)
                            edge_id_ref = edge_id_ref.union(edge_id2)
                            self.strand[self.gene_position.iloc[i, 1]] = "-"
                            current_sequence_edge = edge_id2

                elif method == "graph_free":
                    ### free alignment
                    ref_id = 0
                    if i == ref_id:
                        reverse_bool[contigName[i]] = 0
                        contig2clusters[contigName[i]] = [gene2cluster_dict[gene] for gene in gene_contigs]
                        contig_dic[contigName[i]] = gene_contigs
                    else:
                        reverse_bool[contigName[i]] = 0
                        this_contig_cluster = [gene2cluster_dict[gene] for gene in gene_contigs]
                        contig2clusters[contigName[i]] = this_contig_cluster
                        contig_dic[contigName[i]] = gene_contigs
                        # check if this cluster intersect with previous cluster
                        for j in range(i):
                            if self.gene_position.iloc[j, 0] != self.gene_position.iloc[i, 0]:
                                # print("pair: ", i, j)
                                list_intersect = set(this_contig_cluster).intersection(
                                    self.contig2clusters[contigName[j]]
                                )
                                if len(list_intersect) > 3:
                                    first_elem = list(list_intersect)[
                                        math.floor(9 * len(list_intersect) / 13)
                                    ]  # math.floor(9*len(list_intersect)/13)
                                    gene1_idx = this_contig_cluster.index(first_elem)
                                    gene2_idx = self.contig2clusters[contigName[j]].index(first_elem)
                                    gene1 = gene_contigs[gene1_idx]
                                    gene2 = contig_dic[contigName[j]][gene2_idx]
                                    # print("contig_id: ", i, j, gene1, gene2, gene2cluster_dict[gene1], gene2cluster_dict[gene2], len(list_intersect))
                                    # if gene1[-1] != gene2[-1]:
                                    if (gene1[-1] == "+" and gene2[-1] == "-") or (
                                        gene1[-1] == "-" and gene2[-1] == "+"
                                    ):
                                        if reverse_bool[contigName[j]] == 0:
                                            ## reverse the sequence (if contigName has not been reversed).
                                            reverse_bool[contigName[i]] = 1
                                            gene_contigs.reverse()
                                            print("===Reverse the sequence", "contig_id: ", i, j, gene1, gene2)
                                            break
                                    break

                else:
                    print("Not implemented yet!")
                ### append the weights
                if self.gene_position.iloc[i, 0] != target_genome_id:
                    for j in range(len(gene_contigs) - 1):
                        rows.append(self.gene2cluster_dict[gene_contigs[j]])
                        cols.append(self.gene2cluster_dict[gene_contigs[j + 1]])
                        if edge_weight == "unit":
                            weight_contig.append(1)
                        elif edge_weight == "adjusted":
                            # factor_num = 15.0
                            # similarity_value = float(float(len(current_sequence_edge.intersection(target_genome_edge)))/len(target_genome_edge))
                            # weight_value = math.exp(factor_num*similarity_value*similarity_value)/math.exp(0.8*0.8*factor_num)
                            if j == 0:
                                print("Compute adjusted scheme")
                                # string1 = ''
                                # for elem in current_sequence_edge:
                                #     string1 += '_'+ str(elem)
                                # weight_value = similarity_sequence(string1, string2)
                                weight_value = max_common_subsequence(list(current_sequence_edge), target_genome_edge)
                                if i <= current_genome_id:
                                    normalize_value = weight_value
                                weight_value = weight_value / normalize_value
                            weight_contig.append(weight_value)
                            # print(round(100*weight_value), end =',')
                        else:
                            print("Not impletement yet!!")
                        # weight the edges also by contig.
                        # weight_contig.append(5*i)
                        ## highlight a contig
                        # highlighted_contig = 1
                        # if i == highlighted_contig:
                        #     weight_contig.append(20)
                        # else:
                        #     weight_contig.append(1)
                        # if self.gene_position.iloc[i,0] in highlight_genome_seq:
                        #     weight_contig.append(5*n_samples)
                        # else:
                        #     weight_contig.append(1)
                        #### weight the edges also by sample.
                        # weight_contig.append(5*self.gene_position.iloc[i,0]+1)

            ### Add contig ID
            for gene in gene_contigs:
                self.gene2contigs_dict[gene] = i
            ### Add head and tail of contig
            contigname_ = self.gene_position.iloc[i, 1]
            self.head_contig[contigname_] = self.gene2cluster_dict[gene_contigs[0]]
            self.tail_contig[contigname_] = self.gene2cluster_dict[gene_contigs[-1]]
            # max_len_contigs = min(5, len(gene_contigs))
            max_len_contigs = min(4, len(gene_contigs))
            self.longhead_contig[contigname_] = [self.gene2cluster_dict[ge] for ge in gene_contigs[0:max_len_contigs]]
            self.longtail_contig[contigname_] = [self.gene2cluster_dict[ge] for ge in gene_contigs[-max_len_contigs:]]

        print(
            "Set minimum on number of nucleotides = ", min_nucleotides, "NUMBER OF COMPUTED CONTIGS:", n_computed_contig
        )
        adj_matrix = csr_matrix((np.array(weight_contig), (rows, cols)), shape=(self.n_clusters, self.n_clusters))
        print("Clip the matrix 0.2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        # adj_matrix = adj_matrix>=10 #version 1: quite good.
        adj_matrix = adj_matrix.multiply(adj_matrix >= 0.1 * self.n_samples)
        self.H = nx.from_numpy_matrix(adj_matrix, create_using=nx.DiGraph)
        ## add node info
        mapping = {i: "C-" + str(i) for i in range(self.n_clusters)}
        self.H = nx.relabel_nodes(self.H, mapping)
        return self.H

    def join_contig(self, sample_id, min_weight=1, method="edge_weight", params=None):
        """Join contigs using pangenome graph.
        Parameters
        ----------
        sample_id : integer
            The sample ID (same as in gene_position) we want to join contigs
        Returns
        -------
        H : networkx graph
            The pangenome graph
        """
        self.sample_df = self.gene_position.loc[self.gene_position["SampleID"] == sample_id]
        source_vec = []
        target_vec = []
        weight_vec = []
        if params["method"] == "weight_path_assembly_v2":
            assembly_graph = params["assembly_graph"]

        for i in range(len(self.sample_df.index)):
            for j in range(len(self.sample_df.index)):
                if params["graph"] == "directed":
                    source_node = self.sample_df.iloc[i, 1] + self.strand[self.sample_df.iloc[i, 1]]
                    target_node = self.sample_df.iloc[j, 1] + self.strand[self.sample_df.iloc[j, 1]]
                else:
                    source_node = self.sample_df.iloc[i, 1]
                    target_node = self.sample_df.iloc[j, 1]
                if i != j:
                    source_id = "C-" + str(self.tail_contig[self.sample_df.iloc[i, 1]])
                    target_id = "C-" + str(self.head_contig[self.sample_df.iloc[j, 1]])
                    if params["method"] == "edge_weight":
                        if self.H.has_edge(source_id, target_id):
                            source_vec.append(self.sample_df.iloc[i, 1])
                            target_vec.append(self.sample_df.iloc[j, 1])
                            weight_vec.append(self.H[source_id][target_id]["weight"])
                    elif params["method"] == "weight_path":
                        try:
                            p = nx.shortest_path(self.H, source=source_id, target=target_id)
                        except NetworkXNoPath:
                            p = None
                        if p is not None:
                            weight_p = 0.0
                            for node_p_idx in range(len(p) - 1):
                                weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                            source_vec.append(self.sample_df.iloc[i, 1])
                            target_vec.append(self.sample_df.iloc[j, 1])
                            weight_p = 1.0 + weight_p / float(len(p) * (len(p) - 1))
                            weight_vec.append(weight_p)
                    elif params["method"] == "weight_path_assembly":
                        assembly_graph = params["assembly_graph"]
                        ### only add the edge in the assembly graph
                        # print(self.sample_df.iloc[i,1], self.sample_df.iloc[j,1])
                        if assembly_graph.has_node(source_node) and assembly_graph.has_node(target_node):
                            try:
                                path = nx.shortest_path(assembly_graph, source=source_node, target=target_node)
                            except NetworkXNoPath:
                                path = None
                            if path is not None and len(path) <= params["max_length"]:
                                # print(self.sample_df.iloc[i,1], self.sample_df.iloc[j,1])
                                try:
                                    p = nx.shortest_path(self.H, source=source_id, target=target_id)
                                except NetworkXNoPath:
                                    p = None
                                if p is not None:
                                    # weight_p = self.H[p[0]][p[1]]['weight']
                                    # for node_p_idx in range(1, len(p)-1):
                                    #     weight_p = min(weight_p, self.H[p[node_p_idx]][p[node_p_idx+1]]['weight'])
                                    # weight_p = weight_p * len(p)
                                    weight_p = 0.0
                                    for node_p_idx in range(len(p) - 1):
                                        weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                                    source_vec.append(self.sample_df.iloc[i, 1])
                                    target_vec.append(self.sample_df.iloc[j, 1])
                                    weight_p = 0.05 + weight_p / float(len(p) * (len(p) - 1))  # 1.0 + ...
                                    weight_vec.append(weight_p)
                                    # print("W:", weight_p, end=",")
                    elif params["method"] == "weight_path_assembly_v2":
                        ### only add the edge in the assembly graph
                        path_len = -1
                        path_nucleotides = -1
                        if assembly_graph.has_node(source_node) and assembly_graph.has_node(target_node):
                            try:
                                path = nx.shortest_path(assembly_graph, source=source_node, target=target_node)
                                for idxp in range(1, len(path) - 1):
                                    path_nucleotides = path_nucleotides + int(path[idxp].split("_")[3])
                                path_len = len(path)
                            except NetworkXNoPath:
                                path = None

                        if path_nucleotides <= params["max_length_nucleotides"]:
                            if False:  ## Only use the end gene in the contig.  print("Old code")
                                try:
                                    p = nx.shortest_path(self.H, source=source_id, target=target_id)
                                except NetworkXNoPath:
                                    p = None
                                if p is not None:
                                    weight_p = 0.0
                                    for node_p_idx in range(len(p) - 1):
                                        weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                                    source_vec.append(self.sample_df.iloc[i, 1])
                                    target_vec.append(self.sample_df.iloc[j, 1])
                                    weight_p = 0.05 + weight_p / float(len(p) * (len(p) - 1))
                                    # weight_p = 1.0 + weight_p/(float(math.sqrt(len(p))*(len(p)-1)))
                                    # weight_p = 1.0 + weight_p/(float((len(p)-1)))
                                    weight_vec.append(weight_p)
                                    # print("W:", weight_p, end=",")
                                else:
                                    if path_len > 0:
                                        source_vec.append(self.sample_df.iloc[i, 1])
                                        target_vec.append(self.sample_df.iloc[j, 1])
                                        weight_vec.append(
                                            float(0.2 * self.n_samples / (path_len * path_len))
                                        )  # gr_output_union_opt
                                        # weight_vec.append(float(0.2*self.n_samples/(path_len))) # gr_output_union_opt
                            else:  ## Use long tail (default)
                                path_weight_vec = []
                                for tail_gene in self.longtail_contig[self.sample_df.iloc[i, 1]]:
                                    for head_gene in self.longhead_contig[self.sample_df.iloc[j, 1]]:
                                        source_id = "C-" + str(tail_gene)
                                        target_id = "C-" + str(head_gene)
                                        weight_p = 0.0
                                        try:
                                            p = nx.shortest_path(self.H, source=source_id, target=target_id)
                                        except NetworkXNoPath:
                                            p = None
                                        if p is not None:
                                            for node_p_idx in range(len(p) - 1):
                                                weight_p += self.H[p[node_p_idx]][p[node_p_idx + 1]]["weight"]
                                            if len(p) >= 2:
                                                weight_p = 0.05 + weight_p / float((len(p) - 1) * (len(p) - 1))
                                                path_weight_vec.append(weight_p)
                                            # print("W:", weight_p, end=",")
                                # if len(path_weight_vec)==0 and  path_len > 0:
                                if len(path_weight_vec) == 0 and path_len > 0:
                                    path_weight_vec.append(float(0.2 * self.n_samples / (path_len * path_len)))
                                    # path_weight_vec.append(float(0.3*self.n_samples/(path_len*path_len*path_len)))
                                    # if i <= 5:
                                    #     print('Assign crazy weight to the Pangraph xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
                                    # path_weight_vec.append(1.054321)
                                if len(path_weight_vec) > 0:
                                    weight_p = np.average(np.array(path_weight_vec))
                                    # print("W:", weight_p, end=",")
                                    if weight_p > 0.005:
                                        source_vec.append(self.sample_df.iloc[i, 1])
                                        target_vec.append(self.sample_df.iloc[j, 1])
                                        ## Reduce the cost by 50% if there is no path in the assembly graph.
                                        if i < 2 and j < 2:
                                            print("Reduce the cost by 1/2 if there is no path in the assembly graph")
                                        if path_nucleotides > 1.0:
                                            weight_vec.append(weight_p)
                                        else:
                                            weight_vec.append(0.5 * weight_p)
                            # else:
                            #     path_weight_vec = []
                            #     for tail_gene in self.longtail_contig[self.sample_df.iloc[i,1]]:
                            #         for head_gene in self.longhead_contig[self.sample_df.iloc[j,1]]:
                            #             source_id = 'C-' + str(tail_gene)
                            #             target_id = 'C-' + str(head_gene)
                            #             weight_p = 0.0
                            #             try:
                            #                 p = nx.shortest_path(self.H, source=source_id, target=target_id)
                            #             except NetworkXNoPath:
                            #                 p = None
                            #             if p is not None:
                            #                 for node_p_idx in range(len(p)-1):
                            #                     weight_p += self.H[p[node_p_idx]][p[node_p_idx+1]]['weight']
                            #                 if len(p) >= 2:
                            #                     weight_p = 0.05 + weight_p/float((len(p)-1)*(len(p)-1))
                            #                     path_weight_vec.append(weight_p)
                            #                     # print("W:", weight_p, end=",")
                            #     # if len(path_weight_vec)==0 and  path_len >= 0:
                            #     if len(path_weight_vec)==0:
                            #         # path_weight_vec.append(float(0.2*self.n_samples/(path_len*path_len)))
                            #         # path_weight_vec.append(1.0)
                            #         if self.weighted_CG.has_edge(source_node, target_node):
                            #             edge_weight = self.weighted_CG[source_node][target_node]['weight']
                            #             print("a", end = "")
                            #             path_weight_vec.append(1.0/edge_weight)
                            #     if len(path_weight_vec) > 0:
                            #         weight_p = np.average(np.array(path_weight_vec))
                            #         # print("W:", weight_p, end=",")
                            #         if weight_p > 0.0005:
                            #             source_vec.append(self.sample_df.iloc[i,1])
                            #             target_vec.append(self.sample_df.iloc[j,1])
                            #             # weight_vec.append(2.0)
                            #             ## Reduce the cost by 50% if there is no path in the assembly graph.
                            #             # if i < 2 and j < 2:
                            #             #     print("Reduce the cost by 1/2 if there is no path in the assembly graph")
                            #             if path_nucleotides > -0.5:
                            #                 # source_vec.append(self.sample_df.iloc[i,1])
                            #                 # target_vec.append(self.sample_df.iloc[j,1])
                            #                 weight_vec.append(weight_p)
                            #                 # weight_vec.append(1.0)
                            #             else:
                            #                 weight_vec.append(0.5*weight_p)
                            #                 # weight_vec.append(1.0*weight_p)
                            #                 # weight_vec.append(1.0)

                    else:
                        print("Not implemented yet!")
        # print("if no connection, use eps weight")
        # contigName2ID = {}
        # n_sample_contig = len(self.sample_df.index)
        # for ido in range(len(self.sample_df.index)):
        #     contigName2ID[self.sample_df.iloc[ido,1]] = ido
        # binary_matrix = np.zeros((n_sample_contig, n_sample_contig))
        # for ido in range(len(source_vec)):
        #     row_id = contigName2ID[source_vec[ido]]
        #     col_id = contigName2ID[target_vec[ido]]
        #     binary_matrix[row_id,col_id] = weight_vec[ido]
        # eps = 0.0001
        # for ido in range(n_sample_contig):
        #     for jdo in range(n_sample_contig):
        #         if binary_matrix[ido, jdo] < eps and ido != jdo:
        #             source_vec.append(self.sample_df.iloc[ido,1])
        #             target_vec.append(self.sample_df.iloc[jdo,1])
        #             weight_vec.append(eps)

        # print("XXXXXXX recompute the node, only use the node with multiplicity = 1")
        # source0= []; target0 = []; weight0 = []
        # for i in range(len(source_vec)):
        #     if self.get_node_coverage(source_vec[i]) <= 22.0 and self.get_node_coverage(target_vec[i]) <= 22.0:
        #         source0.append(source_vec[i])
        #         target0.append(target_vec[i])
        #         weight0.append(weight_vec[i])
        # source_vec = source0; target_vec = target0; weight_vec = weight0

        if self.MLR == 1:
            self.multiplicity_bool = True
        else:
            self.multiplicity_bool = False
        if self.multiplicity_bool:
            print("Implement multiplicity")
            self.compute_multiplicity(self.sample_df)  # compute multiplicity for contigs (nodes)
            source0 = []
            target0 = []
            weight0 = []
            for i in range(len(source_vec)):
                source_node = source_vec[i]
                target_node = target_vec[i]
                for j in range(self.node_multiplicity[source_node]):
                    for k in range(self.node_multiplicity[target_node]):
                        source0.append(source_node + "_m" + str(j))
                        target0.append(target_node + "_m" + str(k))
                        weight0.append(weight_vec[i])
            source_vec = source0
            target_vec = target0
            weight_vec = weight0

        if params["maximum_matching"] == "greedy":
            edge_df = pd.DataFrame({"source": source_vec, "target": target_vec, "weight": weight_vec})
            self.edge_df0 = edge_df
            edge_df = edge_df.sort_values(by="weight", ascending=False)

            if self.MLR == 1:
                nodes_set = set(source_vec + target_vec)
                prev = {}
                after = {}
                for node in nodes_set:
                    prev[node] = None
                    after[node] = None
                self.edge_list = []
                count = 0
                self.long_range_dependancy_threshold = 0.04 * (self.n_samples - 1) - 0.001
                print("Set long_range_dependancy_threshold = ", self.long_range_dependancy_threshold)
                while len(edge_df.index) > 0 and edge_df.iloc[0, 2] >= min_weight:
                    source_sol_temp = edge_df.iloc[0, 0]
                    target_sol_temp = edge_df.iloc[0, 1]
                    if top_prev(prev, source_sol_temp) != target_sol_temp:
                        source_id = recent_large_prev(prev, source_sol_temp, 6000)
                        target_id = recent_large_after(after, target_sol_temp, 6000)
                        # print(source_sol_temp, source_id, '-->', target_sol_temp, target_id)
                        break_path = True
                        if source_id != None and target_id != None:
                            if (
                                self.get_value_long_contigs(self.edge_df0, source_id, target_id)
                                < self.long_range_dependancy_threshold
                            ):
                                break_path = False
                        if break_path:
                            self.edge_list.append([source_sol_temp, target_sol_temp])
                            prev[target_sol_temp] = source_sol_temp
                            after[source_sol_temp] = target_sol_temp
                            edge_df.drop(edge_df[edge_df.source == source_sol_temp].index, inplace=True)
                            edge_df.drop(edge_df[edge_df.target == target_sol_temp].index, inplace=True)
                        else:
                            edge_df = edge_df.iloc[1:, :]  ## Remove the first row
                    else:
                        edge_df = edge_df.iloc[1:, :]  ## Remove the first row
                        print("remove a node which creates a loop")

                    count += 1
                    if count > 20000000:
                        break
            else:
                print("NOT long_range_dependancy")
                self.edge_list = []
                count = 0
                while len(edge_df.index) > 0:
                    count += 1
                    source_sol_temp = edge_df.iloc[0, 0]
                    target_sol_temp = edge_df.iloc[0, 1]
                    # print(source_sol_temp, target_sol_temp, len(edge_df.index))
                    if edge_df.iloc[0, 2] >= min_weight:
                        self.edge_list.append([source_sol_temp, target_sol_temp])
                        edge_df.drop(edge_df[edge_df.source == source_sol_temp].index, inplace=True)
                        edge_df.drop(edge_df[edge_df.target == target_sol_temp].index, inplace=True)
                        if count > 200000000:
                            break
                    else:
                        break

        else:
            print("Compute maximum matching")
            edges = [
                (source_vec[i] + "t", target_vec[i] + "h", weight_vec[i])
                for i in range(len(source_vec))
                if weight_vec[i] >= min_weight
            ]
            mm_graph = nx.Graph()
            mm_graph.add_weighted_edges_from(edges)
            # maximum_matching = nx.max_weight_matching(mm_graph)
            maximum_matching = nx.max_weight_matching(mm_graph, maxcardinality=True)
            self.edge_list = []
            for elem in maximum_matching:
                node1 = elem[0]
                node2 = elem[1]
                if node1[-1] == "t":
                    self.edge_list.append([node1[:-1], node2[:-1]])
                else:
                    self.edge_list.append([node2[:-1], node1[:-1]])

        self.contig_graph = nx.DiGraph()
        self.contig_graph.add_edges_from(self.edge_list)
        return self.contig_graph

    def remove_cycle(self, assembly_graph):
        cycle_list = list(nx.simple_cycles(self.contig_graph))
        if len(cycle_list) > 0:
            print("There is a CYCLE - I will remove it")
        if self.multiplicity_bool:
            for i in range(len(cycle_list)):
                cycle_i = cycle_list[i]
                for j in range(len(cycle_i) - 1):
                    if not assembly_graph.has_edge(cycle_i[j][:-3], cycle_i[j + 1][:-3]):
                        # print("Remove edges:", cycle_i[j], cycle_i[j+1])
                        self.contig_graph.remove_edge(cycle_i[j], cycle_i[j + 1])
                        break
        else:
            for i in range(len(cycle_list)):
                cycle_i = cycle_list[i]
                for j in range(len(cycle_i) - 1):
                    if not assembly_graph.has_edge(cycle_i[j], cycle_i[j + 1]):
                        self.contig_graph.remove_edge(cycle_i[j], cycle_i[j + 1])
                        break
        return self.contig_graph

    def run_pangraph_pipeline(
        self,
        data_dir,
        incomplete_sample_name,
        assem_dir,
        fasta_gen,
        output_dir,
        maximum_matching,
        MLR=1,
        SInfer=0,
        min_weight_val=1.0,
    ):
        contig_dir = assem_dir + "/contigs.fasta"
        ### Read the data
        sample_info = pd.read_csv(data_dir + "/samples.tsv", delimiter="\t", header=None)
        sample_info.columns = ["Name", "SampleID"]
        gene_info = pd.read_csv(data_dir + "/gene_info.tsv", delimiter="\t", header=None)
        gene_info.columns = ["GeneName", "SampleID", "clusterID"]
        gene_position = pd.read_csv(data_dir + "/gene_position.tsv", delimiter="\t", header=None)
        gene_position.columns = ["SampleID", "ContigName", "GeneSequence"]
        # sort by length of contigs
        gene_position.sort_values(by="GeneSequence", key=lambda x: x.str.len(), ascending=False, inplace=True)
        n_samples = len(np.unique(gene_position.iloc[:, 0]))
        incomplete_sample_id = sample_info[sample_info.Name == incomplete_sample_name].iloc[0, 1]
        ### Construct pangraph
        self.__init__(sample_info, gene_info, gene_position)
        # H = self.construct_graph(method = "graph_alignment", sample_id_ref = None,  min_nucleotides = 10, min_genes = 0, edge_weight="adjusted",
        H = self.construct_graph(
            method="graph_alignment",
            sample_id_ref=None,
            min_nucleotides=10,
            min_genes=0,
            edge_weight="unit",
            target_genome_id=incomplete_sample_id,
        )
        edge_list_assembly = []
        self.weighted_CG = None
        if 0:
            print("Use simple assembly graph")
            for l, r in getContigsAdjacency(assem_dir):
                edge_list_assembly.append((append_strand(l), append_strand(r)))
                edge_list_assembly.append((append_strand_reverse(r), append_strand_reverse(l)))
        else:
            print("Use modified  assembly graph")
            linkEdges, self.weighted_CG = getContigsAdjacency_v2(assem_dir)
            for l, r in linkEdges:
                # for l,r in getContigsAdjacency_v2(assem_dir):
                edge_list_assembly.append((append_strand(l), append_strand(r)))
                edge_list_assembly.append((append_strand_reverse(r), append_strand_reverse(l)))

        self.gene = read_contigs2dict(contig_dir)
        # print("Use union graph: overlap + assembly: NO")
        edge_list_overlap = buildOverlapEdge(self.gene, 20, "directed")
        print("Use union graph: overlap + assembly: done")
        edge_list_final = edge_list_overlap + edge_list_assembly
        # edge_list_final = edge_list_overlap
        assembly_graph = nx.DiGraph()
        assembly_graph.add_edges_from(edge_list_final)
        self.assembly_graph = assembly_graph

        if SInfer:
            if 0:
                print("Re-infer the contigs strand: YES")
                self.sample_df = self.gene_position.loc[self.gene_position["SampleID"] == incomplete_sample_id]
                self.contigset_list = list(self.sample_df.iloc[:, 1])
                small_contigs = [cg for cg in self.contigset_list if get_node_length(cg) <= 8000]
                contigstrand_dict = {}
                for ct in small_contigs:
                    contigstrand_dict[ct] = []
                des_dist = 2
                for ctg in self.contigset_list:
                    if get_node_length(ctg) > 15000 and assembly_graph.has_node(ctg + self.strand[ctg]):
                        for node in nx.descendants_at_distance(assembly_graph, ctg + self.strand[ctg], des_dist):
                            if node[:-1] in small_contigs:
                                contigstrand_dict[node[:-1]].append(node[-1])
                        reverse_sign = "+" if self.strand[ctg] == "-" else "-"
                        for node in nx.descendants_at_distance(assembly_graph, ctg + reverse_sign, des_dist):
                            if node[:-1] in small_contigs:
                                contigstrand_dict[node[:-1]].append(node[-1])
                n_infer_contigs = 0
                # self.contigstrand_dict = contigstrand_dict
                for key in contigstrand_dict:
                    if len(contigstrand_dict[key]) > 0:
                        n_infer_contigs += 1
                        self.strand[key] = vote_sign(contigstrand_dict[key])
                print("Number of sign-infered contigs: ", n_infer_contigs)
            else:
                print("Re-infer the contigs strand: YES")
                print("Using gene strand information")
                n_infer_contigs = 0
                self.sample_df = self.gene_position.loc[self.gene_position["SampleID"] == incomplete_sample_id]
                for idx_l in range(len(self.sample_df.index)):
                    contig_name__ = self.sample_df.iloc[idx_l, 1]
                    if get_node_length(contig_name__) <= 15000:
                        n_infer_contigs = n_infer_contigs + 1
                        gene_set_sign = self.sample_df.iloc[idx_l, 1].split(";")
                        count_neg = 0
                        for genes in gene_set_sign:
                            if genes[-2:] == "-1":
                                count_neg = count_neg + 1
                        if 2 * count_neg > len(gene_set_sign):
                            self.strand[contig_name__] = "-"
                        else:
                            self.strand[contig_name__] = "+"

                print("Number of sign-infered contigs: ", n_infer_contigs)

        else:
            print("Re-infer the contigs strand: NO")

        params = {
            "method": "weight_path_assembly_v2",
            "assembly_graph": assembly_graph,
            "max_length": 5,
            "maximum_matching": maximum_matching,
            "graph": "directed",
            "max_length_nucleotides": 8000,
            "weighted_CG": self.weighted_CG,
        }
        self.min_weight = min_weight_val
        self.MLR = MLR
        contig_graph = self.join_contig(sample_id=incomplete_sample_id, min_weight=self.min_weight, params=params)
        # print("DO NOT REMOVE CIRCLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        contig_graph = self.remove_cycle(assembly_graph)
        indegree_dict = dict(contig_graph.in_degree())
        adj_list = {}
        for source_node_key in indegree_dict:
            if indegree_dict[source_node_key] == 0:
                # print(source_node_key, '----------')
                adj_list[source_node_key] = []
                next_neighbor_temp = source_node_key
                while 1:
                    next_neighbor_temp = list(contig_graph.neighbors(next_neighbor_temp))
                    if len(next_neighbor_temp) > 0:
                        next_neighbor_temp = next_neighbor_temp[0]
                        adj_list[source_node_key].append(next_neighbor_temp)
                    else:
                        break

        if self.multiplicity_bool:
            adj_list_v2 = {}
            for key in adj_list:
                key_val = [node[:-3] for node in adj_list[key]]
                adj_list_v2[key] = key_val
            adj_list = adj_list_v2
        ## neu ko la adjacent thi bat dau bang contigs moi.
        self.adj_list0 = adj_list
        adj_list_assembly = {}
        for key in adj_list:
            path1 = adj_list[key].copy()
            if self.multiplicity_bool:
                path1.insert(0, key[:-3])
                path2 = [key[:-3] + self.strand[key[:-3]]]
            else:
                path1.insert(0, key)
                path2 = [key + self.strand[key]]
            for i in range(len(path1) - 1):
                src = path1[i] + self.strand[path1[i]]
                dst = path1[i + 1] + self.strand[path1[i + 1]]
                if not assembly_graph.has_node(src):
                    path2.append(src)
                elif not assembly_graph.has_node(dst):
                    path2.append(src)
                    continue
                else:
                    if nx.has_path(assembly_graph, src, dst):
                        paths = [p for p in nx.all_shortest_paths(assembly_graph, src, dst)]
                        for node in paths[0]:
                            path2.append(node)
                    else:
                        path2.append(src)
                        print("Will test this, Ok?")
                        # # construct a new path if they are disconnected on the graph
                        # if len(path2) > 0:
                        #     adj_list_assembly[new_key+self.strand[new_key]] = self.remove_duplicate(path2)
                        # new_key = dst[0:-1]
                        # path2 = []
            if dst != path2[-1]:
                path2.append(dst)
            if len(path2) > 0:
                adj_list_assembly[key] = self.remove_duplicate(path2)

        ## Remove contigs of more than multiplicity.
        print("Remove contigs of more than multiplicity")
        count_multiplicity = {}
        # Initialize the count_multiplicity
        for key in adj_list_assembly:
            for elem in adj_list_assembly[key]:
                count_multiplicity[elem[:-1]] = 0
        ### Compute base coverage
        print("compute base coverage")
        gene_position_sub = self.sample_df.copy()
        nodes_list = list(gene_position_sub.iloc[:, 1].values)
        nodes_len = [int(node.split("_")[3]) for node in nodes_list]
        nodes_coverage = [float(node.split("_")[5]) for node in nodes_list]
        gene_position_sub["length"] = nodes_len
        gene_position_sub["coverage"] = nodes_coverage
        gene_position_sub = gene_position_sub.sort_values(by="length", ascending=False)
        self.basecoverage = np.median(gene_position_sub["coverage"][:5])
        # recompute the adj_list_assembly
        key_list = [key for key in adj_list_assembly]
        key_list.reverse()
        for key in key_list:
            path1 = []
            path2 = []
            for elem in adj_list_assembly[key]:
                if (elem[:-1] not in path1) and count_multiplicity[elem[:-1]] < self.get_multiplicity(elem):
                    path1.append(elem[:-1])
                    path2.append(elem)
                    count_multiplicity[elem[:-1]] = count_multiplicity[elem[:-1]] + 1
            adj_list_assembly[key] = path2

        # gene_origin = generate_fasta_from_dict(gene, adj_list_assembly, 'all')
        if self.MLR == 1:
            # if 0:
            # print("Refining the scaffolds")
            adj_list_assembly_n = self.remove_unsupported_pangraph(adj_list_assembly)
            adj_list_assembly = adj_list_assembly_n

        self.adj_list_assembly = adj_list_assembly
        gene_origin = generate_fasta_from_dict(self.gene, adj_list_assembly, fasta_gen)
        write_fasta(gene_origin, output_dir)

    def RERUN_pangraph_pipeline(
        self, data_dir, incomplete_sample_name, assem_dir, fasta_gen, output_dir, maximum_matching
    ):
        contig_dir = assem_dir + "/contigs.fasta"
        incomplete_sample_id = self.sample_info[self.sample_info.Name == incomplete_sample_name].iloc[0, 1]
        assembly_graph = self.assembly_graph
        # params = {'method': 'weight_path_assembly_v2', 'assembly_graph': self.assembly_graph, 'max_length': 5, 'maximum_matching': maximum_matching, 'graph':'directed', 'max_length_nucleotides': 8000}
        params = {
            "method": "weight_path_assembly_v2",
            "assembly_graph": assembly_graph,
            "max_length": 5,
            "maximum_matching": maximum_matching,
            "graph": "directed",
            "max_length_nucleotides": 8000,
            "weighted_CG": self.weighted_CG,
        }
        contig_graph = self.join_contig(sample_id=incomplete_sample_id, min_weight=self.min_weight, params=params)
        contig_graph = self.remove_cycle(assembly_graph)
        indegree_dict = dict(contig_graph.in_degree())
        adj_list = {}
        for source_node_key in indegree_dict:
            if indegree_dict[source_node_key] == 0:
                adj_list[source_node_key] = []
                next_neighbor_temp = source_node_key
                while 1:
                    next_neighbor_temp = list(contig_graph.neighbors(next_neighbor_temp))
                    if len(next_neighbor_temp) > 0:
                        next_neighbor_temp = next_neighbor_temp[0]
                        adj_list[source_node_key].append(next_neighbor_temp)
                    else:
                        break

        if self.multiplicity_bool:
            adj_list_v2 = {}
            for key in adj_list:
                key_val = [node[:-3] for node in adj_list[key]]
                adj_list_v2[key] = key_val
            adj_list = adj_list_v2
        ## neu ko la adjacent thi bat dau bang contigs moi.
        self.adj_list0 = adj_list
        adj_list_assembly = {}
        for key in adj_list:
            path1 = adj_list[key].copy()
            if self.multiplicity_bool:
                path1.insert(0, key[:-3])
                path2 = [key[:-3] + self.strand[key[:-3]]]
            else:
                path1.insert(0, key)
                path2 = [key + self.strand[key]]
            for i in range(len(path1) - 1):
                src = path1[i] + self.strand[path1[i]]
                dst = path1[i + 1] + self.strand[path1[i + 1]]
                if not assembly_graph.has_node(src):
                    path2.append(src)
                elif not assembly_graph.has_node(dst):
                    path2.append(src)
                    continue
                else:
                    if nx.has_path(assembly_graph, src, dst):
                        paths = [p for p in nx.all_shortest_paths(assembly_graph, src, dst)]
                        for node in paths[0]:
                            path2.append(node)
                    else:
                        path2.append(src)
                        print("Will test this, Ok?")
                        # # construct a new path if they are disconnected on the graph
                        # if len(path2) > 0:
                        #     adj_list_assembly[new_key+self.strand[new_key]] = self.remove_duplicate(path2)
                        # new_key = dst[0:-1]
                        # path2 = []
            if dst != path2[-1]:
                path2.append(dst)
            if len(path2) > 0:
                adj_list_assembly[key] = self.remove_duplicate(path2)

        self.adj_list_assembly2 = adj_list_assembly
        gene_origin = generate_fasta_from_dict(self.gene, adj_list_assembly, fasta_gen)
        write_fasta(gene_origin, output_dir)
