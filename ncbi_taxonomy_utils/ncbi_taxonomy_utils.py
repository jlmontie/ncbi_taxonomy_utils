import time
import sys
import os
import pickle

class ncbi_taxonomy:

    def __init__(self):
        """
        :param resource_dir: directory containing the ncbi taxonomy resources produced by ncbi_taxonomy_utils_constructor.py
        :param canonical_taxa: set of 'canonical' ranks to use with get_lineage_names_ranks
        """
        package_dir = os.path.dirname(__file__)
        with open(os.path.join(package_dir, 'ncbi_taxonomy_merged.dat'), 'rb') as merged_file:
            self.merged = pickle.load(merged_file)  # {old taxid}->new taxid
        with open(os.path.join(package_dir, 'ncbi_taxonomy_nodes_rel.dat'), 'rb') as nodes_rel_file:
            self.nodes_rel = pickle.load(nodes_rel_file)  # {taxid}->parent
        with open(os.path.join(package_dir, 'ncbi_taxonomy_children_nodes.dat'), 'rb') as children_nodes_file:
            self.children_nodes = pickle.load(children_nodes_file)  # {taxid}->[children]
        with open(os.path.join(package_dir, 'ncbi_taxonomy_nodes_rank.dat'), 'rb') as nodes_rank_file:
            self.nodes_rank = pickle.load(nodes_rank_file)  # {taxid}->rank
        with open(os.path.join(package_dir, 'ncbi_taxonomy_names.dat'), 'rb') as names_file:
            self.names = pickle.load(names_file)  # {taxid}->scientific name
        with open(os.path.join(package_dir, 'ncbi_taxonomy_all_names.dat'), 'rb') as all_names_file:
            self.all_names = pickle.load(all_names_file)  # {taxid}->[(name, type)]
        with open(os.path.join(package_dir, 'ncbi_taxonomy_max_taxid.dat'), 'rb') as max_taxid_file:
            self.max_taxid = pickle.load(max_taxid_file)
        with open(os.path.join(package_dir, 'ncbi_taxonomy_canonical_taxa.dat'), 'rb') as canonical_taxa_file:
            self.canonical_taxa = pickle.load(canonical_taxa_file)

    def get_lca(self, taxid1, taxid2):
        if taxid1 == taxid2:
            return taxid1
        path1 = set(self.get_path(taxid1))
        for tx in self.get_path(taxid2):
            if tx in path1:
                return tx
        return 0

    def get_lca_from_list(self, taxids):
        lca = taxids[0]
        for tx in taxids:
            lca = self.get_lca(lca, tx)
        return lca

    def get_max_taxid(self):
        return self.max_taxid

    def get_name(self, taxid):
        return self.names.get(self.merged.get(taxid, taxid), None)

    def get_rank(self, taxid):
        return self.nodes_rank.get(self.merged.get(taxid, taxid), None)

    def get_path(self, taxid):
        """
        :param taxid:  taxid of which to get path to root
        :return: list of taxids from taxid to root
        """
        taxid = self.merged.get(taxid, taxid)
        path = []
        while taxid > 0:
            path.append(taxid)
            taxid = self.nodes_rel.get(taxid, 0)
        return path

    def get_children(self, taxid):
        """
        :param taxid:  taxid for which to get all descendant nodes
        :return: list of descendants
        """
        to_visit = set([x for x in self.children_nodes.get(taxid, [])])
        visited_nodes = set()
        while len(to_visit - visited_nodes) > 0:
            work = to_visit - visited_nodes
            for tx in work:
                visited_nodes.add(tx)
                for child in self.children_nodes.get(tx, []):
                    to_visit.add(child)
        return list(to_visit)

    def get_all_names(self, taxid):
        """
        :param taxid:  taxid for which to get all descendant nodes
        :return: list of descendants
        """
        taxid = self.merged.get(taxid, taxid)
        return self.all_names.get(taxid, [])

    def get_lineage_tx_names_string(self, taxid):
        """
        :param taxid:  taxid from which to compute lineage information
        :return:       string that has the <taxid>:<name>; from leaf to root
        """
        string_info = []
        for tx in self.get_path(taxid):
            string_info.append("%d:%s" % (tx, self.get_name(tx)))
        return ";".join(string_info)

    def get_lineage_lists(self, taxid):
        """
        :param taxid:   taxid from which to compute lineage information
        :return:        two lists: [taxids] [names] [ranks]
        """
        taxids = []
        names = []
        ranks = []
        for tx in self.get_path(taxid):
            taxids.append(tx)
            names.append(self.get_name(tx))
            ranks.append(self.get_rank(tx))
        return taxids, names, ranks

    def get_lineage_names_ranks(self, taxid, canonical=False):
        """
        :param taxid: get lineage name and rank information for taxid
        :param canonical: only use ranks specified by canonical_taxa in constructor
        :return: dictionary {rank}->[name, taxid].  If canonical is True, then include all canonical ranks even if None
        """
        taxid = self.merged.get(taxid, taxid)
        lineage_names_ranks = {}
        if canonical is True:
            for k in self.canonical_taxa:
                lineage_names_ranks[k] = ["NONE", -1]
        while taxid > 0:
            rank = self.nodes_rank.get(taxid, "NONE")
            if canonical is True:
                if rank in self.canonical_taxa:
                    lineage_names_ranks[rank] = [self.names.get(taxid, "NONE"), taxid]
            else:
                lineage_names_ranks[rank] = [self.names.get(taxid, "NONE"), taxid]
            taxid = self.nodes_rel.get(taxid, 0)
        return lineage_names_ranks

    def create_subset_tri(self, input_taxids_file, output_file):
        """
        :param input_taxids_file:  input file with 1 taxid per line
        :param output_file:  output tri file including paths necessary for input_taxids_file
        :return: Nothing
        """
        use_taxids = set()
        for line in open(input_taxids_file, 'r'):
            taxid = int(line.strip())
            path = self.get_path(taxid)
            for tx in path:
                use_taxids.add(tx)

        use_taxids.add(1)
        out = open(output_file, 'w')
        for tx in use_taxids:
            parent = self.nodes_rel.get(tx, 0)
            if parent == 0 and tx != 1:
                parent = 1
                sys.stderr.write("WARNING: taxid %d has no parent.  Parent assigned 1\n" % (tx))
            out.write("%d\t%d\n" % (tx, parent))
        out.close()

    def get_species_taxid(self, taxid):
        names_rank = self.get_lineage_names_ranks(taxid, canonical=True)
        if names_rank['species'][0] == "NONE":
            return None
        else:
            return names_rank['species'][1]

    def get_species_taxid_if_exists(self, taxid):
        names_rank = self.get_lineage_names_ranks(taxid, canonical=True)
        if names_rank['species'][0] == "NONE":
            return taxid
        else:
            return names_rank['species'][1]
