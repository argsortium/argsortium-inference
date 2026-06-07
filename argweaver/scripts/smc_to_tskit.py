#!/usr/bin/env python3
"""
Convert an ARGweaver .smc(.gz) file to a tskit TreeSequence (SINGER-style).

This conversion is similar to how SINGER converts to tree sequences:
  * node identity PERSISTS across recombination breakpoints (one node id,
    edges spanning multiple local trees)
  * recombination is left IMPLICIT (a child simply switches parents at the
    breakpoint position) -- there are NO recombination-event nodes,
  * the tree sequence therefore stores recombination POSITIONS but not
    recombination TIMES.

If recombination times are needed, read them back from the .smc SPR lines
(each SPR records pos / recomb_node / recomb_time / coal_node / coal_time).

Technical details:
ArgWeaver produces tree sequences differing by 1 SPR. Tree sequence and SPR is defined in the .smc.gz file.
Internal node names are recycled, for example if Node 10 recombines above Node 17. The parent of Node 10 connects directly to the child of Node 10. Now the Node 10 label gets recycled.
This script adds a new internal node every SPR event and keeps track of its identity using the parent_child function.
Since Argweaver has discrete time, many coalescent nodes can have the same time. In cases where a parent has same age as child, the parent age is nudged bythe next float point value
Since recombination nodes are not stored, sometimes can end up with cycles for e.g.
[0,x) A is a parent of B
[x,L) B is a parent of A
In these rare cases (<2% of internal nodes) Node A is split into A1 and A2.
"""

import gzip
import bisect
import collections
import numpy as np
import dendropy
import tskit


# ---------------------------------------------------------------------------
# Reading the SMC file
# ---------------------------------------------------------------------------
def read_smc(filename):
    """Read an SMC file and return (names, region, tree_records, spr).

    names        : list of leaf names
    region       : [start, end] (1-indexed, as written in the file)
    tree_records : list of (start, end, newick_string)
    spr          : list of dicts with pos / recomb_node / recomb_time /
                   coal_node / coal_time (node labels refer to the tree
                   immediately to the SPR's left)
    """
    names = []
    region = []
    tree_records = []
    spr = []
    open_fn = gzip.open if filename.endswith(".gz") else open
    with open_fn(filename, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            tag = fields[0]
            if tag == "NAMES":
                names = fields[1:]
            elif tag == "REGION":
                region = [int(fields[2]), int(fields[3])]
            elif tag == "TREE":
                tree_records.append((int(fields[1]), int(fields[2]), fields[3]))
            elif tag == "SPR":
                spr.append({
                    "pos": int(fields[1]),
                    "recomb_node": int(fields[2]),
                    "recomb_time": float(fields[3]),
                    "coal_node": int(fields[4]),
                    "coal_time": float(fields[5]),
                })
    return names, region, tree_records, spr


def read_sites(filename):
    """Read an ARGweaver .sites file and return (names, region, sites).

    names  : list of leaf names (must match the SMC NAMES line)
    region : [start, end] (1-indexed, as written in the file)
    sites  : dict {position (int, 1-based) : allele_string} for VARIANT sites
             only (a site is variant if it has more than one distinct allele).
             The allele string has one character per leaf, in `names` order.
    """
    names = []
    region = []
    sites = {}
    open_fn = gzip.open if filename.endswith(".gz") else open
    with open_fn(filename, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            tag = fields[0]
            if tag == "NAMES":
                names = fields[1:]
            elif tag == "REGION":
                region = [int(fields[2]), int(fields[3])]
            else:
                alleles = fields[1]
                if len(set(alleles)) > 1:          # keep variant sites only
                    sites[int(fields[0])] = alleles
    return names, region, sites


# ---------------------------------------------------------------------------
# Parsing a single local tree
# ---------------------------------------------------------------------------
def get_dendropy_node_id(node):
    """Return the SMC integer label of a dendropy node."""
    if node._label is None:
        return int(node.taxon._label)
    return int(node._label)


def parse_topology(nwk, taxa):
    """Parse a newick string into (parent, time, is_leaf, root).

    All dicts are keyed by the SMC integer label.
    parent[child_label] = parent_label   (root has no entry)
    """
    tree = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=taxa,
                             extract_comment_metadata=True, rooting="force-rooted")
    tree.calc_node_ages()
    parent = {}
    time = {}
    is_leaf = {}
    root = None
    for node in tree.preorder_node_iter():
        label = get_dendropy_node_id(node)
        time[label] = node.age
        is_leaf[label] = len(node.child_nodes()) == 0
        if node.parent_node is None:
            root = label
        else:
            parent[label] = get_dendropy_node_id(node.parent_node)
    return parent, time, is_leaf, root


def edge_pairs(parent, mapping):
    """Return the set of (child_tskit_id, parent_tskit_id) for a tree.

    `parent` is the SMC parent dict; `mapping` translates SMC labels to
    persistent tskit node ids.
    """
    pairs = set()
    for child_label, parent_label in parent.items():
        pairs.add((mapping[child_label], mapping[parent_label]))
    return pairs


# ---------------------------------------------------------------------------
# Building the tskit tables
# ---------------------------------------------------------------------------
def add_initial_nodes(tables, tree0):
    """Add the nodes of tree 0; leaves first so samples get ids 0..n-1.

    Returns (cur, leaf_labels):
      cur         : SMC label -> persistent tskit id
      leaf_labels : sorted list of leaf SMC labels
    """
    parent, time, is_leaf, root = tree0
    leaf_labels = sorted(label for label in is_leaf if is_leaf[label])
    internal_labels = sorted(label for label in is_leaf if not is_leaf[label])
    cur = {}
    for label in leaf_labels:
        cur[label] = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=time[label])
    for label in internal_labels:
        cur[label] = tables.nodes.add_row(flags=0, time=time[label])
    return cur, leaf_labels


def topological_order(tables):
    """Return node ids ordered children-before-parents (Kahn's algorithm)."""
    children = collections.defaultdict(set)
    parents = collections.defaultdict(set)
    for edge in tables.edges:
        children[edge.parent].add(edge.child)
        parents[edge.child].add(edge.parent)
    indegree = {}
    for node_id in range(tables.nodes.num_rows):
        indegree[node_id] = len(children[node_id])
    queue = collections.deque()
    for node_id in indegree:
        if indegree[node_id] == 0:
            queue.append(node_id)
    order = []
    while queue:
        node_id = queue.popleft()
        order.append(node_id)
        for parent_id in parents[node_id]:
            indegree[parent_id] -= 1
            if indegree[parent_id] == 0:
                queue.append(parent_id)
    return order, children


def find_cycle_nodes(tables):
    """Return the set of node ids that lie in a directed cycle.

    ARGweaver's discretized times let several coalescences share one time;
    across loci their ancestry can flip, producing a directed cycle in the
    union of all marginal-tree edges. Such a cycle has no valid global time
    ordering, so tskit cannot store it. We find them as strongly-connected
    components (size > 1) of the union child->parent graph.
    """
    children_of = collections.defaultdict(list)   # parent -> children (child->parent reversed)
    for edge in tables.edges:
        children_of[edge.parent].append(edge.child)

    # iterative Tarjan strongly-connected components
    index_of = {}
    lowlink = {}
    on_stack = set()
    stack = []
    counter = [0]
    cycle_nodes = set()

    def strongconnect(start):
        work = [(start, 0)]
        while work:
            node, child_i = work[-1]
            if child_i == 0:
                index_of[node] = counter[0]
                lowlink[node] = counter[0]
                counter[0] += 1
                stack.append(node)
                on_stack.add(node)
            recursed = False
            kids = children_of[node]
            while child_i < len(kids):
                child = kids[child_i]
                if child not in index_of:
                    work[-1] = (node, child_i + 1)
                    work.append((child, 0))
                    recursed = True
                    break
                if child in on_stack:
                    lowlink[node] = min(lowlink[node], index_of[child])
                child_i += 1
            if recursed:
                continue
            if lowlink[node] == index_of[node]:
                component = []
                while True:
                    member = stack.pop()
                    on_stack.discard(member)
                    component.append(member)
                    if member == node:
                        break
                if len(component) > 1:
                    cycle_nodes.update(component)
            work.pop()
            if work:
                parent_node = work[-1][0]
                lowlink[parent_node] = min(lowlink[parent_node], lowlink[node])

    for node in range(tables.nodes.num_rows):
        if node not in index_of:
            strongconnect(node)
    return cycle_nodes


def break_tied_cycles(tables):
    """Split cycle nodes by their edge-interval boundaries to remove cycles.

    A cycle node is duplicated into one copy per genomic sub-interval over
    which its set of edges is constant. Within any single position only one
    copy is present, so every marginal tree is unchanged; across positions the
    conflicting ancestral roles are carried by different ids, which makes the
    union graph acyclic. Copies keep the original node's time and flags.
    """
    cycle_nodes = find_cycle_nodes(tables)
    if not cycle_nodes:
        return  # already acyclic (e.g. cycle-free inputs)

    edges = list(tables.edges)

    # boundaries[node] = sorted coordinates where node's edge set changes
    boundaries = collections.defaultdict(set)
    for edge in edges:
        if edge.parent in cycle_nodes:
            boundaries[edge.parent].add(edge.left)
            boundaries[edge.parent].add(edge.right)
        if edge.child in cycle_nodes:
            boundaries[edge.child].add(edge.left)
            boundaries[edge.child].add(edge.right)
    for node in cycle_nodes:
        boundaries[node] = sorted(boundaries[node])

    # one id per sub-interval; first sub-interval reuses the original id
    times = tables.nodes.time
    flags = tables.nodes.flags
    piece_id = {}  # (node, left_coordinate_of_piece) -> tskit id
    for node in cycle_nodes:
        coords = boundaries[node]
        for k in range(len(coords) - 1):
            left = coords[k]
            if k == 0:
                piece_id[(node, left)] = node
            else:
                piece_id[(node, left)] = tables.nodes.add_row(
                    flags=int(flags[node]), time=float(times[node]))

    def id_at(node, pos):
        coords = boundaries[node]
        k = bisect.bisect_right(coords, pos) - 1
        return piece_id[(node, coords[k])]

    # rebuild edges, clipping at the boundaries of any cycle endpoint
    new_edges = []
    for edge in edges:
        parent_split = edge.parent in cycle_nodes
        child_split = edge.child in cycle_nodes
        if not parent_split and not child_split:
            new_edges.append((edge.left, edge.right, edge.parent, edge.child))
            continue
        cuts = {edge.left, edge.right}
        if parent_split:
            cuts.update(b for b in boundaries[edge.parent] if edge.left < b < edge.right)
        if child_split:
            cuts.update(b for b in boundaries[edge.child] if edge.left < b < edge.right)
        cuts = sorted(cuts)
        for left, right in zip(cuts, cuts[1:]):
            parent_id = id_at(edge.parent, left) if parent_split else edge.parent
            child_id = id_at(edge.child, left) if child_split else edge.child
            new_edges.append((left, right, parent_id, child_id))

    tables.edges.clear()
    for left, right, parent_id, child_id in new_edges:
        tables.edges.add_row(left, right, parent_id, child_id)


def break_zero_length_branches(tables):
    """Nudge each parent strictly above its children.

    ARGweaver's discretized times produce zero-length branches
    (parent.time == child.time), which tskit rejects. np.nextafter gives the
    smallest representable increase at any magnitude. Processing children
    before parents lets the fix cascade up chains of tied nodes.
    """
    order, children = topological_order(tables)
    times = np.array(tables.nodes.time, dtype=np.float64)
    for node_id in order:
        for child_id in children[node_id]:
            if times[node_id] <= times[child_id]:
                times[node_id] = np.nextafter(times[child_id], np.inf)
    tables.nodes.time = times


def build_tree_sequence(tree_records, spr, names, taxa):
    """SMC trees + SPRs -> tskit TreeSequence with persistent node identity.

    Identity rule (argweaver's get_local_node_mapping): across a breakpoint
    every label keeps its node id EXCEPT the broken node (parent of the
    recomb branch), whose recycled label denotes a brand-new coalescent node.
    """
    trees = [parse_topology(nwk, taxa) for _, _, nwk in tree_records]
    n_trees = len(trees)

    # ARGweaver's REGION end is 1-based inclusive; use end+1 so a variant sitting
    # exactly at the inclusive end coordinate falls inside the half-open [0, L)
    # sequence and is spanned by the last tree (relate does the same).
    region_end = tree_records[-1][1]
    tables = tskit.TableCollection(sequence_length=region_end + 1)
    cur, leaf_labels = add_initial_nodes(tables, trees[0])
    assert len(leaf_labels) <= len(names)

    # open an edge for every parent-child pair in tree 0
    active = {}  # child_id -> (parent_id, left_coordinate)
    first_left = tree_records[0][0] - 1
    for child_id, parent_id in edge_pairs(trees[0][0], cur):
        active[child_id] = (parent_id, first_left)

    # walk the breakpoints: each SPR turns tree i into tree i+1
    for i in range(n_trees - 1):
        pos = spr[i]["pos"]
        broken = trees[i][0][spr[i]["recomb_node"]]   # parent of recomb branch
        new_cur = dict(cur)                           # copy current mapping
        new_cur[broken] = tables.nodes.add_row(flags=0, time=spr[i]["coal_time"])
        desired = edge_pairs(trees[i + 1][0], new_cur)

        # close edges that no longer exist in the next tree
        for child_id in list(active):
            parent_id, left = active[child_id]
            if (child_id, parent_id) not in desired:
                tables.edges.add_row(left, pos, parent_id, child_id)
                del active[child_id]
        # open edges that newly appear in the next tree
        for child_id, parent_id in desired:
            if child_id not in active:
                active[child_id] = (parent_id, pos)

        cur = new_cur

    # close everything still open at the end of the sequence (extend the last
    # tree to region_end + 1 so it spans a variant at the inclusive end coord)
    last_right = region_end + 1
    for child_id, (parent_id, left) in active.items():
        tables.edges.add_row(left, last_right, parent_id, child_id)

    break_tied_cycles(tables)            # split tied-time cycles (discretized ages)
    break_zero_length_branches(tables)   # strict parent>child via nextafter
    tables.sort()
    return tables.tree_sequence(), leaf_labels


# ---------------------------------------------------------------------------
# Verification: each marginal tree must match the SMC tree
# ---------------------------------------------------------------------------
def smc_clades(nwk, taxa):
    """Set of clades (frozensets of leaf labels) in an SMC newick tree."""
    tree = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=taxa,
                             rooting="force-rooted")
    clades = set()
    for node in tree.preorder_internal_node_iter():
        clade = frozenset(int(leaf.taxon.label) for leaf in node.leaf_iter())
        clades.add(clade)
    return clades


def ts_clades(tree, leaf_labels, label_to_tskit):
    """Set of clades for one tskit marginal tree, in SMC-leaf-label terms."""
    clades = set()
    for node_id in tree.nodes():
        if not tree.is_internal(node_id):
            continue
        clade = set()
        for leaf_label in leaf_labels:
            if tree.is_descendant(label_to_tskit[leaf_label], node_id):
                clade.add(leaf_label)
        clades.add(frozenset(clade))
    return clades


def verify_marginal_trees(ts, tree_records, leaf_labels, taxa):
    """Return True if every marginal tree matches its SMC tree's clades."""
    label_to_tskit = {}
    for index, leaf_label in enumerate(leaf_labels):
        label_to_tskit[leaf_label] = index   # leaves were added first (ids 0..n-1)
    for start, end, nwk in tree_records:
        midpoint = (start - 1 + end) // 2
        if ts_clades(ts.at(midpoint), leaf_labels, label_to_tskit) != smc_clades(nwk, taxa):
            return False
    return True


# ---------------------------------------------------------------------------
# Mutations: map each variant site onto a branch of its local tree
# ---------------------------------------------------------------------------
def is_monophyletic(tree, carriers):
    """Return True if the carrier samples form a clade in `tree`.

    carriers: iterable of sample node ids carrying the allele of interest.
    """
    carriers = set(carriers)
    if len(carriers) <= 1:
        return True                            # 0 or 1 sample is trivially a clade
    mrca = tree.mrca(*carriers)                 # tskit accepts multiple node ids
    if mrca == tskit.NULL:
        return False                           # carriers span more than one root
    return set(tree.samples(mrca)) == carriers  # exactly these ids, no extras


def _branch_age(times, node, tree):
    """Midpoint age of the branch above `node` (the allele age).

    tskit needs mutation.time strictly below the parent node's time. On
    epsilon-perturbed (1-ulp) branches the midpoint rounds up to the parent's
    time, so we fall back to the node's own time there.
    """
    parent = tree.parent(node)
    if parent == tskit.NULL:
        return float(times[node])
    lo = float(times[node])
    hi = float(times[parent])
    mid = lo + (hi - lo) / 2.0
    return mid if lo < mid < hi else lo


def map_mutations_to_branches(ts, sites, names, site_names):
    """Map variant sites onto branches of the tree sequence.

    `names`      : SMC leaf names in tskit-sample-id order (sample i -> names[i]).
    `site_names` : leaf names in the sites file's allele-string order. The two
                   lists hold the SAME leaves but may be in a DIFFERENT order, so
                   alleles are looked up by name, not by position.

    For each variant site (binary ancestral-vs-derived, like arg-summarize):
      * if the derived-allele carriers form a clade -> ONE mutation on the
        branch above their MRCA (ancestral = major allele);
      * else if the ancestral-allele carriers form a clade -> polarization is
        flipped: ONE mutation on the branch above the non-carrier MRCA, with the
        ancestral state set to the (minor) derived allele;
      * else (neither monophyletic, or multiallelic) -> tskit's Hartigan
        parsimony `tree.map_mutations`, which may place several mutations.
    Mutation time is the branch midpoint (allele age).

    Returns (ts_with_mutations, n_monophyletic, n_flipped, n_hartigan), where
    n_monophyletic counts carrier-clade placements (default polarization),
    n_flipped counts non-carrier-clade placements (ancestral state flipped),
    and n_hartigan counts sites resolved by Hartigan parsimony.
    """
    tables = ts.dump_tables()
    times = ts.tables.nodes.time
    n_samples = ts.num_samples
    n_mono = 0
    n_flipped = 0
    n_hartigan = 0

    # allele-string position for each tskit sample id (map leaves by name)
    sites_index = {name: j for j, name in enumerate(site_names)}
    perm = [sites_index[names[i]] for i in range(n_samples)]

    for pos in sorted(sites):
        alleles = sites[pos]
        site_pos = pos                          # keep the VCF/.sites bp position
                                                # (matches relate, threads, and the simulation truth;
                                                #  a prior `pos - 1` shifted every site 1 bp left)
        if site_pos < 0 or site_pos >= ts.sequence_length:
            continue
        tree = ts.at(site_pos)
        if tree.num_edges == 0:
            continue                            # position not spanned by a tree

        counts = collections.Counter(alleles)
        distinct = sorted(counts, key=lambda a: counts[a], reverse=True)
        ancestral = distinct[0]                 # most common = ancestral (default)
        derived = distinct[1] if len(distinct) > 1 else ancestral
        carriers = [i for i in range(n_samples) if alleles[perm[i]] != ancestral]
        non_carriers = [i for i in range(n_samples) if alleles[perm[i]] == ancestral]

        placed = False
        if len(distinct) == 2:
            if is_monophyletic(tree, carriers):
                node = carriers[0] if len(carriers) == 1 else tree.mrca(*carriers)
                site_id = tables.sites.add_row(position=site_pos, ancestral_state=ancestral)
                tables.mutations.add_row(site=site_id, node=node,
                                         derived_state=derived,
                                         time=_branch_age(times, node, tree))
                n_mono += 1
                placed = True
            elif is_monophyletic(tree, non_carriers):
                # flip polarization: the major allele arose once on the non-carrier clade
                node = non_carriers[0] if len(non_carriers) == 1 else tree.mrca(*non_carriers)
                site_id = tables.sites.add_row(position=site_pos, ancestral_state=derived)
                tables.mutations.add_row(site=site_id, node=node,
                                         derived_state=ancestral,
                                         time=_branch_age(times, node, tree))
                n_flipped += 1
                placed = True

        if not placed:
            # Hartigan parsimony (handles non-monophyletic and multiallelic sites)
            allele_list = distinct
            index = {allele: k for k, allele in enumerate(allele_list)}
            genotypes = np.fromiter((index[alleles[perm[i]]] for i in range(n_samples)),
                                    dtype=np.int8, count=n_samples)
            anc_state, mutations = tree.map_mutations(genotypes, alleles=tuple(allele_list))
            site_id = tables.sites.add_row(position=site_pos, ancestral_state=anc_state)
            for mutation in mutations:
                tables.mutations.add_row(site=site_id, node=mutation.node,
                                         derived_state=mutation.derived_state,
                                         time=_branch_age(times, mutation.node, tree))
            n_hartigan += 1

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence(), n_mono, n_flipped, n_hartigan


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def reorder_samples(ts, current_sample_names, desired_names):
    """Permute the sample nodes so that sample k corresponds to desired_names[k].

    ARGweaver's `arg-sample` rewrites the NAMES line of its .smc output in a
    SCRAMBLED order (not the input VCF/.sites order). build_tree_sequence keeps
    that SMC order, so tskit sample i ends up holding haplotype `current_sample_names[i]`
    (the SMC NAMES order) rather than the canonical VCF haplotype i. Left
    unfixed, every per-sample / per-pair metric compares the wrong haplotypes.

    We remap the samples to `desired_names` (the .sites/VCF order). Internal
    nodes are kept and simply reindexed by tskit; topology, edges and mutations
    are unchanged.
    """
    assert set(current_sample_names) == set(desired_names), \
        "sample name sets differ; cannot reorder"
    name_to_node = {nm: int(node) for node, nm in zip(ts.samples(), current_sample_names)}
    new_sample_order = [name_to_node[nm] for nm in desired_names]
    keep = set(new_sample_order)
    non_samples = [n for n in range(ts.num_nodes) if n not in keep]
    return ts.subset(new_sample_order + non_samples)


def convert(filename, sites_file=None, add_mutations=True, verify=False):
    """Convert an SMC file to a tskit TreeSequence.

    If `add_mutations` is True and a `sites_file` is given, variant sites are
    mapped onto branches after the node/edge tables are built. The sites file's
    leaf names must exactly match the SMC NAMES line.

    The returned tree sequence's samples are ordered to match the `sites_file`
    (i.e. the input VCF haplotype order), not ARGweaver's scrambled SMC NAMES
    order. See `reorder_samples`.
    """
    names, region, tree_records, spr = read_smc(filename)
    print(f"names: {len(names)} | region: {region} | "
          f"trees: {len(tree_records)} | sprs: {len(spr)}")
    taxa = dendropy.TaxonNamespace()
    ts, leaf_labels = build_tree_sequence(tree_records, spr, names, taxa)
    print(f"built TS: {ts.num_nodes} nodes, {ts.num_edges} edges, {ts.num_trees} trees")
    if verify:
        ok = verify_marginal_trees(ts, tree_records, leaf_labels, taxa)
        print("ALL MARGINAL TREES MATCH" if ok else "MISMATCH FOUND")

    if add_mutations and sites_file is not None:
        site_names, site_region, sites = read_sites(sites_file)
        # same leaves must be present (order may differ; alleles are mapped by name)
        assert set(site_names) == set(names) and len(site_names) == len(names), (
            "leaf identity differs between sites file and SMC file")
        ts, n_mono, n_flipped, n_hartigan = map_mutations_to_branches(
            ts, sites, names, site_names)
        print(f"mutations: {ts.num_sites} variant sites, {ts.num_mutations} mutations "
              f"| monophyletic={n_mono} flipped={n_flipped} hartigan={n_hartigan}")

        # tskit sample s holds SMC haplotype names[leaf_labels[s]]; reorder to
        # the .sites/VCF haplotype order so sample k == site_names[k].
        current_sample_names = [names[leaf_labels[s]] for s in range(ts.num_samples)]
        ts = reorder_samples(ts, current_sample_names, site_names)
        print(f"reordered samples to .sites/VCF haplotype order ({ts.num_samples} samples)")
    return ts


def main():
    import sys
    if len(sys.argv) < 2:
        print("usage: python smc_to_tskit.py <input.smc[.gz]> [output.trees] [sites_file]")
        sys.exit(1)
    sites_file = sys.argv[3] if len(sys.argv) > 3 else None
    ts = convert(sys.argv[1], sites_file=sites_file)
    if len(sys.argv) > 2:
        ts.dump(sys.argv[2])
        print(f"written {sys.argv[2]}")


if __name__ == "__main__":
    main()
