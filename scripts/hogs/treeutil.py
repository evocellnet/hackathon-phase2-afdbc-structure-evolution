"""Minimal rooted-tree utilities (newick I/O, pruning, MRCA, clades) for the panel tree and HOG selection."""
from __future__ import annotations

import re


class Node:
    __slots__ = ("name", "length", "children", "parent", "age", "id")

    def __init__(self, name="", length=None):
        self.name = name; self.length = length; self.children = []; self.parent = None; self.age = None
        self.id = None

    def add(self, c):
        c.parent = self; self.children.append(c); return c

    def is_leaf(self):
        return not self.children

    def preorder(self):
        st = [self]
        while st:
            n = st.pop(); yield n
            st.extend(reversed(n.children))

    def postorder(self):
        st = [(self, False)]
        while st:
            n, done = st.pop()
            if done:
                yield n
            else:
                st.append((n, True)); st.extend((c, False) for c in reversed(n.children))

    def leaves(self):
        return [n for n in self.preorder() if n.is_leaf()]


_TOK = re.compile(r"\s*('(?:[^']|'')*'|[(),:;]|[^(),:;']+)")


def parse_newick(s: str) -> Node:
    toks = _TOK.findall(s.strip())
    root = Node(); cur = root; i = 0; expect_len = False; stack = []
    cur_node = None
    while i < len(toks):
        t = toks[i].strip()
        if t == "(":
            n = Node(); cur.add(n)
            stack.append(cur); cur = n; cur_node = None
        elif t == ",":
            cur_node = None
        elif t == ")":
            cur_node = cur; cur = stack.pop()
        elif t == ":":
            expect_len = True
        elif t == ";":
            break
        elif t:
            if expect_len:
                (cur_node if cur_node is not None else cur).length = float(t); expect_len = False
            else:
                name = t[1:-1].replace("''", "'") if t.startswith("'") else t
                if cur_node is not None:   # label after ')'
                    cur_node.name = name
                else:
                    cur_node = cur.add(Node(name))
        i += 1
    # the outer dummy root holds one child = the real root
    r = root.children[0] if len(root.children) == 1 and not root.name else root
    r.parent = None
    return r


def _q(name: str) -> str:
    if name and re.search(r"[\s(),:;'\[\]]", name):
        return "'" + name.replace("'", "''") + "'"
    return name or ""


def to_newick(root: Node, lengths=True, internal_names=True) -> str:
    out = {}
    for n in root.postorder():
        s = ""
        if n.children:
            s = "(" + ",".join(out.pop(id(c)) for c in n.children) + ")"
            if internal_names:
                s += _q(n.name)
        else:
            s = _q(n.name)
        if lengths and n.length is not None and n.parent is not None:
            s += f":{n.length:.6g}"
        out[id(n)] = s
    return out[id(root)] + ";"


def prune_to(root: Node, keep: set) -> Node:
    """Induced subtree on leaves named in ``keep`` (unary nodes suppressed, lengths summed)."""
    keepn = {}
    for n in root.postorder():
        if n.is_leaf():
            if n.name in keep:
                c = Node(n.name, n.length); keepn[id(n)] = c
        else:
            kids = [keepn.pop(id(c)) for c in n.children if id(c) in keepn]
            if not kids:
                continue
            if len(kids) == 1:
                k = kids[0]
                if k.length is not None and n.length is not None:
                    k.length += n.length
                keepn[id(n)] = k
            else:
                m = Node(n.name, n.length)
                for k in kids:
                    m.add(k)
                keepn[id(n)] = m
    r = keepn[id(root)]
    r.parent = None
    return r


def index_tree(root: Node):
    """Assign preorder ids; return (nodes list, leafset dict id->frozenset of leaf names)."""
    nodes = list(root.preorder())
    for i, n in enumerate(nodes):
        n.id = i
    leafset = {}
    for n in root.postorder():
        leafset[n.id] = frozenset([n.name]) if n.is_leaf() else frozenset().union(*(leafset[c.id] for c in n.children))
    return nodes, leafset


def mrca(nodes, leafset, names) -> Node:
    """Deepest node whose leaf set contains all ``names`` (preorder scan from the leaves up)."""
    s = set(names)
    best = nodes[0]
    for n in nodes:
        if s <= leafset[n.id] and len(leafset[n.id]) < len(leafset[best.id]):
            best = n
    return best
