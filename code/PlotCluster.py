import sys
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Circle

def lattice_pos(L):
    return {i: (i % L, i // L) for i in range(L * L)}

def read_graph(filename):
    import math
    import networkx as nx

    node_data = {}
    edge_data = []
    max_node = 0

    # First pass: gather all node statuses
    with open(filename, 'r') as f:
        # skip first line (comment/header)
        next(f)

        lines = [line.strip() for line in f if line.strip()]
        for line in lines:
            i, j, J, A, B = line.split()
            i = int(i)
            A, B = int(A), int(B)

            if i not in node_data:
                node_data[i] = {'active': False, 'boundary': False}

            if A == 1:
                node_data[i]['active'] = True
            elif B == 1:
                node_data[i]['boundary'] = True

            j = int(j)
            if j not in node_data:
                node_data[j] = {'active': False, 'boundary': False}

            max_node = max(max_node, i, j)

    # Second pass: add full graph + subgraph edges
    G = nx.Graph()
    Gsub = nx.Graph()

    for line in lines:
        i, j, J, A, B = line.split()
        i, j = int(i), int(j)
        J = float(J)

        G.add_edge(i, j, J=J)

        if (node_data[i]['active'] and node_data[j]['active']):
            Gsub.add_edge(i, j, J=J)

    L = int(math.sqrt(max_node + 1))
    assert L * L == max_node + 1, f"Lattice must be square, max node is {max_node}"

    return G, Gsub, node_data, L


def detect_frustrated_loops(G, node_data, type, start_node=None, seed=None):
    if type == 'graph':
        G_active = G
    elif type == 'subgraph':
        G_active = G.subgraph([
            n for n in G.nodes
            if node_data.get(n, {}).get('active', False)
        ])

    if seed is not None:
        np.random.seed(seed)
        import random
        random.seed(seed)

    # Determine search scope
    if start_node is not None:
        if start_node not in G_active.nodes:
            raise ValueError(f"Start node {start_node} not found in active graph")
        for comp_nodes in nx.connected_components(G_active):
            if start_node in comp_nodes:
                search_graph = G_active.subgraph(comp_nodes)
                break
    else:
        search_graph = G_active

    def shortest_cycle_through_edge(graph, u, v):
        """Find shortest cycle using edge (u,v) via BFS from u to v avoiding direct edge."""
        # BFS from u to v in graph without direct (u,v) edge
        from collections import deque
        visited = {u: None}
        queue = deque([u])
        while queue:
            node = queue.popleft()
            for neighbor in graph.neighbors(node):
                if node == u and neighbor == v:
                    continue  # skip the direct edge we're testing
                if neighbor not in visited:
                    visited[neighbor] = node
                    if neighbor == v:
                        # Reconstruct path
                        path = []
                        cur = v
                        while cur is not None:
                            path.append(cur)
                            cur = visited[cur]
                        return list(reversed(path))
                    queue.append(neighbor)
        return None  # no cycle through this edge

    def is_chordless(cycle, graph):
        """Check no shortcut exists between non-adjacent cycle nodes."""
        n = len(cycle)
        cycle_set = set(cycle)
        for i in range(n):
            for j in range(i + 2, n):
                if i == 0 and j == n - 1:
                    continue  # adjacent via wrap-around
                if graph.has_edge(cycle[i], cycle[j]):
                    return False
        return True

    # Collect unique minimal chordless cycles
    seen_cycles = set()
    frustrated_edges = set()
    frustrated_loops = []

    for u, v in search_graph.edges():
        cycle_nodes = shortest_cycle_through_edge(search_graph, u, v)
        if cycle_nodes is None or len(cycle_nodes) < 3:
            continue
        if not is_chordless(cycle_nodes, search_graph):
            continue

        # Canonical key: frozenset of nodes (dedup rotations/reflections)
        key = frozenset(cycle_nodes)
        if key in seen_cycles:
            continue
        seen_cycles.add(key)

        edges = [(cycle_nodes[i], cycle_nodes[(i + 1) % len(cycle_nodes)])
                 for i in range(len(cycle_nodes))]

        # Check all edges exist and compute frustration
        product = 1
        valid = True
        for a, b in edges:
            if not G_active.has_edge(a, b):
                valid = False
                break
            product *= np.sign(G_active[a][b]['J'])

        if not valid:
            continue

        if product < 0:
            frustrated_loops.append(edges)
            for e in edges:
                frustrated_edges.add(frozenset(e))

    return frustrated_edges, frustrated_loops
    
    


#Draw the full graph
def draw_graph(G, L, frustrated_edges, frustrated_loops):
    pos = {i: (i % L, i // L) for i in range(L * L)}
    fig, ax = plt.subplots(figsize=(8, 8))

    # Draw all edges
    for u, v in G.edges():
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        J = G[u][v]['J']
        edge_color = 'tab:green' if J > 0 else 'tab:red'

        dx = abs(x1 - x2)
        dy = abs(y1 - y2)

        if dx > 1:  # horizontal wrap
            if x1 < x2:
                ax.plot([x1, x1 - 1], [y1, y1], color=edge_color, lw=2.5, zorder=1)
                ax.plot([x2, x2 + 1], [y2, y2], color=edge_color, lw=2.5, zorder=1)
            else:
                ax.plot([x1, x1 + 1], [y1, y1], color=edge_color, lw=2.5, zorder=1)
                ax.plot([x2, x2 - 1], [y2, y2], color=edge_color, lw=2.5, zorder=1)
        elif dy > 1:  # vertical wrap
            if y1 < y2:
                ax.plot([x1, x1], [y1, y1 - 1], color=edge_color, lw=2.5, zorder=1)
                ax.plot([x2, x2], [y2, y2 + 1], color=edge_color, lw=2.5, zorder=1)
            else:
                ax.plot([x1, x1], [y1, y1 + 1], color=edge_color, lw=2.5, zorder=1)
                ax.plot([x2, x2], [y2, y2 - 1], color=edge_color, lw=2.5, zorder=1)
        else:
            ax.plot([x1, x2], [y1, y2], color=edge_color, lw=2.5, zorder=1)


    # Draw all nodes
    for n, (x, y) in pos.items():
        ax.add_patch(Circle((x, y), 0.35, facecolor='white', edgecolor='tab:gray', zorder=2))
        ax.text(x, y, str(n), color='black', ha='center', va='center', fontsize=4, zorder=3)


    #Draw plaquettes
    for x in range(L):
        for y in range(L):
            n0 = x + y * L
            n1 = ((x + 1) % L) + y * L
            n2 = ((x + 1) % L) + ((y + 1) % L) * L
            n3 = x + ((y + 1) % L) * L
            edges = [frozenset((n0, n1)), frozenset((n1, n2)),
                     frozenset((n2, n3)), frozenset((n3, n0))]

            # All 4 edges must exist in Gsub AND be in frustrated_edges
            if all(
                (lambda u, v: G.has_edge(u, v) or G.has_edge(v, u))(*tuple(e)) and
                e in frustrated_edges
                for e in edges
            ):
                xc, yc = (x + 0.5) % L, (y + 0.5) % L
                ax.plot([xc - 0.1, xc + 0.1], [yc - 0.1, yc + 0.1], 'r', lw=2.5)
                ax.plot([xc - 0.1, xc + 0.1], [yc + 0.1, yc - 0.1], 'r', lw=2.5)


    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='tab:green', lw=2.5, label='J>0 bond'),
        Line2D([0], [0], color='tab:red', lw=2.5, label='J<0 bond'),
        Line2D([0], [0], marker='x', color='tab:red', label='Frustrated loop',
               linestyle='None', markersize=10, markeredgewidth=2),
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=3,
        frameon=True
    )

    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()


#Draw the subgraph
def draw_subgraph(G, Gsub, node_data, L, frustrated_edges, frustrated_loops):
    pos = {i: (i % L, i // L) for i in range(L * L)}
    fig, ax = plt.subplots(figsize=(8, 8))

    #Draw edges
    for u, v in G.edges():
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        J = G[u][v]['J']
        u_data = node_data.get(u, {})
        v_data = node_data.get(v, {})

        if not u_data.get('active') and not v_data.get('active'):
            edge_color = 'tab:gray'
            linestyle = '-'
            alpha = 0.3
            zorder = 0
        else:
            edge_color = 'tab:green' if J > 0 else 'tab:red'
            linestyle = '--' if ( u_data.get('boundary') or v_data.get('boundary') ) else '-'
            alpha = 1.0
            zorder = 100

        dx = abs(x1 - x2)
        dy = abs(y1 - y2)

        if dx > 1:  # horizontal wrap
            if x1 < x2:
                ax.plot([x1, x1 - 1], [y1, y1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
                ax.plot([x2, x2 + 1], [y2, y2], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
            else:
                ax.plot([x1, x1 + 1], [y1, y1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
                ax.plot([x2, x2 - 1], [y2, y2], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
        elif dy > 1:  # vertical wrap
            if y1 < y2:
                ax.plot([x1, x1], [y1, y1 - 1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
                ax.plot([x2, x2], [y2, y2 + 1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
            else:
                ax.plot([x1, x1], [y1, y1 + 1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
                ax.plot([x2, x2], [y2, y2 - 1], color=edge_color, lw=2.5, zorder=zorder, linestyle=linestyle, alpha=alpha)
        else:
            ax.plot([x1, x2], [y1, y2], color=edge_color, lw=2.5, zorder=1, linestyle=linestyle, alpha=alpha)


    # Draw nodes
    for n, (x, y) in pos.items():
        d = node_data.get(n, {'active': False, 'boundary': False})
        
        if d['boundary']:
            face = 'tab:gray'
            edge = 'black'
        elif d['active']:
            face = 'tab:blue'
            edge = 'black'
        else:
            face = 'white'
            edge = 'tab:gray'
        ax.add_patch(plt.Circle((x, y), 0.35, facecolor=face, edgecolor=edge, zorder=2))
        ax.text(x, y, str(n), color='black', ha='center', va='center', fontsize=4, zorder=3)


    #Draw plaquettes
    for x in range(L):
        for y in range(L):
            n0 = x + y * L
            n1 = ((x + 1) % L) + y * L
            n2 = ((x + 1) % L) + ((y + 1) % L) * L
            n3 = x + ((y + 1) % L) * L
            edges = [frozenset((n0, n1)), frozenset((n1, n2)),
                     frozenset((n2, n3)), frozenset((n3, n0))]

            # All 4 edges must exist in Gsub AND be in frustrated_edges
            if all(
                (lambda u, v: Gsub.has_edge(u, v) or Gsub.has_edge(v, u))(*tuple(e)) and
                e in frustrated_edges
                for e in edges
            ):
                xc, yc = (x + 0.5) % L, (y + 0.5) % L
                ax.plot([xc - 0.1, xc + 0.1], [yc - 0.1, yc + 0.1], 'r', lw=2.5)
                ax.plot([xc - 0.1, xc + 0.1], [yc + 0.1, yc - 0.1], 'r', lw=2.5)


    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='black', label='Internal node',
               markerfacecolor='tab:blue', markersize=10),
        Line2D([0], [0], marker='o', color='black', label='Boundary node',
               markerfacecolor='tab:gray', markersize=10),
        Line2D([0], [0], marker='o', color='tab:gray', label='Inactive node',
               markerfacecolor='white', markersize=10),
        Line2D([0], [0], color='tab:green', lw=2.5, label='J>0 bond'),
        Line2D([0], [0], color='tab:red', lw=2.5, label='J<0 bond'),
        Line2D([0], [0], marker='x', color='tab:red', label='Frustrated loop',
               linestyle='None', markersize=10, markeredgewidth=2),
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, 1.05),
        ncol=3,
        frameon=True
    )

    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()


def print_average_connectivity(G):
    degrees = [deg for _, deg in G.degree()]
    avg_conn = sum(degrees) / len(degrees) if degrees else 0
    print(f"\nAverage connectivity ⟨C⟩: {avg_conn:.2f}")




#-------------------------------------------------------------------------------
#------------------------------------ MAIN -------------------------------------
#-------------------------------------------------------------------------------

if len(sys.argv) != 1:
    infile = sys.argv[1]
 
else:
    sys.stderr.write('\n SYNTAX ERROR!: program.py <graph.conf> \n\n')
    sys.exit(1)

G, Gsub, node_data, L = read_graph(infile)

#Graph info
print_average_connectivity(G)

frustrated_edges_all, frustrated_loops_all = detect_frustrated_loops(G, node_data, 'graph', start_node=None, seed=None)
print(f"\nNumber of frustrated loops (full graph): {len(frustrated_loops_all)}\n")

S = sum(1 for v in node_data.values() if ( v['active'] and not v['boundary']) )
B = sum(1 for v in node_data.values() if v['boundary'])
print(f"Subgraph size (S): {S}")
print(f"Boundary size (B): {B}")
frustrated_edges, frustrated_loops = detect_frustrated_loops(Gsub, node_data, 'subgraph', start_node=None, seed=None)
print(f"Number of frustrated loops (subgraph): {len(frustrated_loops)}\n")
#print(frustrated_loops)

draw_graph(G, L, frustrated_edges_all, frustrated_loops_all)
draw_subgraph(G, Gsub, node_data, L, frustrated_edges, frustrated_loops)
plt.show()













