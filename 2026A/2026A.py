# # Computational Mathematics
# ## Project 2026A
# ## Candidate Number: 1105769

# General references and documentation used:
# - https://docs.python.org/3.12/library/stdtypes.html
# - https://numpy.org/doc/stable/reference/index.html

from itertools import combinations
import numpy as np
from pprint import pprint
import zipfile

# Helper functions to load the complexes:

# Print out files inside the .zip
with zipfile.ZipFile("complexes.zip") as z:
    filenames = z.namelist()
print(f"Available filenames: {filenames}")
print()

def load_complex_npz(file):
    data = np.load(file)
    cells = data["cells"]
    coordinates = data["coordinates"]
    return (cells, coordinates)

# Load a complex from within the .zip
def load_complex_zip(zipname, complexname):
    with zipfile.ZipFile(zipname) as z:
        with z.open(complexname) as f:
            return load_complex_npz(f)


# ***
# ### Question A.1.
#
# > Write a function `build_complex(cells)` that builds a data
# > structure complex that is a list of dictionaries. Each dictionary
# > `complex[dim]` should map a tuple of vertices of length `dim` with
# > vertices listed in ascending order to an id number. The id numbers
# > should be natural numbers beginning from zero. The vertices should
# > be plain Python integers (not e.g. `np.int32`). For example, for the
# > complex of Figure A.4 that is stored in `triangle.npz`, the output
# > should be
# >
# > ```
# > # complex[dim] = {sorted_vertex_tuple: simplex_id}
# > complex = [
# >     {(0,): 0, (1,): 1, (2,): 2}, # dimension 0 (vertices)
# >     {(0, 1): 0, (1, 2): 1, (0, 2): 2}, # dimension 1 (edges)
# >     {(0, 1, 2): 0} # dimension 2 (triangles)
# > ]
# > ```
# >
# > Only simplices actually present in the complex should be listed in the
# > dictionary: for example, if vertices `v3` and `v8` are not connected by an
# > edge, `(3, 8)` must not be a key in `complex[1]`. The map encoded in
# > `complex[0]` should always be the identity.
# > Apply your code to the complexes in `triangle.npz` and
# > `triangle_minus_triangle.npz`. Run the following diagnostic
# > code on the two complexes:
# >
# > ```
# > from pprint import pprint
# > for (dim, c) in enumerate(complex):
# >     print(f"Printing complex of dimension {dim}")
# >     pprint(c)
# > ```
# >
# > [5 marks]

def pprint_complex(complex):
    for (dim, c) in enumerate(complex):
        print(f"Printing complex of dimension {dim}")
        pprint(c)

# build_complex() constructs an ordered basis for a simplicial complex K.
def build_complex(cells):
    # > The first records the connections between the highest-dimensional
    # > cells and vertices, with each row corresponding to one cell.
    # Therefore each cell will be the same dimension
    dim = cells.shape[1]

    # Using sets so I don't have to deal with deduplicating the sub-simplices
    simplices = [set() for _ in range(dim)]
    for cell in cells.tolist():
        # Ensure each simplex is represented uniquely with vertices in ascending order
        cell = sorted(cell)
        for k in range(dim):
            # Adds all possible sub-simplices of dimension k from the current cell
            simplices[k].update(combinations(cell, k + 1))

    # Assign each simplex to a unique id in lexicographic order
    complex = [{simplex: id for id, simplex in enumerate(sorted(k_simplices))}
               for k_simplices in simplices]
    return complex

# Some additional comments about my solution:
# > The map encoded in `complex[0]` should always be the identity
#
# My solution above would encode an identity map iff the vertices in the cells are
# always labelled with `0, 1, ..., n`, which I think is very reasonable to assume.

# Even though the sample output for the 1-simplices (edges) provided is
# > {(0, 1): 0, (1, 2): 1, (0, 2): 2}
#
# but there is no criteria provided for deciding if $(i_1, j_1)$ has a smaller
# id than $(i_2, j_2)$, apart from that $i_1 < j_1$ and $i_2 < j_2$.
# Therefore I simply choose lexicographic ordering to assign id.
# so I would instead have the edge ids be `{(0, 1): 0, (0, 2): 1, (1, 2): 2}`

def question_a1():
    print("Question A.1")

    for filename in ["triangle.npz", "triangle_minus_triangle.npz"]:
        print()
        print(f"The cells of {filename} are:")
        cells, _ = load_complex_zip("complexes.zip", filename)
        print(cells)
        print()

        complex = build_complex(cells)
        pprint_complex(complex)

question_a1()


# ***
# ### Question A.2.
#
# > Write a function `boundary_matrix(complex, k)`
# > that builds the matrix representation of the boundary operator
# > $\delta_k : C_k \to C_{k - 1}$ as a numpy array.
# >
# > Print all boundary operators for the complexes stored in
# > `triangle.npz` and `triangle_minus_triangle.npz`.
# >
# > [5 marks]

def boundary_matrix(complex, k):
    # Handle the case k = 0 seperately.
    # We shall represent delta_0 as a matrix with zero rows, since the codomain
    # is C_{-1} = \varnothing.
    if k == 0:
        return np.zeros((0, len(complex[k])), dtype=np.int32)

    # delta_k : C_k \to C_{k-1}
    # Domain of delta_k has a basis consisting of complex[k]
    # Codomain of delta_k has a basis consisting of complex[k - 1]
    op = np.zeros((len(complex[k - 1]), len(complex[k])), dtype=np.int32)
    for simplex, simplex_id in complex[k].items():
        for i in range(k + 1):
            # Removing i-th index to get boundary sub-simplex
            boundary = simplex[:i] + simplex[i + 1:]
            # Find corresponding boundary id
            boundary_id = complex[k - 1][boundary]
            # Update the matrix representation of the operator.
            op[boundary_id][simplex_id] = (-1)**i

    return op

def question_a2():
    print("Question A.2")

    for filename in ["triangle.npz", "triangle_minus_triangle.npz"]:
        print()
        print("=" * 10)
        print(filename)
        print("=" * 10)
        cells, _ = load_complex_zip("complexes.zip", filename)
        complex = build_complex(cells)
        pprint_complex(complex)
        dim = len(complex)
        print(f"The boundary operators of {filename} are:")
        for k in range(dim):
            print()
            print(f"delta_{k} = ")
            print(boundary_matrix(complex, k))
    print()

question_a2()


# ***
# ### Question A.3.
#
# > For each complex in `complexes.zip`, build all
# > boundary operators, and verify that (A.1.11) holds.
# >
# > $$\delta_k \circ \delta_{k + 1} = 0. \tag{A.1.11}$$
# >
# > [2 marks]

def question_a3():
    print("Question A.3")

    for filename in filenames:
        cells, _ = load_complex_zip("complexes.zip", filename)
        complex = build_complex(cells)
        dim = len(complex)

        # Compute all boundary operators once beforehand
        boundary_ops = [boundary_matrix(complex, k) for k in range(dim)]

        print()
        for k in range(dim - 1):
            # Compute the representation of delta_k \circ delta_{k + 1}
            boundary_of_boundary = boundary_ops[k] @ boundary_ops[k + 1]
            # Check if it is the zero map
            is_zero = (boundary_of_boundary == 0).all()
            print(f"{filename} delta_{k} circ delta_{k + 1} is zero? {is_zero}")
            # Used assert because if it is for whatever reason False, something went wrong for sure.
            assert is_zero, f"\t{filename} delta_{k} circ delta_{k + 1} is not zero!"
    print()

question_a3()


# ***
# ### Question A.4.
#
# > Using the rank-nullity theorem, derive a simpler
# > expression for $b_k$ in terms of matrix ranks from (A.2.8).
# >
# > $$b_k = \dim Z_k - \dim B_k = \dim \ker(\delta_k) - \dim \operatorname{im}(\delta_{k + 1}) \tag{A.2.8}$$
# >
# > [3 marks]

# Since $\delta_k : C_k \to C_{k - 1}$, rank-nullity theorem tells us that
# $$\dim C_k = \dim \ker(\delta_k) + \dim \operatorname{im}(\delta_k).$$
# Therefore
# \begin{align*}
#     b_k &= \dim \ker(\delta_k) - \dim \operatorname{im}(\delta_{k + 1}) \\
#         &= \dim C_k - \dim \operatorname{im}(\delta_k) - \dim \operatorname{im}(\delta_{k + 1}),
# \end{align*}
# where $\dim \operatorname{im}(\delta_k)$ and $\dim \operatorname{im}(\delta_{k + 1})$
# are the ranks of the matrix representation of the $\delta_k$ and $\delta_{k + 1}$ operator respectively.


# ***
# ### Question A.5.
#
# > Using the expression you derived in Question A.4,
# > implement a function `betti_numbers(cells)` that computes the
# > Betti numbers for a given `cells` array.
# >
# > Compute and print the Betti numbers for each complex stored in
# > `complexes.zip`.
# >
# > [5 marks]

def betti_numbers(cells):
    complex = build_complex(cells)

    dim = len(complex)
    rank_delta_k = np.zeros((dim + 1,), dtype=np.int32)

    # Note that delta_0 has rank 0 since delta_0: C_0 \to C_{-1} = \varnothing.
    # Hence we can start calculating from k = 1 to dim.
    for k in range(1, dim):
        delta_k = boundary_matrix(complex, k)
        # We only need the rank of the boundary operators to compute the Betti numbers
        rank_delta_k[k] = np.linalg.matrix_rank(delta_k)

    bk = [0] * dim

    for k in range(dim):
        # Using the formula we derived in A.4.
        # Also cast type from np.int32 to int.
        bk[k] = int(len(complex[k]) - rank_delta_k[k] - rank_delta_k[k + 1])

    return bk

def question_a5():
    print("Question A.5")
    print()
    for filename in filenames:
        cells, _ = load_complex_zip("complexes.zip", filename)
        print(f"The Betti numbers for {filename} are", betti_numbers(cells))
    print()

question_a5()

# ***
# Some additional test cases, where we try to compute the Betti numbers of
# triangulations of certain manifolds.
# 
# Since the computation of the Betti numbers relies on our answer to A.1, A.2 and A.4,
# checking that the Betti numbers computed for these manifolds are correct is 
# evidence that these questions have been correctly implemented.

# Expected Betti numbers for
# - Sphere $S^2$: $b_0 = 1, b_1 = 0, b_2 = 1$.
# - Torus $S^1 \times S^1$: $b_0 = 1, b_1 = 2, b_2 = 1$.
# - Klein bottle: $b_0 = 1, b_1 = 1, b_2 = 0$.

# Visualisations of the triangulations used are available in `./triangulation-diagrams/`

def additional_tests():
    test_cases = {
        "sphere": np.array([
            [0, 4, 1],
            [1, 5, 2],
            [2, 6, 3],
            [1, 4, 5],
            [2, 5, 6],
            [3, 6, 2],
            [4, 7, 5],
            [5, 8, 6],
            [6, 9, 2],
            [5, 7, 8],
            [6, 8, 9],
            [2, 9, 1],
            [7, 10, 8],
            [8, 7, 9],
            [9, 4, 1],
            [8, 10, 7],
            [9, 7, 4]
        ], dtype=np.int32),
        "torus" : np.array([
            [0, 3, 1],
            [1, 4, 2],
            [2, 5, 0],
            [1, 3, 4],
            [2, 4, 5],
            [0, 5, 3],
            [3, 6, 4],
            [4, 7, 5],
            [5, 8, 3],
            [4, 6, 7],
            [5, 7, 8],
            [3, 8, 6],
            [6, 0, 7],
            [7, 1, 8],
            [8, 2, 6],
            [7, 0, 1],
            [8, 1, 2],
            [6, 2, 0]
        ], dtype=np.int32),
        "klein_bottle": np.array([
            [0, 3, 1],
            [1, 4, 2],
            [2, 5, 0],
            [1, 3, 4],
            [2, 4, 5],
            [0, 5, 3],
            [3, 6, 4],
            [4, 7, 5],
            [5, 8, 3],
            [4, 6, 7],
            [5, 7, 8],
            [3, 8, 6],
            [6, 9, 10],
            [7, 10, 11],
            [8, 11, 9],
            [6, 10, 7],
            [7, 11, 8],
            [8, 9, 6],
            [9, 2, 10],
            [10, 1, 11],
            [11, 0, 9],
            [10, 2, 1],
            [11, 1, 0],
            [9, 0, 2]
        ], dtype=np.int32)
    }

    ground_truth_betti_numbers = {
        "sphere": [1, 0, 1],
        "torus": [1, 2, 1],
        "klein_bottle": [1, 1, 0]
    }

    for complex_name, cells in test_cases.items():
        found_betti_numbers = betti_numbers(cells)
        assert found_betti_numbers == ground_truth_betti_numbers[complex_name], f"Error: {complex_name} should have Betti numbers {ground_truth_betti_numbers}, but computed {found_betti_numbers} instead!"
        print(f"Correctly found that {complex_name} has Betti numbers {found_betti_numbers}")

additional_tests()

