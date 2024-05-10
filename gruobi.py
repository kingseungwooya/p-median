from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

def generate_grid_points(grid_shape):
    grid_size = grid_shape[0] * grid_shape[1]
    grid = list(product(list(range(grid_shape[0])), repeat=2))
    print(grid)
    grid = np.array(grid)
    return grid, grid_size

def generate_distance_matrix(ST_points, grid_points):
    distance_matrix = {(i, j): np.hypot(ST_points[i][0] - grid_points[j][0], ST_points[i][1] - grid_points[j][1])
                       for i in range(len(ST_points)) for j in range(len(grid_points))}
    # print(distance_matrix)
    return distance_matrix

def draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode):
    plt.figure(figsize=(12, 5))
    plt.scatter(ST_points[:, 0], ST_points[:, 1], color='blue')

    for i in range(len(ST_points)):
        plt.annotate(f'$d_{i}={flow[i]}$', (ST_points[i, 0] - 0.5, ST_points[i, 1] - 5))

    for n in supernodesInVertex:
        plt.scatter(grid_points[n, 0], grid_points[n, 1], color='green', marker='D', s=150)

    for n in node_to_supernode:
        i, j = n
        plt.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]])

    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.title("Nodes")
    plt.grid()
    plt.show()

def solve_P_Median(ST_points, grid_points, flow, P):
    ST_index = range(len(ST_points))
    vertex_index = range(len(grid_points))
    position = [(i, j) for i in ST_index for j in vertex_index]
    print(position)
    distance_matrix = generate_distance_matrix(ST_points, grid_points)

    model = Model('P-Median')

    z = model.addVars(position, vtype=GRB.BINARY, name='X')
    y = model.addVars(vertex_index, vtype=GRB.BINARY, name='Y')

    model.setObjective(
        quicksum(distance_matrix[i, j] * z[i, j] for i, j in position)
        , GRB.MINIMIZE
    )

    # s.t 10
    model.addConstrs(
        quicksum(z[i, j] for j in vertex_index) == 1 for i in ST_index
    )
    # s.t 11
    model.addConstrs(
        z[i, j] - y[j] <= 0 for i in ST_index for j in vertex_index
    )
    # # s.t 12
    # model.addConstrs(
    #
    # )
    # # s.t 13
    # model.addConstrs(
    #
    # )
    # s.t 14
    model.addConstr(
        quicksum(y[j] for j in vertex_index) <= P
    )

    model.optimize()

    supernodesInVertex = [k for k in vertex_index if y[k].x == 1]
    node_to_supernode = [k for k in position if z[k].x == 1]
    print(supernodesInVertex)
    print(node_to_supernode)
    return supernodesInVertex, node_to_supernode



def read_coords(path):
    coords = []
    with open(path, "r") as f:
        for line in f.readlines():
            line = [float(x.replace("\n", "")) for x in line.split(" ")]
            coords.append(line)
    return np.array(coords)

def main():

    N = 10
    # ST_points = np.random.randint(low=0, high=10, size=(N, 2))

    ST_points = read_coords("coord2.txt")

    # flow
    flow = list(np.random.randint(low=1, high=5, size=N))


    grid_shape = (10, 10)
    grid_points, grid_size = generate_grid_points(grid_shape)

    # P value
    P = 2

    supernodesInVertex, node_to_supernode = solve_P_Median(ST_points, grid_points, flow, P)


    draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode)

if __name__ == "__main__":
    main()
