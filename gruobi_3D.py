from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from dotenv import load_dotenv
import os

# load .env
load_dotenv()

GUROBI_LICENSE_KEY = os.environ.get('GRUOBI_LICENSE')
Env(GUROBI_LICENSE_KEY)

def generate_grid_points(grid_shape):
    grid_size = grid_shape[0] * grid_shape[1] * grid_shape[2]
    grid = list(product(list(range(grid_shape[0])), repeat=3))
    grid = np.array(grid)
    # print(grid)
    return grid, grid_size

def generate_distance_matrix(ST_points, grid_points):
    distance_matrix = {}
    for i in range(len(ST_points)):
        for j in range(len(grid_points)):
            distance = np.sqrt((ST_points[i][0] - grid_points[j][0])**2 +
                               (ST_points[i][1] - grid_points[j][1])**2 +
                               (ST_points[i][2] - grid_points[j][2])**2)
            distance_matrix[(i, j)] = distance
    return distance_matrix


def draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # ST points 그리기
    ax.scatter(ST_points[:, 0], ST_points[:, 1], ST_points[:, 2], color='blue')

    # Supernodes 그리기
    for n in supernodesInVertex:
        ax.scatter(grid_points[n, 0], grid_points[n, 1], grid_points[n, 2], color='green', marker='D', s=150)

    # Node_to_supernode 연결선 그리기
    for n in node_to_supernode:
        i, j = n
        ax.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]],
                [ST_points[i, 2], grid_points[j, 2]], color='black')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Nodes')
    ax.grid()

    plt.show()

def solve_P_Median(ST_points, grid_points, flow, P):
    ST_index = range(len(ST_points))
    vertex_index = range(len(grid_points))
    position = [(i, j) for i in ST_index for j in vertex_index]

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

    ST_points = read_coords("coord3D.txt")
    print(ST_points)

    # flow
    flow = list(np.random.randint(low=1, high=5, size=N))


    grid_shape = (10, 10, 10)
    grid_points, grid_size = generate_grid_points(grid_shape)

    # P value
    P = 2

    supernodesInVertex, node_to_supernode = solve_P_Median(ST_points, grid_points, flow, P)


    draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode)

if __name__ == "__main__":
    main()
