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
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # 처음 10개의 점을 파란색으로 그리기
    ax.scatter(ST_points[:10, 0], ST_points[:10, 1], ST_points[:10, 2], color='blue')

    # 나머지 점을 빨간색으로 그리기
    ax.scatter(ST_points[10:, 0], ST_points[10:, 1], ST_points[10:, 2], color='red')

    # Supernodes 그리기
    for n in supernodesInVertex:
        # ax.scatter(grid_points[n, 0], grid_points[n, 1], grid_points[n, 2], color='green', marker='o', s=150)
        ax.scatter(grid_points[n, 0], grid_points[n, 1], 9, color='green', marker='o', s=150)

    # Node_to_supernode 연결선 그리기
    for n in node_to_supernode:
        i, j = n
        # ax.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]],
        #         [ST_points[i, 2], grid_points[j, 2]], color='black')
        ax.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]],
                [ST_points[i, 2], 9], color='black')

    # Super Node 끼리 연결선 그리기
    print("supernode index들")
    print(supernodesInVertex)
    for i in range(len(supernodesInVertex)):

        try:
            index1 = i
            index2 = -i
            src_sn_point = grid_points[supernodesInVertex[i]]
            sink_sn_point = grid_points[supernodesInVertex[i+2]]
            print(src_sn_point[0])

            ax.plot([src_sn_point[0], sink_sn_point[0]],[src_sn_point[1], sink_sn_point[1]],
                    [9, 9], color='green', lw=5)
        except IndexError:
            continue

    # ax.view_init(0, 60)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    # 축 범위 설정
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_zlim(0, 10)
    ax.set_title('Tray Way Structure')
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

    N = 20
    # ST_points = np.random.randint(low=0, high=10, size=(N, 2))

    ST_points = read_coords("coord3D.txt")
    print(ST_points)

    # flow
    flow = list(np.random.randint(low=1, high=5, size=N))


    grid_shape = (10, 10, 10)
    grid_points, grid_size = generate_grid_points(grid_shape)

    # P value
    P = 4

    supernodesInVertex, node_to_supernode = solve_P_Median(ST_points, grid_points, flow, P)


    draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode)

if __name__ == "__main__":
    main()
