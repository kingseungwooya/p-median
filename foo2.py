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
            distance = np.sqrt((ST_points[i][0] - grid_points[j][0]) ** 2 +
                               (ST_points[i][1] - grid_points[j][1]) ** 2 +
                               (ST_points[i][2] - grid_points[j][2]) ** 2)
            distance_matrix[(i, j)] = distance
    return distance_matrix


def draw_graph(ST_points, flow, grid_points, supernodesInVertex, node_to_supernode):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    for point, f in zip(ST_points, flow):
        # ax.scatter(point[0], point[1], point[2], color='blue')
        ax.text(point[0], point[1], point[2], f'd={f}', color='red')
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


def solve_P_Median(ST_points, grid_points, flow, super_node_count, Ft):
    ST_index = range(len(ST_points))
    vertex_index = range(len(grid_points))
    super_node_index = range(super_node_count)

    position22 = [(i, j) for i in ST_index for j in vertex_index]

    distance_matrix = generate_distance_matrix(ST_points, grid_points)

    model = Model('P-Median')
    decision_var_z = [
        (l, i, j) for l in super_node_index for i in ST_index for j in vertex_index
    ]
    decision_var_y = [
        (l, j) for l in super_node_index for j in vertex_index
    ]

    z = model.addVars(decision_var_z, vtype=GRB.BINARY, name='X')
    y = model.addVars(decision_var_y, vtype=GRB.BINARY, name='Y')

    model.setObjective(
        quicksum(distance_matrix[i, j] * z[l, i, j] * abs(flow[i]) for l, i, j in decision_var_z)
        , GRB.MINIMIZE
    )

    # s.t 10
    model.addConstrs(
        quicksum(z[l, i, j] for l in super_node_index for j in vertex_index) == 1
        # 모든 i에 대해..
        for i in ST_index
    )
    # s.t 11
    model.addConstrs(
        z[l, i, j] - y[l, j] <= 0
        for l in super_node_index for i in ST_index for j in vertex_index
    )
    # s.t 12
    model.addConstrs(
        quicksum(y[l, j] for l in super_node_index) <= 1
        # 모든 j에 대해..
        for j in vertex_index
    )
    # s.t 13
    model.addConstrs(
        quicksum(y[l, j] for j in vertex_index) == 1
        # 모든 l에 대해..
        for l in super_node_index
    )
    # s.t 14
    model.addConstr(
        quicksum(y[l, j] for l in super_node_index for j in vertex_index) <= super_node_count
    )
    # 15 초노드에 도달하는 총 flow가 해당 용량 Ft보다 작아야한다.
    model.addConstrs(
        quicksum(abs(flow[i]) * z[l, i, j] for l in super_node_index for i in ST_index) <= Ft * y[l2, j] for l2 in super_node_index
        # 모든 j에 대해.
        for j in vertex_index
    )

    model.optimize()
    if model.status == GRB.OPTIMAL:
        supernodesInVertex = [k for k in decision_var_y if y[k].x > 0.5]
        node_to_supernode = [k for k in decision_var_z if z[k].x > 0.5]

        print("Supernodes in Vertex:")
        print(supernodesInVertex)
        print("Node to Supernode mapping:")
        print(node_to_supernode)

        # supernodesInVertex와 node_to_supernode를 출력하여 확인
        for l, j in supernodesInVertex:
            print(f'Supernode {l} at Vertex {j}')

        for l, i, j in node_to_supernode:
            print(f'Node {i} assigned to Supernode {l} at Vertex {j}')
    else:
        print("Optimal solution was not found.")

    supernodesInVertex = [k for k in decision_var_y if y[k].x > 0.5]
    node_to_supernode = [k for k in decision_var_z if z[k].x > 0.5]

    print(supernodesInVertex)
    print(node_to_supernode)

    for idx in range(0, len(supernodesInVertex) - 1, 2):  # idx가 0, 2, 4, ... 인 경우
        j1 = supernodesInVertex[idx]  # 홀수 번째
        print("j1")
        print(j1)
        print(quicksum(z[i, j1] * flow[i] for i in ST_index))
        j2 = supernodesInVertex[idx + 1]  # 그 다음 인덱스
        print("j2")
        print(j2)
        print( (quicksum(z[i, j2] * -flow[i] for i in ST_index)))
        model.addConstr(
            quicksum(z[i, j1] * flow[i] for i in ST_index) == (quicksum(z[i, j2] * -flow[i] for i in ST_index))
        )

    model.optimize()
    supernodesInVertex = [k for k in decision_var_y if y[k].x > 0.5]
    node_to_supernode = [k for k in decision_var_z if z[k].x > 0.5]

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


def generate_flows(N):
    np.random.seed(1004)
    src_flow = list(np.random.randint(low=1, high=5, size=int(N / 2)))
    sink_flow = [-x for x in src_flow]

    return src_flow, sink_flow


def main():
    # ST_points = np.random.randint(low=0, high=10, size=(N, 2))

    ST_points = read_coords("coord3D.txt")
    N = len(ST_points)

    src_flow, sink_flow = generate_flows(N)
    flows = src_flow + sink_flow
    print("flow:", flows)
    cetaT = sum(src_flow)
    print(cetaT)
    # maximum number of services of type t that can be jointly routed in the same cable tray
    Ft = 10
    q = int(cetaT / Ft)
    print(q)
    L_length = 2 * q
    print(L_length)

    grid_shape = (10, 10, 10)
    grid_points, grid_size = generate_grid_points(grid_shape)


    supernodesInVertex, node_to_supernode = solve_P_Median(ST_points, grid_points, flows, L_length, Ft)
    #
    #
    draw_graph(ST_points, flows, grid_points, supernodesInVertex, node_to_supernode)


if __name__ == "__main__":
    main()

