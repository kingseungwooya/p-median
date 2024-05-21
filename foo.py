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

    # Super Node 끼리 연결선 그리기
    # print("supernode index들")
    # print(supernodesInVertex)
    # for i in range(len(supernodesInVertex)):
    #
    #     try:
    #         index1 = i
    #         index2 = -i
    #         src_sn_point = grid_points[supernodesInVertex[i]]
    #         sink_sn_point = grid_points[supernodesInVertex[i + 2]]
    #         print(src_sn_point[0])
    #
    #         ax.plot([src_sn_point[0], sink_sn_point[0]], [src_sn_point[1], sink_sn_point[1]],
    #                 [9, 9], color='green', lw=5)
    #     except IndexError:
    #         continue

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


def solve_P_Median(ST_points, grid_points, flow, P, Ft):
    ST_index = range(len(ST_points))
    vertex_index = range(len(grid_points))
    position = [(i, j) for i in ST_index for j in vertex_index]

    distance_matrix = generate_distance_matrix(ST_points, grid_points)

    model = Model('P-Median')

    z = model.addVars(position, vtype=GRB.BINARY, name='X')
    y = model.addVars(vertex_index, vtype=GRB.BINARY, name='Y')

    model.setObjective(
        quicksum(distance_matrix[i, j] * z[i, j] * abs(flow[i]) for i, j in position)
        , GRB.MINIMIZE
    )

    # s.t 10
    model.addConstrs(
        quicksum(z[i, j] for j in vertex_index) == 1 for i in ST_index
    )
    # s.t 11
    model.addConstrs(
        z[i, j, l] - y[j] <= 0 for i in ST_index for j in vertex_index
    )

    # s.t 14
    model.addConstr(
        quicksum(y[j] for j in vertex_index) <= P
    )
    # # 15 초노드에 도달하는 총 flow가 해당 용량 Ft보다 작아야한다.
    model.addConstrs(
        quicksum(abs(flow[i]) * z[i, j] for i in ST_index) <= Ft * y[j] for j in vertex_index
    )
    # #  16 홀수 번째 z들의 합과 짝수 번째 z들의 합이 같도록 제약 조건 추가
    # odd_sum = quicksum(z[i, j] * flow[i] for idx, (i, j) in enumerate(position) if idx % 2 == 0)
    # even_sum = quicksum(z[i, j] * flow[i] for idx, (i, j) in enumerate(position) if idx % 2 != 0)
    # model.addConstr(odd_sum == -even_sum)

    # L이 홀수일 때 L번째 z[i, j]는 L+1번째 z[i, j]와 같다
    # for idx in range(1, len(position), 2):  # idx가 1, 3, 5, ... 인 경우
    #     i1, j1 = position[idx - 1]  # L번째 (idx가 1-based이므로 실제로는 idx-1)
    #     i2, j2 = position[idx]  # L+1번째
    #     print("----------------------" )
    #     print(idx)
    #     print(quicksum(z[i1, j1] * flow[i] for i in ST_index))
    #     even = -quicksum(z[i2, j2] * flow[i] for i in ST_index)
    #     print(even)
    #     model.addConstr(quicksum(z[i1, j1] * flow[i] for i in ST_index) == -quicksum(z[i2, j2] * flow[i] for i in ST_index))

    model.optimize()

    # z[k].x == 1인 것들 중에서 홀수 번째 인덱스의 z와 그 다음 인덱스의 z가 같다는 제약 조건 추가
    print()

    # active_positions = [(i, j) for i in ST_index for j in vertex_index if z[i, j].x > 0.5]
    # print('---------------------')
    # print(active_positions)
    # for idx in range(0, len(active_positions) - 1, 2):  # idx가 0, 2, 4, ... 인 경우
    #     i1, j1 = active_positions[idx]  # 홀수 번째
    #     i2, j2 = active_positions[idx + 1]  # 그 다음 인덱스
    #     model.addConstr(quicksum(z[i1, j1] * flow[i] for i in ST_index) == -(quicksum(z[i2, j2] * flow[i] for i in ST_index)))


    print('---------------------')
    supernodesInVertex = [k for k in vertex_index if y[k].x == 1]
    node_to_supernode = [k for k in position if z[k].x == 1]

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
    print('어떻게 바뀔까?')
    supernodesInVertex = [k for k in vertex_index if y[k].x > 0.5]
    node_to_supernode = [k for k in position if z[k].x > 0.5]

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

    # P value
    P = 4
    #
    supernodesInVertex, node_to_supernode = solve_P_Median(ST_points, grid_points, flows, L_length, Ft)
    #
    #
    draw_graph(ST_points, flows, grid_points, supernodesInVertex, node_to_supernode)


if __name__ == "__main__":
    main()

# flow: [3, 4, 4, 2, 4, 4, 4, 3, 4, 1, -3, -4, -4, -2, -4, -4, -4, -3, -4, -1]


# ---------------------
# [16, 45, 60, 93, 934, 995]
# [(0, 93), (1, 16), (2, 60), (3, 60), (4, 16), (5, 45), (6, 60), (7, 45), (8, 93), (9, 16), (10, 995), (11, 995), (12, 934), (13, 934), (14, 934), (15, 995), (16, 995), (17, 995), (18, 934), (19, 934)]
# 어떻게 바뀔까?

# [15, 51, 55, 82, 924, 985]
# [(0, 82), (1, 15), (2, 51), (3, 51), (4, 15), (5, 55), (6, 51), (7, 55), (8, 82), (9, 15), (10, 985), (11, 985), (12, 924), (13, 985), (14, 924), (15, 985), (16, 985), (17, 985), (18, 924), (19, 924)]