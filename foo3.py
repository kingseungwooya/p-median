from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from dotenv import load_dotenv
import os
import random
from matplotlib.lines import Line2D

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

# draw_graph(src_points, sink_points, src_flow, sink_flow, grid_points, supernodesInVertex, src_to_supernode, sink_to_supernode)
def draw_graph(src_points, sink_points, src_flow, sink_flow, grid_points, supernodesInVertex, src_to_supernode, sink_to_supernode):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    for idx, (point, f) in enumerate(zip(src_points, src_flow)):
        ax.scatter(point[0], point[1], point[2], color='red')
        ax.text(point[0], point[1], point[2], f'index={idx} : of={f}', color='red')

    for idx, (point, f) in enumerate(zip(sink_points, sink_flow)):
        ax.scatter(point[0], point[1], point[2], color='blue')
        ax.text(point[0], point[1], point[2], f'index={idx} : if={f}', color='blue')
    # Create custom legend
    custom_lines = [Line2D([0], [0], color='red', marker='o', linestyle='None', markersize=10, label='Source Node'),
                    Line2D([0], [0], color='blue', marker='o', linestyle='None', markersize=10, label='Sink Node')]
    # Supernodes 그리기
    for i in range(0, len(supernodesInVertex) - 1, 2):
        l1, j1 = supernodesInVertex[i]
        l2, j2 = supernodesInVertex[i + 1]

        print("super node info")
        print("l1 = ", l1, ", j1 = ", j1)
        print("l2 = ", l2, ", j2 = ", j2)
        color = random.random(), random.random(), random.random()
        ax.scatter(grid_points[j1, 0], grid_points[j1, 1], 9, color=color, marker='o', s=150)
        ax.scatter(grid_points[j2, 0], grid_points[j2, 1], 9, color=color, marker='o', s=150)

        # Super node 잇기
        ax.plot([grid_points[j1, 0], grid_points[j2, 0]],
                [grid_points[j1, 1], grid_points[j2, 1]],
                [9, 9],
                color=color, lw=5)

    # Node_to_supernode 연결선 그리기
    for n in src_to_supernode:
        l, i, j = n
        # ax.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]],
        #         [ST_points[i, 2], grid_points[j, 2]], color='black')
        ax.plot([src_points[i, 0], grid_points[j, 0]], [src_points[i, 1], grid_points[j, 1]],
                [src_points[i, 2], 9], color='black')
    for n in sink_to_supernode:
        l, i, j = n
        # ax.plot([ST_points[i, 0], grid_points[j, 0]], [ST_points[i, 1], grid_points[j, 1]],
        #         [ST_points[i, 2], grid_points[j, 2]], color='black')
        ax.plot([sink_points[i, 0], grid_points[j, 0]], [sink_points[i, 1], grid_points[j, 1]],
                [sink_points[i, 2], 9], color='black')

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

class CostCallback:
    def __init__(self):
        self.iterations = []
        self.costs = []

    def __call__(self, model, where):
        if where == GRB.Callback.MIP:
            iteration = model.cbGet(GRB.Callback.MIP_NODCNT)
            cost = model.cbGet(GRB.Callback.MIP_OBJBST)
            if iteration % 1 == 0:  # 예를 들어 10회마다 기록
                self.iterations.append(iteration)
                self.costs.append(cost)

def solve_P_Median(src_points, sink_points, grid_points, src_flow, sink_flow, super_node_count, Ft):
    src_index = range(len(src_points))
    sink_index = range(len(sink_points))
    vertex_index = range(len(grid_points))
    super_node_index = range(super_node_count)

    # position22 = [(i, j) for i in ST_index for j in vertex_index]

    src_distance_matrix = generate_distance_matrix(src_points, grid_points)
    sink_distance_matrix = generate_distance_matrix(sink_points, grid_points)


    model = Model('P-Median')
    decision_var_z_src = [
        (l, i, j) for l in super_node_index for i in src_index for j in vertex_index
    ]
    decision_var_z_sink = [
        (l, i, j) for l in super_node_index for i in sink_index for j in vertex_index
    ]
    decision_var_y = [
        (l, j) for l in super_node_index for j in vertex_index
    ]

    z_src = model.addVars(decision_var_z_src, vtype=GRB.BINARY, name='Z_1')
    z_sink = model.addVars(decision_var_z_sink, vtype=GRB.BINARY, name='Z_2')
    y = model.addVars(decision_var_y, vtype=GRB.BINARY, name='Y')

    # 어짜피 src와 sink의 index는 동일하게 움직여야한다. 고로 같은 index를 써도 상관없지 않을까? 그대신 j는 달라져야함 하지만 각 j의 index는 같다.
    model.setObjective(
        quicksum(
            src_distance_matrix[i, j] * z_src[l, i, j] * abs(src_flow[i]) for l, i, j in decision_var_z_src
        )
        + quicksum(
            sink_distance_matrix[i, j] * z_sink[l, i, j] * abs(sink_flow[i]) for l, i, j in decision_var_z_sink
        )
        , GRB.MINIMIZE
    )

    # s.t 10
    model.addConstrs(
        quicksum(z_src[l, i, j] for l in super_node_index for j in vertex_index) == 1
        # 모든 i에 대해..
        for i in src_index
    )
    # s.t 10-1
    model.addConstrs(
        quicksum(z_sink[l, i, j] for l in super_node_index for j in vertex_index) == 1
        # 모든 i에 대해..
        for i in sink_index
    )
    # s.t 11
    model.addConstrs(
        z_src[l, i, j] - y[l, j] <= 0
        for l in super_node_index for i in src_index for j in vertex_index
    )
    # s.t 11-1
    model.addConstrs(
        z_sink[l, i, j] - y[l, j] <= 0
        for l in super_node_index for i in sink_index for j in vertex_index
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
        quicksum(y[l, j] for l in super_node_index for j in vertex_index) == super_node_count
    )
    # 15 초노드에 도달하는 총 flow가 해당 용량 Ft보다 작아야한다.
    model.addConstrs(
        quicksum(abs(src_flow[i]) * z_src[l, i, j] for i in src_index) <= Ft * y[l, j]
        # 모든 j에 대해.
        for l in super_node_index for j in vertex_index
    )
    # 15-1
    model.addConstrs(
        quicksum(abs(sink_flow[i]) * z_sink[l, i, j] for i in sink_index) <= Ft * y[l, j]
        # 모든 j에 대해.
        for l in super_node_index for j in vertex_index
    )

    # (16)
    model.addConstrs(
        (quicksum(sink_flow[i] * z_sink[l, i, j] for j in vertex_index) ==
         -quicksum(src_flow[i] * z_src[l + 1, i, j] for j in vertex_index))
         for l in range(0, len(super_node_index) - 1, 2) for i in sink_index
    )
    # j는 따로여도 같은 i에 대해서 각 z와 z+1번째에 매칭되어야함
    model.addConstrs(
        (quicksum(z_src[l, i, j] for j in vertex_index)) == quicksum(z_sink[l+1, i, j] for j in vertex_index)
         for l in range(0, len(super_node_index) - 1, 2) for i in sink_index
    )


    # # (17) and z constraint
    # for l in super_node_index:
    #     for j in vertex_index:
    #         model.addConstr(y[l, j] <= 1, name=f'c17_y_{l}_{j}')
    # for l in super_node_index:
    #     for i in ST_index:
    #         for j in vertex_index:
    #             model.addConstr(z[l, i, j] <= 1, name=f'c17_z_{l}_{i}_{j}')
    # 제약 조건 추가
    # 만약 j번쨰 vertexindex에 flow가 1 3 5들어온다면 ?
    # 만약 2번째 슈퍼노드이며, 2번째 STindex가 j가 connected 되어있다면 1 그러면서 flow가
    # 이터레이션을 l -> j
    # for j in vertex_index:
    #     model.addConstr(
    #         quicksum(abs(flow[i]) * z[l, i, j] for i in ST_index for l in super_node_index) <=
    #         quicksum(Ft * y[l, j] for l in super_node_index)
    #     )

    callback = CostCallback()
    model.optimize(callback)
    # 최적화 과정의 Cost 그래프 그리기
    plt.plot(callback.iterations, callback.costs, marker='o')
    plt.xlabel('Iteration')
    plt.ylabel('Cost')
    plt.title('Optimization Cost Over Iterations')
    plt.grid(True)
    plt.show()
    supernodesInVertex = [k for k in decision_var_y if y[k].x > 0.5]
    src_to_supernode = [k for k in decision_var_z_src if z_src[k].x > 0.5]
    sink_to_supernode = [k for k in decision_var_z_sink if z_sink[k].x > 0.5]

    print(supernodesInVertex)
    print(src_to_supernode)
    print(sink_to_supernode)


    return supernodesInVertex, src_to_supernode, sink_to_supernode


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

    src_points = read_coords("coord_src.txt")
    sink_points = read_coords("coord_sink.txt")
    N = len(src_points) + len(sink_points)

    src_flow, sink_flow = generate_flows(N)
    flows = src_flow + sink_flow
    print("flow:", flows)
    cetaT = sum(src_flow)
    print(cetaT)
    # maximum number of services of type t that can be jointly routed in the same cable tray
    Ft = 15
    q = int(cetaT / Ft) + 1
    print(q)
    L_length = 2 * q
    print(L_length)

    grid_shape = (10, 10, 10)
    grid_points, grid_size = generate_grid_points(grid_shape)

    supernodesInVertex, src_to_supernode, sink_to_supernode = solve_P_Median(src_points, sink_points, grid_points, src_flow, sink_flow, L_length, Ft)

    draw_graph(src_points, sink_points, src_flow, sink_flow, grid_points, supernodesInVertex, src_to_supernode, sink_to_supernode)


if __name__ == "__main__":
    main()
